%% Master Script: Multicellular and Unicellular Evolution
clearvars
clc
tic

% ---- User-defined Settings ----
runNo = 1;
fitnessMethod = 'geometric';  % Options: 'arithmetic' (default), 'geometric'
burnInSteps = 200;
% --------------------------------

outputFolder = sprintf('simData_runNo_%d_%s', runNo, fitnessMethod);
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

L_values = [100];
K_values = [4,8,16];
delta_values = [0.15,0.75,1.25,1.75];
selectionPressureSets = {[0.5, 0.5],[0.01, 0.99]};
numSteps = 400;
numSimulations = 50;

precomputedGenomes = generateInitialGenomes(numSimulations, max(K_values), L_values(1));
referenceGenome = precomputedGenomes{1};
environmentStore = prepareEnvironments(L_values(1), delta_values, 1.25);

paramCombinations = {};
for L = L_values
    for K = K_values
        for delta = delta_values
            for selIdx = 1:length(selectionPressureSets)
                selPress = selectionPressureSets{selIdx};
                paramCombinations{end+1, 1} = L;
                paramCombinations{end,   2} = K;
                paramCombinations{end,   3} = delta;
                paramCombinations{end,   4} = selPress;
            end
        end
    end
end

for mode = ["multicell"]
    for idx = 1:size(paramCombinations, 1)
        L = paramCombinations{idx,1};
        K = paramCombinations{idx,2};
        delta = paramCombinations{idx,3};
        selPress = paramCombinations{idx,4};

        runSimulations(mode, L, K, delta, selPress, ...
            numSteps, numSimulations, ...
            precomputedGenomes, referenceGenome, environmentStore, outputFolder, fitnessMethod, burnInSteps);
    end
end


%% Helper Functions

function initialGenomes = generateInitialGenomes(numSim, maxK, L)
    initialGenomes = cell(numSim, 1);
    for i = 1:numSim
        rng(1000 + i);
        initialGenomes{i} = double(rand(L, maxK) < (1 / maxK));
    end
end

function envMap = prepareEnvironments(L, deltaVals, refDelta)
    refSeed = 200 + round(refDelta * 100);
    E_ref = createEnv(L, refDelta, refSeed);
    envMap = containers.Map('KeyType','double','ValueType','any');
    envMap(refDelta) = E_ref;
    for d = setdiff(deltaVals, refDelta)
        envMap(d) = exaggerateEnv(E_ref, d);
    end
end


function runSimulations(mode, L, K, delta, selPress, ...
                        numSteps, numSim, ...
                        precomputedGenomes, referenceGenome, ...
                        envStore, outFolder, fitnessMethod, burnInSteps)

    if ~envStore.isKey(delta)
        error('Missing environment for delta=%.2f', delta);
    end
    targetEnv = envStore(delta);
    deltaE = norm(targetEnv(:,1) - targetEnv(:,2));
    numEnv = size(targetEnv, 2);

    initialGenomes = cell(numSim, 1);
    for i = 1:numSim
        initialGenomes{i} = precomputedGenomes{i}(:, 1:K);
    end

    if ~isequal(initialGenomes{1}, referenceGenome(:, 1:K))
        error('Initial genome mismatch.');
    end

    allMeanFitness = NaN(numSteps, numSim);
    allFitnessHistory = cell(numSim, 1);
    allGenomes = cell(numSim, 1);
    allTradeoffIndices = NaN(numSteps, numSim);
    allModularityIndices = NaN(numSteps, numSim);
    allOptimalCoefficients = cell(numSim, 1);
    allGmats = cell(numSteps, numSim);
    allCondEvol = cell(numSteps, numSim);
    allAutonomy = cell(numSteps, numSim);

    parfor i = 1:numSim
        simResult = [];
        genome = initialGenomes{i};

        if strcmp(mode, 'multicell')
            simResult = simulateMulticellEvolution(genome, L, K, numSteps, numEnv, targetEnv, selPress, fitnessMethod, burnInSteps);
            allOptimalCoefficients{i} = simResult.optimalCoefficients;
            allGmats(:, i) = simResult.Gmut;
            allCondEvol(:, i) = simResult.conditionalEvolvability;
            allAutonomy(:, i) = simResult.autonomy;
        elseif strcmp(mode, 'unicell')
            simResult = simulateUnicellEvolution(genome, L, K, numSteps, numEnv, targetEnv);
            allOptimalCoefficients{i} = [];
        else
            error('Unknown mode: %s', mode);
        end

        allMeanFitness(:, i) = simResult.totalFitness;
        allFitnessHistory{i} = simResult.tissueFitness;
        allGenomes{i} = simResult.genomeHistory;
        allTradeoffIndices(:, i) = simResult.tradeoffIndex;
        allModularityIndices(:, i) = simResult.modularityIndex;
    end

    params.numSteps = numSteps;
    params.L = L;
    params.K = K;
    params.delta = delta;
    params.selectionPressure = selPress;
    params.targetEnvironments = targetEnv;
    params.deltaE = deltaE;
    params.numEnvironments = numEnv;
    params.numSimulations = numSim;
    params.fitnessMethod = fitnessMethod;
    params.burnInSteps = burnInSteps;

    results.allMeanFitness = allMeanFitness;
    results.allFitnessHistory = allFitnessHistory;
    results.allGenomes = allGenomes;
    results.allTradeoffIndices = allTradeoffIndices;
    results.allModularityIndices = allModularityIndices;
    results.allOptimalCoefficients = allOptimalCoefficients;
    results.initialGenomes = initialGenomes;
    results.allGmats = allGmats;
    results.allCondEvol = allCondEvol;
    results.allAutonomy = allAutonomy;

    selStr = sprintf('%.2f-', selPress); selStr = selStr(1:end-1);
    fname = sprintf('%s_L%d_K%d_deltaE%.2f_selPress%s.mat', ...
        capitalize(mode), L, K, deltaE, selStr);

    save(fullfile(outFolder, fname), 'params', 'results');
    fprintf('Saved: %s\n', fullfile(outFolder, fname));
end

function E_new = exaggerateEnv(E_ref, newDelta)
    envA = E_ref(:,1); envB = E_ref(:,2);
    logA = log(envA); logB = log(envB);
    logMean = 0.5 * (logA + logB);
    logA_prime = logMean + (newDelta/2) * (logA - logB);
    logB_prime = logMean - (newDelta/2) * (logA - logB);
    envA_prime = exp(logA_prime) / norm(exp(logA_prime));
    envB_prime = exp(logB_prime) / norm(exp(logB_prime));
    E_new = [envA_prime, envB_prime];
end

function str = capitalize(s)
    str = [upper(s(1)), s(2:end)];
end
