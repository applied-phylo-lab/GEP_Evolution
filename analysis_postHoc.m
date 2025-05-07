%% Script: Extract GEP Usage and Structural Metrics
clearvars
clc

dataFolder = 'simData_runNo_1_geometric_burnin400';
saveFolder = 'gepMetrics_runNo_1_geometric_burnin400';
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end

matFiles = dir(fullfile(dataFolder, '*.mat'));

for k = 1:length(matFiles)
    load(fullfile(matFiles(k).folder, matFiles(k).name), 'results', 'params');

    numSteps = params.numSteps;
    numReps = params.numSimulations;
    numEnv = params.numEnvironments;
    L = params.L;

    usageOverlap = cell(numSteps, numReps);
    mutationBurden = cell(numSteps, numReps);
    jaccardMatrix = cell(numSteps, numReps);
    deploymentCorr = cell(numSteps, numReps);

    for rep = 1:numReps
        genomes = results.allGenomes{rep};
        coeffs  = results.allOptimalCoefficients{rep};

        for t = 1:numSteps
            G = genomes{t};    % [L x K]
            A = coeffs{t};     % [K x E]
            if isempty(G) || isempty(A), continue; end
            A = A';            % [E x K]
            K = size(A, 2);
            G = G(:, 1:K);

            % --- Usage Overlap ---
            O = NaN(1, K);
            for i = 1:K
                a = A(:, i);
                numer = (sum(a))^2;
                denom = numEnv * sum(a.^2);
                if denom > 0
                    O(i) = numer / denom;
                end
            end
            usageOverlap{t, rep} = O;

            % --- Mutation Burden ---
            if t > 1
                prevG = genomes{t-1};
                prevG = prevG(:, 1:K);
                mu = sum(abs(G - prevG), 1);  % bit-flip mutation count per GEP
                mutationBurden{t, rep} = mu;
            else
                mutationBurden{t, rep} = zeros(1, K);
            end

            % --- Jaccard Similarity Matrix ---
            Gbin = G > 0;
            J = NaN(K, K);
            for i = 1:K
                for j = i+1:K
                    inter = sum(Gbin(:,i) & Gbin(:,j));
                    union = sum(Gbin(:,i) | Gbin(:,j));
                    if union > 0
                        J(i,j) = inter / union;
                        J(j,i) = J(i,j);
                    end
                end
            end
            jaccardMatrix{t, rep} = J;

            % --- Deployment Correlation Matrix ---
            R = corr(A);  % correlation between GEPs by usage pattern
            deploymentCorr{t, rep} = R;
        end
    end

    gepMetrics.usageOverlap = usageOverlap;
    gepMetrics.mutationBurden = mutationBurden;
    gepMetrics.jaccardMatrix = jaccardMatrix;
    gepMetrics.deploymentCorr = deploymentCorr;

    saveName = fullfile(saveFolder, ['gepMetrics_' matFiles(k).name]);
    save(saveName, 'gepMetrics', '-v7.3');
    fprintf('Saved GEP metrics: %s\n', saveName);
end
