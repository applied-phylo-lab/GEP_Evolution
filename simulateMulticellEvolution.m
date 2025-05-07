function [simResults] = simulateMulticellEvolution(genome, L, K, numSteps, numEnvironments, targetEnvironments, selectionPressure, fitnessMethod, burnInSteps)
    if nargin < 9
        burnInSteps = 0;
    end
    if nargin < 8
        fitnessMethod = 'arithmetic';
    end

    simResults.totalFitness = NaN(numSteps, 1);
    simResults.tissueFitness = NaN(numSteps, numEnvironments);
    simResults.genomeHistory = cell(numSteps, 1);
    simResults.tradeoffIndex = NaN(numSteps, 1);
    simResults.modularityIndex = NaN(numSteps, 1);
    simResults.optimalCoefficients = cell(numSteps, 1);
    simResults.Gmut = cell(numSteps, 1);
    simResults.conditionalEvolvability = cell(numSteps, 1);
    simResults.autonomy = cell(numSteps, 1);

    simResults.genomeHistory{1} = genome;
    [simResults.tissueFitness(1, :), a_opt] = computeTissueFitness(genome, targetEnvironments);
    simResults.optimalCoefficients{1} = a_opt;
    currentFitness = computeTotalFitness(simResults.tissueFitness(1, :), selectionPressure, fitnessMethod);
    simResults.totalFitness(1) = currentFitness;

    % Generate 1-mutant neighbors for the initial genome
    oneMutantNeighbors = cell(L * K, 1);
    mutantFitnesses_initial = zeros(L * K, numEnvironments);
    index = 1;
    for row = 1:L
        for col = 1:K
            mutantGenome = genome;
            mutantGenome(row, col) = 1 - mutantGenome(row, col);
            oneMutantNeighbors{index} = mutantGenome;
            mutantFitnesses_initial(index, :) = computeTissueFitness(mutantGenome, targetEnvironments);
            index = index + 1;
        end
    end

    [Gmut, condEvolvability, autonomy] = computeGmut(mutantFitnesses_initial, simResults.tissueFitness(1, :));
    simResults.Gmut{1} = Gmut;
    simResults.conditionalEvolvability{1} = condEvolvability;
    simResults.autonomy{1} = autonomy;

    timestep = 2;
    while timestep <= numSteps
        isBurnIn = (timestep <= burnInSteps + 1);

        oneMutantNeighbors = cell(L * K, 1);
        mutantFitnesses = zeros(L * K, numEnvironments);
        mutantCoefficients = cell(L * K, 1);
        index = 1;

        for row = 1:L
            for col = 1:K
                mutantGenome = genome;
                mutantGenome(row, col) = 1 - mutantGenome(row, col);
                oneMutantNeighbors{index} = mutantGenome;
                [mutantFitnesses(index, :), mutantCoefficients{index}] = ...
                    computeTissueFitness(mutantGenome, targetEnvironments);
                index = index + 1;
            end
        end

        % Compute and store mutational metrics (G, condEvolvability, autonomy)
        [Gmut, condEvolvability, autonomy] = computeGmut(mutantFitnesses, simResults.tissueFitness(timestep-1, :));
        simResults.Gmut{timestep} = Gmut;
        simResults.conditionalEvolvability{timestep} = condEvolvability;
        simResults.autonomy{timestep} = autonomy;

        % Always compute tradeoff and modularity
        deltaF_A = mutantFitnesses(:, 1) - simResults.tissueFitness(timestep-1, 1);
        deltaF_B = mutantFitnesses(:, 2) - simResults.tissueFitness(timestep-1, 2);
        beneficialMask = deltaF_A > 0 | deltaF_B > 0;
        chiNumer = -sum(deltaF_A(beneficialMask) .* deltaF_B(beneficialMask));
        chiDenom = sqrt(sum(deltaF_A(beneficialMask).^2) * sum(deltaF_B(beneficialMask).^2));
        simResults.tradeoffIndex(timestep) = chiNumer / chiDenom;

        rho = sqrt(deltaF_A.^2 + deltaF_B.^2);
        phi = atan2(deltaF_B, deltaF_A);
        modularityNumer = sum(rho .* abs(sin(2 * phi)));
        modularityDenom = sum(rho);
        simResults.modularityIndex(timestep) = 1 - (modularityNumer / modularityDenom);

        if isBurnIn
            chosenIndex = randi(L * K);  % purely random mutation
        else
            totalFitnesses = computeTotalFitness(mutantFitnesses, selectionPressure, fitnessMethod);
            selectionCoeffs = totalFitnesses - currentFitness;
            fixationProb = 1 - exp(-2 .* selectionCoeffs);
            beneficialIndices = find(fixationProb > 0);
            if isempty(beneficialIndices)
                break;
            end
            beneficialFixationProbs = fixationProb(beneficialIndices);
            weights = beneficialFixationProbs / sum(beneficialFixationProbs);
            chosenIndex = beneficialIndices(randsample(length(weights), 1, true, weights));
        end

        genome = oneMutantNeighbors{chosenIndex};
        simResults.genomeHistory{timestep} = genome;
        simResults.tissueFitness(timestep, :) = mutantFitnesses(chosenIndex, :);
        simResults.optimalCoefficients{timestep} = mutantCoefficients{chosenIndex};
        currentFitness = computeTotalFitness(simResults.tissueFitness(timestep, :), selectionPressure, fitnessMethod);
        simResults.totalFitness(timestep) = currentFitness;

        timestep = timestep + 1;
    end

    if timestep <= numSteps
        simResults.genomeHistory(timestep:end) = {NaN};
        simResults.tissueFitness(timestep:end, :) = NaN;
        simResults.totalFitness(timestep:end) = NaN;
        simResults.tradeoffIndex(timestep:end) = NaN;
        simResults.modularityIndex(timestep:end) = NaN;
        simResults.optimalCoefficients(timestep:end) = {NaN};
        simResults.Gmut(timestep:end) = {NaN};
        simResults.conditionalEvolvability(timestep:end) = {NaN};
        simResults.autonomy(timestep:end) = {NaN};
    end
end
