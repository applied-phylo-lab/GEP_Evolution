function [fitnessMatrix, optimalCoefficients] = computeTissueFitness(genome, targetEnvironments)
    % computeTissueFitness - Calculate the fitness of a genome across multiple environments
    %                        and record the optimal coefficients (a_opt) for GEPs.
    %
    % Syntax: [fitnessMatrix, optimalCoefficients] = computeTissueFitness(genome, targetEnvironments)
    %
    % Inputs:
    %   genome - Binary matrix representing the genome (L x K)
    %   targetEnvironments - Matrix of target environments (L x numEnvironments)
    %
    % Outputs:
    %   fitnessMatrix - Vector of fitness values across environments (1 x numEnvironments)
    %   optimalCoefficients - Optimal coefficients (a_opt) for GEPs (K x numEnvironments)

    numEnvironments = size(targetEnvironments, 2); % Number of environments
    fitnessMatrix = zeros(1, numEnvironments); % Initialize fitness matrix
    optimalCoefficients = zeros(size(genome, 2), numEnvironments); % Initialize a_opt matrix

    % Loop over each environment to compute fitness and a_opt
    for e = 1:numEnvironments
        currentEnv = targetEnvironments(:, e); % Retrieve current environment

        % Optimal expression levels using non-negative least squares
        a_opt = lsqnonneg(genome, currentEnv);

        % Store the optimal coefficients
        optimalCoefficients(:, e) = a_opt;

        % Calculate the phenotype based on the optimal expression levels
        optimalPhenotype = genome * a_opt;

        % Compute fitness as the negative Euclidean distance to the target
        fitnessMatrix(e) = -norm(currentEnv - optimalPhenotype);
    end
end
