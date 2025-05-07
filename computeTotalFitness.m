
function totalFitness = computeTotalFitness(fitnessMatrix, weights, method)
    % computeTotalFitness - Calculate the total fitness using either arithmetic or geometric mean
    %
    % Syntax: totalFitness = computeTotalFitness(fitnessMatrix, weights, method)
    %
    % Inputs:
    %   fitnessMatrix - Matrix (N x numEnvironments) of fitness values
    %   weights       - Row vector of weights (1 x numEnvironments)
    %   method        - String, either 'arithmetic' or 'geometric' (default: 'arithmetic')
    %
    % Output:
    %   totalFitness - Column vector (N x 1) of total fitness values

    if nargin < 3
        method = 'arithmetic';
    end

    [N, numEnvironments] = size(fitnessMatrix);
    if length(weights) ~= numEnvironments
        error('Length of weights must match the number of environments.');
    end

    if abs(sum(weights) - 1) > 1e-6
        error('Weights must sum to 1.');
    end

    switch lower(method)
        case 'arithmetic'
            totalFitness = fitnessMatrix * weights(:);
        case 'geometric'
            if any(fitnessMatrix >= 0, 'all')
                error('Fitness values must be negative for geometric mean.');
            end
            transformed = -fitnessMatrix;  % shift from [-1, 0] to [0, 1]
            logTerms = log(transformed) * weights(:);
            totalFitness = -exp(logTerms);
        otherwise
            error('Unknown fitness aggregation method: %s', method);
    end
end
