function targetEnvironments = createEnv(L, delta, rngSeed)
    % createEnv generates two environment vectors with tunable divergence in log space
    % based on Tikhonov et al.'s model description.
    %
    % Inputs:
    %   L - Dimensionality of the phenotype space (integer)
    %   delta - Divergence parameter (0 results in identical environments,
    %           larger values increase divergence)
    %   rngSeed - Random seed for reproducibility (integer)
    %
    % Output:
    %   targetEnvironments - A matrix containing two environment vectors 
    %                        with controlled divergence, size [L, 2]

    rng(rngSeed);  % Set the random seed for reproducibility
    
    % Step 1: Generate two initial random environment vectors
    envA = abs(rand(L, 1));  % Generate a random vector with positive components
    envB = abs(rand(L, 1));  % Another random vector with positive components

    % Step 2: Calculate the geometric mean of envA and envB in log space
    logEnvA = log(envA);
    logEnvB = log(envB);
    logEnvMean = 0.5 * (logEnvA + logEnvB);

    % Step 3: Apply the logarithmic transformation to parameterize divergence
    logEnvA_prime = logEnvMean + (delta / 2) * (logEnvA - logEnvB);
    logEnvB_prime = logEnvMean - (delta / 2) * (logEnvA - logEnvB);

    % Step 4: Convert back from log space to get the parameterized environment vectors
    envA_prime = exp(logEnvA_prime);
    envB_prime = exp(logEnvB_prime);

    % Step 5: Normalize the resulting environment vectors to unit length
    envA_prime = envA_prime / norm(envA_prime);
    envB_prime = envB_prime / norm(envB_prime);

    % Store the environment vectors in a matrix
    targetEnvironments = [envA_prime, envB_prime];
end
