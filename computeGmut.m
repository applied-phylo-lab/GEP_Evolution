function [Gmut, condEvolvability, autonomy] = computeGmut(mutantFitnesses, wtFitness)
    % computeGmut - Compute mutational G-matrix, conditional evolvability, and autonomy
    %
    % Inputs:
    %   mutantFitnesses - Matrix (N x 2) of fitnesses for mutants across two environments
    %   wtFitness       - 1 x 2 vector of wildtype fitnesses
    %
    % Outputs:
    %   Gmut               - 2x2 mutational G-matrix
    %   condEvolvability   - Conditional evolvability of trait 1 given 2 and vice versa
    %   autonomy           - Autonomy of trait 1 and 2

    if size(mutantFitnesses, 2) ~= 2 || length(wtFitness) ~= 2
        error('This function assumes exactly 2 environments.');
    end

    deltas = mutantFitnesses - wtFitness;  % N x 2
    Gmut = cov(deltas);  % 2 x 2 matrix

    % Conditional evolvability
    e1_given_2 = Gmut(1,1) - (Gmut(1,2)^2 / Gmut(2,2));
    e2_given_1 = Gmut(2,2) - (Gmut(1,2)^2 / Gmut(1,1));

    % Autonomy
    alpha1 = e1_given_2 / Gmut(1,1);
    alpha2 = e2_given_1 / Gmut(2,2);

    condEvolvability = [e1_given_2, e2_given_1];
    autonomy = [alpha1, alpha2];
end
