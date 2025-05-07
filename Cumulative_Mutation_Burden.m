%% Script: Plot Cumulative Mutation Burden (Fastest vs. Slowest on Same Plot)
clearvars
clc

% ---- Configuration ----
dataFolder = 'gepMetrics_runNo_1_geometric_burnin400';
targetFile = 'gepMetrics_MULTICELL_L100_K8_deltaE1.36_selPress0.50-0.50.mat';
replicateThreshold = 50;
% ------------------------

load(fullfile(dataFolder, targetFile), 'gepMetrics');

numSteps = size(gepMetrics.mutationBurden, 1);
numReps = size(gepMetrics.mutationBurden, 2);

% Infer maximum K
K_max = 0;
for t = 1:numSteps
    for rep = 1:numReps
        mu = gepMetrics.mutationBurden{t, rep};
        if ~isempty(mu)
            K_max = max(K_max, numel(mu));
        end
    end
end

% Accumulate burden
cumulativeAll = NaN(numSteps, K_max, numReps);
for rep = 1:numReps
    validT = find(cellfun(@(c) ~isempty(c), gepMetrics.mutationBurden(:, rep)));
    if isempty(validT), continue; end
    maxT = max(validT);
    cumulative = zeros(maxT, K_max);
    for t = 2:maxT
        mu = gepMetrics.mutationBurden{t, rep};
        if isempty(mu), continue; end
        prev = cumulative(t-1, :);
        cumulative(t, 1:numel(mu)) = prev(1:numel(mu)) + mu;
    end
    cumulativeAll(1:maxT, :, rep) = cumulative;
end

% Rank and extract
fastestAll = NaN(numSteps, numReps);
slowestAll = NaN(numSteps, numReps);
for rep = 1:numReps
    traj = cumulativeAll(:,:,rep);
    if all(isnan(traj), 'all'), continue; end
    lastValidRow = find(any(~isnan(traj),2), 1, 'last');
    if isempty(lastValidRow), continue; end
    finalBurden = traj(lastValidRow, :);
    if all(isnan(finalBurden)), continue; end
    [~, sortedIdx] = sort(finalBurden, 'descend', 'MissingPlacement', 'last');
    fastestIdx = sortedIdx(1);
    slowestIdx = sortedIdx(find(~isnan(finalBurden(sortedIdx)), 1, 'last'));
    fastestAll(:, rep) = traj(:, fastestIdx);
    slowestAll(:, rep) = traj(:, slowestIdx);
end

% Thresholding
x = (1:numSteps)';
validFast = ~isnan(fastestAll);
validSlow = ~isnan(slowestAll);
repCountsFast = sum(validFast, 2);
repCountsSlow = sum(validSlow, 2);
meanFast = nanmean(fastestAll, 2);
meanSlow = nanmean(slowestAll, 2);
meanFast(repCountsFast < replicateThreshold) = NaN;
meanSlow(repCountsSlow < replicateThreshold) = NaN;

% Plot
figure('Position', [100, 100, 800, 400]); hold on;
for rep = 1:numReps
    if ~all(isnan(fastestAll(:, rep)))
        plot(x, fastestAll(:, rep), 'Color', [1 0 0 0.2]);
    end
    if ~all(isnan(slowestAll(:, rep)))
        plot(x, slowestAll(:, rep), 'Color', [0 0 1 0.2]);
    end
end
plot(x, meanFast, 'r:', 'LineWidth', 2);
plot(x, meanSlow, 'b:', 'LineWidth', 2);
xlabel('Mutation Steps');
ylabel('Cumulative Mutation Burden');
title('Fastest vs. Slowest Evolving GEPs');
ylim([0, max([meanFast; meanSlow], [], 'omitnan') * 1.2]);
legend({'Fastest (rep)', 'Slowest (rep)', 'Fastest (mean)', 'Slowest (mean)'}, ...
       'Location', 'northwest', 'Box', 'off');
