clearvars
clc

replicateThreshold = 10;

rootFolder = 'simData_runNo_1_geometric_burnin0';
figFolder = 'fig_simData_runNo_1_geometric_burnin0';
if ~exist(figFolder, 'dir')
    mkdir(figFolder);
end

matFiles = dir(fullfile(rootFolder, '*L*_K*_deltaE*_selPress*.mat'));

selLabels = {};
deltaLabels = {};
for k = 1:length(matFiles)
    selToken = regexp(matFiles(k).name, 'selPress([\d\.\-_]+)', 'tokens');
    deltaToken = regexp(matFiles(k).name, 'deltaE([\d\.]+)', 'tokens');
    if ~isempty(selToken)
        selLabels{end+1} = strrep(selToken{1}{1}, '-', '_');
    end
    if ~isempty(deltaToken)
        deltaLabels{end+1} = deltaToken{1}{1};
    else
        deltaLabels{end+1} = 'NA';
    end
end
uniqueSelLabels = unique(selLabels);

for selIdx = 1:length(uniqueSelLabels)
    selLabel = uniqueSelLabels{selIdx};
    selectedFiles = matFiles(contains({matFiles.name}, ['selPress' strrep(selLabel, '_', '-')]));
    deltaForFiles = deltaLabels(contains(selLabels, selLabel));
    numDeltas = length(selectedFiles);

    fig = figure('Position', [100, 100, 2200, 350 * numDeltas]);

    for row = 1:numDeltas
        filename = fullfile(selectedFiles(row).folder, selectedFiles(row).name);
        load(filename, 'results');

        numSteps = size(results.allCondEvol, 1);
        numReps = size(results.allCondEvol, 2);
        x = (1:numSteps)';

        condMat1 = nan(numSteps, numReps);
        condMat2 = nan(numSteps, numReps);
        autonomyMat = nan(numSteps, numReps);
        for i = 1:numSteps
            for j = 1:numReps
                val = results.allCondEvol{i, j};
                if ~isempty(val) && all(~isnan(val))
                    condMat1(i,j) = val(1);
                    condMat2(i,j) = val(2);
                end
                val2 = results.allAutonomy{i, j};
                if ~isempty(val2) && all(~isnan(val2))
                    autonomyMat(i,j) = val2(1);
                end
            end
        end

        rowBase = (row - 1) * 4;

        subplot(numDeltas, 4, rowBase + 1); hold on;
        for j = 1:numReps
            semilogy(x, condMat1(:,j), 'Color', [0 0 1 0.2]);
            semilogy(x, condMat2(:,j), 'Color', [1 0 0 0.2]);
        end
        count1 = sum(~isnan(condMat1), 2);
        count2 = sum(~isnan(condMat2), 2);
        mean1 = nanmean(condMat1, 2);
        mean2 = nanmean(condMat2, 2);
        mean1(count1 < replicateThreshold) = NaN;
        mean2(count2 < replicateThreshold) = NaN;
        h1 = semilogy(x, mean1, ':', 'Color', [0 0 1], 'LineWidth', 2);
        h2 = semilogy(x, mean2, ':', 'Color', [1 0 0], 'LineWidth', 2);
        ylabel(['\DeltaE = ', deltaForFiles{row}]);
        if row == numDeltas, xlabel('Mutation Steps'); end
        legend([h1 h2], {'Trait 1', 'Trait 2'}, ...
            'Orientation', 'horizontal', ...
            'Location', 'northeast', ...
            'Box', 'off');
        title('Conditional Evolvability');
        set(gca, 'YScale', 'log');
        ylim([1e-7, 1e-1]);

        subplot(numDeltas, 4, rowBase + 2); hold on;
        for j = 1:numReps
            plot(x, autonomyMat(:,j), 'Color', [0 0 0 0.2]);
        end
        countA = sum(~isnan(autonomyMat), 2);
        meanA = nanmean(autonomyMat, 2);
        meanA(countA < replicateThreshold) = NaN;
        plot(x, meanA, ':k', 'LineWidth', 2);
        if row == numDeltas, xlabel('Mutation Steps'); end
        title('Trait Autonomy');
        ylim([0,1]);

        subplot(numDeltas, 4, rowBase + 3); hold on;
        if iscell(results.allModularityIndices)
            M = cellfun(@(c) c(min(end,1)), results.allModularityIndices);
        else
            M = results.allModularityIndices;
        end
        for j = 1:numReps
            plot(x, M(:,j), 'Color', [0 0 0 0.2]);
        end
        countM = sum(~isnan(M), 2);
        meanM = nanmean(M, 2);
        meanM(countM < replicateThreshold) = NaN;
        plot(x, meanM, ':k', 'LineWidth', 2);
        if row == numDeltas, xlabel('Mutation Steps'); end
        title('Mutational Modularity');
        ylim([0,1]);

        subplot(numDeltas, 4, rowBase + 4); hold on;
        if iscell(results.allTradeoffIndices)
            T = cellfun(@(c) c(min(end,1)), results.allTradeoffIndices);
        else
            T = results.allTradeoffIndices;
        end
        for j = 1:numReps
            plot(x, T(:,j), 'Color', [0 0 0 0.2]);
        end
        countT = sum(~isnan(T), 2);
        meanT = nanmean(T, 2);
        meanT(countT < replicateThreshold) = NaN;
        plot(x, meanT, ':k', 'LineWidth', 2);
        if row == numDeltas, xlabel('Mutation Steps'); end
        title('Mutational Tradeoff');
        ylim([-1,1]);
    end

    saveas(fig, fullfile(figFolder, sprintf('Multicell_conditional_%s.png', selLabel)));
    close(fig);
end
