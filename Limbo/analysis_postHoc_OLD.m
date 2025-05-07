clearvars
clc

dataFolder = 'simData_runNo_1_geometric_burnin0';
saveFolder = 'structuralMetrics_runNo_1_geometric_burnin0';
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end

matFiles = dir(fullfile(dataFolder, '*.mat'));

for k = 1:length(matFiles)
    load(fullfile(matFiles(k).folder, matFiles(k).name), 'results', 'params');

    numSteps = params.numSteps;
    numReps = params.numSimulations;
    numEnv = params.numEnvironments;

    meanUsageOverlap = NaN(numSteps, numReps);
    meanJaccardSimilarity = NaN(numSteps, numReps);

    for rep = 1:numReps
        genomes = results.allGenomes{rep};                  % cell of [L x K]
        coeffs  = results.allOptimalCoefficients{rep};      % cell of [K x E] (transposed!)

        for t = 1:numSteps
            G = genomes{t};  % [L x K]
            A = coeffs{t};   % [K x E] ? transpose needed

            if isempty(G) || isempty(A)
                continue;
            end

            A = A';  % Now A is [E x K] as required

            K_current = size(A, 2);
            G = G(:, 1:K_current);  % Match in case K < expected

            % --- Usage Overlap (per GEP) ---
            overlap = NaN(1, K_current);
            for i = 1:K_current
                a = A(:, i);
                numer = (sum(a))^2;
                denom = numEnv * sum(a.^2);
                if denom > 0
                    overlap(i) = numer / denom;
                end
            end
            meanUsageOverlap(t, rep) = nanmean(overlap);

            % --- Jaccard Similarity ---
            Gbin = G > 0;
            simVals = NaN(K_current, K_current);
            for i = 1:K_current
                for j = i+1:K_current
                    inter = sum(Gbin(:,i) & Gbin(:,j));
                    union = sum(Gbin(:,i) | Gbin(:,j));
                    if union > 0
                        simVals(i,j) = inter / union;
                    end
                end
            end
            meanJaccardSimilarity(t, rep) = nanmean(simVals(~isnan(simVals)));
        end
    end

    structResults.meanUsageOverlap = meanUsageOverlap;
    structResults.meanJaccardSimilarity = meanJaccardSimilarity;

    saveName = fullfile(saveFolder, ['structureMetrics_' matFiles(k).name]);
    save(saveName, 'structResults', '-v7.3');
    fprintf('Saved structural metrics: %s\n', saveName);
end

%% Plot Figures

clearvars
clc

dataFolder = 'structuralMetrics_runNo_1_geometric_burnin0';
figFolder = 'fig_multicell_structuralmetrics';
if ~exist(figFolder, 'dir')
    mkdir(figFolder);
end

matFiles = dir(fullfile(dataFolder, 'structureMetrics_*.mat'));

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

    fig = figure('Position', [100, 100, 1400, 350 * numDeltas]);

    for row = 1:numDeltas
        load(fullfile(selectedFiles(row).folder, selectedFiles(row).name), 'structResults');

        x = (1:size(structResults.meanUsageOverlap, 1))';
        numReps = size(structResults.meanUsageOverlap, 2);
        rowBase = (row - 1) * 2;

        % --- Column 1: Mean Usage Overlap ---
        subplot(numDeltas, 2, rowBase + 1); hold on;
        for j = 1:numReps
            plot(x, structResults.meanUsageOverlap(:, j), 'Color', [0 0 1 0.2]);
        end
        meanVal = nanmean(structResults.meanUsageOverlap, 2);
        plot(x, meanVal, ':', 'Color', [0 0 1], 'LineWidth', 2);
        ylabel(['\DeltaE = ', deltaForFiles{row}]);
        if row == numDeltas, xlabel('Step'); end
        title('Mean Usage Overlap'); ylim([0, 1]);

        % --- Column 2: Mean Jaccard Similarity ---
        subplot(numDeltas, 2, rowBase + 2); hold on;
        for j = 1:numReps
            plot(x, structResults.meanJaccardSimilarity(:, j), 'Color', [0 0 0 0.2]);
        end
        meanVal = nanmean(structResults.meanJaccardSimilarity, 2);
        plot(x, meanVal, ':k', 'LineWidth', 2);
        if row == numDeltas, xlabel('Step'); end
        title('Mean Jaccard Similarity'); ylim([0, 1]);
    end

    saveas(fig, fullfile(figFolder, sprintf('Multicell_structuralMetrics_%s.png', selLabel)));
    close(fig);
end

disp('Structural metric plots complete.');
