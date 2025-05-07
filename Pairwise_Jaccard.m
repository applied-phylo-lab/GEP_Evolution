clearvars
clc

dataFolder = 'simData_runNo_1_geometric_burnin0';
saveFolder = 'pairwiseJaccardMetrics_runNo_1';
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end

matFiles = dir(fullfile(dataFolder, '*.mat'));

for k = 1:length(matFiles)
    load(fullfile(matFiles(k).folder, matFiles(k).name), 'results', 'params');

    numSteps = params.numSteps;
    numReps = params.numSimulations;

    allPairJaccards = cell(numSteps, numReps);

    for rep = 1:numReps
        genomes = results.allGenomes{rep};

        for t = 1:numSteps
            G = genomes{t};  % [L x K]
            if isempty(G), continue; end
            Gbin = G > 0;
            K_current = size(G, 2);

            Jmat = NaN(K_current, K_current);
            for i = 1:K_current
                for j = i+1:K_current
                    inter = sum(Gbin(:,i) & Gbin(:,j));
                    union = sum(Gbin(:,i) | Gbin(:,j));
                    if union > 0
                        Jmat(i,j) = inter / union;
                        Jmat(j,i) = Jmat(i,j);  % symmetry
                    end
                end
            end
            allPairJaccards{t, rep} = Jmat;
        end
    end

    saveName = fullfile(saveFolder, ['pairwiseJaccards_' matFiles(k).name]);
    save(saveName, 'allPairJaccards', 'params', '-v7.3');
    fprintf('Saved Jaccard trajectories: %s\n', saveName);
end


%%

clearvars
clc

replicateThreshold = 10;

dataFolder = 'pairwiseJaccardMetrics_runNo_1';
figFolder = 'fig_multicell_pairwiseJaccard_meanOnly';
if ~exist(figFolder, 'dir')
    mkdir(figFolder);
end

matFiles = dir(fullfile(dataFolder, 'pairwiseJaccards_*.mat'));

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

    fig = figure('Position', [100, 100, 1600, 350 * numDeltas]);

    for row = 1:numDeltas
        load(fullfile(selectedFiles(row).folder, selectedFiles(row).name), 'allPairJaccards', 'params');

        numSteps = params.numSteps;
        numReps = params.numSimulations;
        K = params.K;
        x = (1:numSteps)';

        pairTrajectories = containers.Map();

        for rep = 1:numReps
            for t = 1:numSteps
                J = allPairJaccards{t, rep};
                if isempty(J), continue; end
                [K_actual, ~] = size(J);
                for i = 1:K_actual
                    for j = i+1:K_actual
                        key = sprintf('%d_%d', i, j);
                        if ~isKey(pairTrajectories, key)
                            pairTrajectories(key) = NaN(numSteps, numReps);
                        end
                        mat = pairTrajectories(key);
                        mat(t, rep) = J(i,j);
                        pairTrajectories(key) = mat;
                    end
                end
            end
        end

        subplot(numDeltas, 1, row); hold on;
        keys = pairTrajectories.keys;
        for kidx = 1:length(keys)
            mat = pairTrajectories(keys{kidx});
            count = sum(~isnan(mat), 2);
            mval = nanmean(mat, 2);
            mval(count < replicateThreshold) = NaN;
            plot(x, mval, '-', 'LineWidth', 1.2);
        end
        ylabel(['\DeltaE = ', deltaForFiles{row}]);
        title('Mean Pairwise GEP Jaccard Overlap (Per GEP Pair)');
        if row == numDeltas, xlabel('Mutation Steps'); end
        ylim([0, 0.2]);
    end

    saveas(fig, fullfile(figFolder, sprintf('Multicell_PairwiseJaccard_meanOnly_%s.png', selLabel)));
    close(fig);
end
