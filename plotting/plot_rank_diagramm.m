% Load processed results
loaded = load('results/processed/fig4/processed_results.mat');
processed = loaded.processed;

% Access BOUND rank frequency
scenario = processed.scenarios.n6;
rank_freq = scenario.methods.BOUND.metrics.rank_frequency;
k_list = scenario.k_list;

% Determine all unique ranks across all K values
all_ranks = [];
for k_idx = 1:length(k_list)
    data = rank_freq{k_idx};
    all_ranks = union(all_ranks, data.ranks);
end
all_ranks = sort(all_ranks);

% Build frequency matrix: rows = ranks, columns = K values
freq_matrix = zeros(length(all_ranks), length(k_list));
for k_idx = 1:length(k_list)
    data = rank_freq{k_idx};
    for r = 1:length(data.ranks)
        rank_idx = find(all_ranks == data.ranks(r));
        freq_matrix(rank_idx, k_idx) = data.frequencies(r);
    end
end

% Create grouped bar chart
figure('Position', [100, 100, 700, 500]);
bar(all_ranks, freq_matrix);
xlabel('Rank of X_{opt}', 'FontSize', 12);
ylabel('Relative frequency', 'FontSize', 12);
grid on;

% Create legend with K values
legend_labels = arrayfun(@(k) sprintf('%d users', k), k_list, 'UniformOutput', false);
legend(legend_labels, 'Location', 'northeast');

% Set x-axis to show integer ranks only
xticks(all_ranks);

% Optional: Save figure
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
output_dir = fullfile(project_root, 'results', 'figures', processed.meta.experiment_name);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
saveas(gcf, fullfile(output_dir, 'rank_distribution.png'));
fprintf('Figure saved to: %s\n', fullfile(output_dir, 'rank_distribution.png'));