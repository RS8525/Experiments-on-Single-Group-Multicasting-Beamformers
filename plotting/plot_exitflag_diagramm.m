% Load processed results
loaded = load('results/processed/fig5/processed_results.mat');
processed = loaded.processed;

% Access SQP exitflag frequency
scenario = processed.scenarios.n4;
exitflag_freq = scenario.methods.sqp_beamformer.metrics.exitflag_frequency;
k_list = scenario.k_list;

% Determine all unique exitflags across all K values
all_flags = [];
for k_idx = 1:length(k_list)
    data = exitflag_freq{k_idx};
    all_flags = union(all_flags, data.flags);
end
all_flags = sort(all_flags);

% Build frequency matrix: rows = exitflags, columns = K values
freq_matrix = zeros(length(all_flags), length(k_list));
for k_idx = 1:length(k_list)
    data = exitflag_freq{k_idx};
    for f = 1:length(data.flags)
        flag_idx = find(all_flags == data.flags(f));
        freq_matrix(flag_idx, k_idx) = data.frequencies(f);
    end
end

% Create grouped bar chart
figure('Position', [100, 100, 700, 500]);
bar(all_flags, freq_matrix);
xlabel('Exitflag', 'FontSize', 12);
ylabel('Relative frequency', 'FontSize', 12);
grid on;

% Create legend with K values
legend_labels = arrayfun(@(k) sprintf('%d users', k), k_list, 'UniformOutput', false);
legend(legend_labels, 'Location', 'northeast');

% Set x-axis to show integer exitflags only
xticks(all_flags);

% Add title with exitflag meanings
title('SQP Beamformer Convergence (1=success, 0=iter limit, <0=failed)', 'FontSize', 10);

% Optional: Save figure
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
output_dir = fullfile(project_root, 'results', 'figures', processed.meta.experiment_name);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
saveas(gcf, fullfile(output_dir, 'exitflag_distribution.png'));
fprintf('Figure saved to: %s\n', fullfile(output_dir, 'exitflag_distribution.png'));
