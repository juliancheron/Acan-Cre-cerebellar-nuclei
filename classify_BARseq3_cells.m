%% BARseq3 cell classification pipeline
% Acan-P2A-Cre cerebellar nuclei dataset
%
% This script classifies cells into five categories based on
% marker gene expression thresholds.
%
% Categories:
% 1 = A cells
% 2 = B cells
% 3 = Inhibitory neurons
% 4 = Non-neuronal cells
% 5 = Unassigned
%
% Input:
%   BARseq3_expression_matrix.csv
%
% Output:
%   BARseq3_cells_with_categories.csv

clear;
clc;

%% Load dataset
data = readtable('BARseq3_expression_matrix.csv');

%% Classification thresholds

thresholds = struct();

% A cells
thresholds.A_Fnbp1l = 1;
thresholds.A_Slc6a1 = 1;
thresholds.A_Sv2c_max = 2;

% B cells
thresholds.B_Acan_min = 1;
thresholds.B_Sv2c_low = 3;
thresholds.B_Sv2c_high = 4;

% Common excitatory markers
thresholds.Common_Slc17a6_A = 3;
thresholds.Common_Slc17a6_B = 4;

% Inhibitory neurons
thresholds.GAD1 = 4;
thresholds.Inhibitory_Slc6a1 = 3;

% Non-neuronal markers
thresholds.Astrocyte_Gfap = 21;
thresholds.Microglia_Tmem119 = 8;

% OPC markers
thresholds.OPC_Pdgfra = 1;
thresholds.OPC_Slc17a6_max = 1;
thresholds.OPC_Sv2c_max = 1;

%% Initialize output vector

cell_categories = zeros(size(data, 1), 1);

%% Cell classification

for i = 1:size(data, 1)

    gene_values = table2array(data(i, ...
        {'Fnbp1l','Acan','Slc6a1','Sv2c','Slc17a6','Gad1','Gfap','Tmem119'}));

    pdgfra_value = data.Pdgfra(i);

    %% 1. NON-NEURONAL CELLS

    is_non_neuronal = false;

    % OPC-like cells
    if pdgfra_value >= thresholds.OPC_Pdgfra && ...
       gene_values(5) <= thresholds.OPC_Slc17a6_max && ...
       gene_values(4) <= thresholds.OPC_Sv2c_max && ...
       gene_values(6) < thresholds.GAD1 && ...
       gene_values(8) < thresholds.Microglia_Tmem119

        is_non_neuronal = true;
    end

    % Astrocyte-like cells
    if gene_values(7) >= thresholds.Astrocyte_Gfap && ...
       gene_values(8) < thresholds.Microglia_Tmem119

        is_non_neuronal = true;
    end

    % Microglia-like cells
    if gene_values(8) >= thresholds.Microglia_Tmem119 && ...
       gene_values(5) < 1 && ...
       gene_values(4) < thresholds.B_Sv2c_high

        is_non_neuronal = true;
    end

    if is_non_neuronal
        cell_categories(i) = 4;
        continue;
    end

    %% 2. A CELLS

    if gene_values(1) >= thresholds.A_Fnbp1l && ...
       gene_values(5) >= thresholds.Common_Slc17a6_A && ...
       gene_values(4) <= thresholds.A_Sv2c_max && ...
       gene_values(8) < thresholds.Microglia_Tmem119

        cell_categories(i) = 1;
        continue;
    end

    %% 3. B CELLS

    if gene_values(5) >= thresholds.Common_Slc17a6_B

        if (gene_values(4) >= thresholds.B_Sv2c_low && ...
            gene_values(2) >= thresholds.B_Acan_min) || ...
           (gene_values(4) >= thresholds.B_Sv2c_high) || ...
           (gene_values(2) >= 2 && gene_values(4) >= 1)

            cell_categories(i) = 2;
            continue;
        end
    end

    %% 4. INHIBITORY NEURONS

    if gene_values(6) >= thresholds.GAD1 && ...
       gene_values(3) >= thresholds.Inhibitory_Slc6a1 && ...
       gene_values(5) <= thresholds.Common_Slc17a6_A && ...
       gene_values(8) < thresholds.Microglia_Tmem119

        cell_categories(i) = 3;
        continue;
    end

    %% 5. UNASSIGNED

    cell_categories(i) = 5;
end

%% GFP analysis

gfp_values = table2array(data(:, 'GFP'));

fprintf('\n--- GFP proportions (threshold = 15) ---\n');

for category = 1:5

    mask = (cell_categories == category);
    n = sum(mask);

    if n > 0

        gfp_positive = sum(gfp_values(mask) >= 15);

        fprintf('Category %d: %.2f%% (%d/%d)\n', ...
            category, ...
            100 * gfp_positive / n, ...
            gfp_positive, ...
            n);
    end
end

%% GFP threshold curves

gfp_range = 0:1:100;

figure;
hold on;

colors = lines(5);

labels = { ...
    'A cells', ...
    'B cells', ...
    'Inhibitory neurons', ...
    'Non-neuronal cells', ...
    'Unassigned'};

for category = 1:5

    mask = (cell_categories == category);
    gfp_subset = gfp_values(mask);

    if isempty(gfp_subset)
        continue;
    end

    counts = zeros(size(gfp_range));

    for t = 1:length(gfp_range)
        counts(t) = sum(gfp_subset >= gfp_range(t));
    end

    % Normalize to proportions
    counts = counts / length(gfp_subset);

    plot(gfp_range, counts, ...
        'LineWidth', 2, ...
        'Color', colors(category,:));
end

xlabel('GFP threshold');
ylabel('Proportion of GFP-positive cells');

title('GFP threshold curves by cell category');

legend(labels, 'Location', 'northeast');

grid on;

% Selected GFP threshold
xline(15, '--k', 'Threshold = 15');

hold off;

%% Save classified dataset

data.CellCategory = cell_categories;

writetable(data, 'BARseq3_cells_with_categories.csv');

fprintf('\nClassified dataset saved successfully.\n');