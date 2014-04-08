function vars = lut_data(validate_data, visualize_data)
%
% Luttenberg (1965)
% Human developmental data on corpus callosum

    if ~exist('validate_data', 'var'), validate_data = true; end;
    if ~exist('visualize_data', 'var'), visualize_data = false; end;


    %% Collect data
    lut_table1_age = [[10 11 13 13.5 14 15 17 19 21 26 30 36 40]*7 (9*30+[3 5]*30)];
    lut_table1_crlen = [53 62 89 102 110 124 149 175 200 250 290 330 365 nan nan];
    lut_table1_density = [24 32 37 44 49 52 62 65 65 65 62 62 61 59 59
                24 31 33 42 49 53 62 65 65 65 62 61 61 59 59
                24 33 36 42 50 52 62 65 65 65 62 61 60 59 59];
    %lut_table2_area = [nan nan nan 0.5 1.275 5.108 4.97 6.161 11.25 20.093 24.462 48.608 57.87 53.37 66.78
    %          nan nan nan 0.66 3.150 3.784 6.334 15.634 13.265 26.575 34.434 44.527 60.086 70.684 70.394
    %          nan nan nan .425 2.5515 0.945 4.16 5.338 9.158 14.669 26.664 24.304 49.66 44.44 57.78
    %          .05408 .097344 1.563 2.5 10.204 13.23 19.39 35.552 40.816 75.548 106.656 145.824 206.6 217.756 244.42];
    %
    %correct by /0.7?
    %lut_table2_naxons = [nan nan nan 0.220 0.629 2.640 3.063 4.005 7.307 13.011 15.228 29.991 35.110 31.424 39.373
    %                   nan nan nan 0.660 3.150 3.784 6.334 15.634 13.265 26.575 34.434 44.527 60.086 70.684 70.394
    %                   nan nan nan 0.179 1.271 0.489 2.560 3.466 5.958 9.516 16.505 14.825 30.005 26.357 34.159
    %                   0.013 0.031 0.554 1.059 5.05 6.913 11.957 23.105 26.530 49.102 66.167 89.343 125.200 128.465 143.926];

    lut_table2_areas = [nan nan nan 0.05408;
                    nan nan nan 0.097344;
                    nan nan nan 1.563;
                    0.5 1.575 0.425 2.5;
                    1.275 6.3775 2.5515 10.204;
                    5.108 7.177 0.945 13.230;
                    4.97 10.26 4.16 19.39;
                    6.161 24.053 5.338 35.552;
                    11.25 20.408 9.158 40.816;
                    20.093 40.84 14.669 75.548
                    24.462 55.53 26.664 106.656;
                    48.608 72.912 24.304 145.824;
                    57.87 99.07 49.66 206.6;
                    53.37 119.946 44.44 217.756;
                    66.78 119.86 57.78 244.42];
    lut_table2_nfibers = 1E6*[nan nan nan 0.013;
                          nan nan nan 0.031;
                          nan nan nan 0.554;
                          .220 .660 .1798 1.059;
                          .629 3.150 1.271 5.050;
                          2.640 3.784 0.489 6.913;
                          3.063 6.334 2.560 11.957;
                          4.005 15.634 3.466 23.105;
                          7.307 13.265 5.958 26.530;
                          13.011 26.575 9.516 49.102;
                          15.228 34.434 16.505 66.167;
                          29.991 44.527 14.825 89.343;
                          35.110 60.086 30.005 125.200;
                          31.424 70.684 26.357 128.465;
                          39.373 70.394 34.159 143.926];

    lut_table3_diameters = [0.5 1.0 1.5 2.0 2.5 3.0];
    lut_table3_genu = [66.5 30.5 3 0 0 0;
                   62.5 23.6 11.5 2.4 0 0;
                   77.5 16.1 5 1.4 0 0;
                   52.2 31.7 13.2 2.9 0 0;
                   47.9 41.2 8.1 2.8 0 0;
                   52.7 39.6 6.5 1.2 0 0 ;
                   45 49.9 4.6 0.5 0 0;
                   46.7 50.8 2.2 0.3 0 0;
                   38.1 49.5 11.9 0.5 0 0;
                   38.5 41.5 16.1 2.9 1 0;
                   37.7 39.1 18.4 4 0.8 0;
                   42.4 46.8 9.2 1.2 0.4 0;
                   28.4 33.2 27.1 8.2 2.8 0.3;
                   35.4 36.3 15.2 8.1 3.5 1.5;
                   37.7 40.4 15.1 5 1.3 0.5];
    lut_table3_truncus = [82 16 2 0 0 0;
                      74.3 20.2 5.5 0 0 0;
                      76.4 18.1 4.7 0.8 0 0;
                      55.6 38.4 4.5 1.5 0 0;
                      42 45.4 11.9 0.7 0 0;
                      62.6 35.4 1.6 0.4 0 0;
                      49.1 48 2.7 0.2 0 0;
                      26.2 64 8.9 0.9 0 0;
                      38.6 49.7 10.7 1 0 0;
                      38 40.9 15.2 4.9 1 0 ;
                      39.7 43.3 13.6 2.7 0.7 0;
                      37.7 43.3 12 5.1 1.9 0;
                      53 38.2 7.5 1.1 0.2 0;
                      26.7 51.3 16.5 4.9 0.6 0;
                      44.5 37.9 10.8 5.3 1 0.5];
    lut_table3_splenium = [78.4 15.8 5.8 0 0 0;
                       75.5 21.7 2.8 0 0 0;
                       72.7 17.4 8.1 1.8 0 0;
                       65 30.8 4 0.2 0 0;
                       56.5 39.6 3.2 0.7 0 0;
                       nan nan nan nan 0 0;
                       50.1 45.2 3.9 0.8 0 0;
                       34.7 57.4 7.8 0.1 0 0;
                       33.8 44 19.8 2.4 0 0;
                       37.3 39.6 18.3 4 0.8 0;
                       40.3 42.9 12.5 3.1 1.2 0;
                       39.9 45.1 12 2 1 0;
                       29.6 33.4 19.2 10.8 4.2 2.8;
                       nan nan nan nan nan nan;
                       37.8 33.4 23.2 4.2 1.1 0.3];


    %% validate data
    if validate_data
        fprintf('All data taken from tables; no data to validate!\n');
    end;

    %% Visualize data
    if visualize_data
        % Figure 2
        lut_subj_idx = [2 9 13];
        lut_area_proportions = lut_table2_areas(:,1:3)./repmat(lut_table2_areas(:,4),[1 3]);
        lut_area_proportions(isnan(lut_area_proportions)) = 1/3;
        lut_distn = lut_table3_genu.*repmat(lut_area_proportions(:,1),[1 size(lut_table3_genu,2)]) ...
                + lut_table3_truncus.*repmat(lut_area_proportions(:,2),[1 size(lut_table3_truncus,2)]) ...
                + lut_table3_splenium.*repmat(lut_area_proportions(:,3),[1 size(lut_table3_splenium,2)]);

        lut_density = sum(lut_area_proportions.*lut_table1_density', 2);

        % Reproduce figure 2
        figure; set(gca, 'FontSize', 14);
        title('Luttenberg (1965), Figure 2');
        hold on;
        plot([0 lut_table3_diameters], [zeros(length(lut_subj_idx),1) lut_distn(lut_subj_idx,:)]', 'LineWidth', 2);
        legend({'11 weeks','21 weeks','40 weeks (newborn)'})
        xlabel('axon diameter ({\mu}m)'); ylabel('percent');

        % Show all, as a waterfall plot
        [Y,X] = meshgrid([0 lut_table3_diameters], lut_table1_age);
        figure;
        waterfall(X,Y,[zeros(size(lut_distn,1),1) lut_distn]);
    end;


    %% Construct outputs
    varnames = who('lut_*');
    varvals = cellfun(@eval, varnames, 'UniformOutput', false);
    vars = cell2struct(varvals, varnames);
