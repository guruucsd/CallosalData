function vars = tr_data(validate_data, visualize_data)
%
% Dataset:
%   Tarpley & Ridgway (1994)
%
% Data:
%
% Figures:
%    Figure 2: brain mass (g) vs. CCA (mm^2)
%
% Notes:
%   Generally processed fresh, "some brains" list 1.5-7.8% mass over time.

    if ~exist('validate_data', 'var'), validate_data = true; end;
    if ~exist('visualize_data', 'var'), visualize_data = false; end;

    RID_dirpath = fileparts(which(mfilename));
    RID_dirname = guru_fileparts(RID_dirpath, 'name');
    RID_img_dirpath = fullfile(RID_dirpath, '..', '..', 'img', RID_dirname);


    %% Fig 2
    [data_ypix,data_xpix] = get_pixels_by_color(fullfile(RID_img_dirpath, 'Fig2_reddot.png'), 'r', -0.05);
    [img] = scrub_image(fullfile(RID_img_dirpath, 'Fig2.png'), -0.05);

    % reorder in xy
    [data_xpix,idx1] = sort(data_xpix); data_ypix = data_ypix(idx1);

    % Get the (b&w) axes
    [yaxis_idx,yaxis_width]=get_groups(sum(img,1)>0.5*size(img,1));
    [xaxis_idx,xaxis_width]=get_groups(sum(img,2)>0.5*size(img,2));

    if visualize_data
        figure;
        imshow(img); hold on;
    end;

    % Grab the ticks
    xticks     = get_groups(img(floor(xaxis_idx-xaxis_width/2)-2,:)); xticks = xticks(3:end);
    yticks = get_groups(img(:,floor(yaxis_idx+yaxis_width/2)+2,:)); %yticks = yticks([1:3 5])

    if visualize_data
        figure;
        imshow(img); hold on;
        plot(yaxis_idx(1), yticks, 'r*');
        plot(xticks, xaxis_idx, 'b*');
    end;

    % Convert the datapoints from pixels to data(x,y) values
    tr_fig2_brain_weights = 10.^(2.6+(data_xpix - xticks(1))/mean(diff(xticks))*0.2)'; % average over cca and wmv; should have same xval!
    tr_fig2_ccas = 10.^(2+(yticks(end)-data_ypix)/mean(diff(yticks))*0.1)';


    %% Validate data
    if validate_data
        fprintf('Data validation NYI\n'); %keyboard
    end;

    if visualize_data
        % Re-create the plot!
        p_cca = polyfit(log10(tr_fig2_brain_weights), log10(tr_fig2_ccas), 1);
        figure; set(gcf, 'Position', [49         290        1194         394]);
        subplot(1,2,1);
        imshow(~img); hold on;
        plot(data_xpix, data_ypix, 'g*');

        subplot(1,2,2);
        plot(log10(tr_fig2_brain_weights), log10(tr_fig2_ccas), 'o');
        hold on; plot(log10([0.1; tr_fig2_brain_weights]), p_cca(1)*log10([0.1; tr_fig2_brain_weights])+ p_cca(2));

        set(gca, 'FontSize', 14, 'xlim', [2.6 4], 'ylim', [2 3]);
        xlabel('log(brain weight) (g)'); ylabel('log(cca) (mm^2)');
        title(sprintf('Regressions: %4.2fX + %4.2f', p_cca));
    end;


    %% Construct outputs
    varnames = who('tr_*');
    varvals = cellfun(@eval, varnames, 'UniformOutput', false);
    vars = cell2struct(varvals, varnames);

