function vars = ab_data(validate_data)
%
%
% Figures:
%   Figure 1

    if ~exist('validate_data', 'var'), validate_data = true; end;
    if ~exist('visualize_data', 'var'), visualize_data = false; end;

    AB_dirpath = fileparts(which(mfilename));
    AB_dirname = guru_fileparts(AB_dirpath, 'name');
    AB_img_dirpath = fullfile(AB_dirpath, '..', '..', 'img', AB_dirname);


    %% Create diameter => weight mapping function
    [img,orig_img] = scrub_image(fullfile(AB_img_dirpath, 'Fig1b-1.png'), -0.1);

    % Find x axes
    [gx,gxs] = get_groups(sum(img,2) > size(img,2)*0.75, 'right'); xaxis_last = gx(2); xaxis_first = gx(2)-gxs(2);
    [gy,gys] = get_groups(sum(img,1) > size(img,1)*0.75, 'left');  yaxis_first = gy(2); yaxis_last = gy(2)+gys(2);

    data_pix = img(gx(1)+2:xaxis_first-2, gy(1)+gys(1)+1:yaxis_first-1);

    % Get x ticks after that
    % find x and y ticks
    yticks = get_groups(data_pix(:,end-2));
    xticks = get_groups(data_pix(end, :))-1;

    % Pre-plot image
    if visualize_data
        figure; imshow(data_pix); hold on;
        plot(xticks,(size(data_pix,1)-1)*ones(size(xticks)),'*g');
        plot(size(data_pix,2)-1, yticks, '*g');
    end;
    
    dens_pix = zeros(size(xticks));
    for xi = 1:length(xticks)
        xsliver = data_pix(:,max(1,floor(xticks(xi)-3)):min(size(data_pix,2),floor(xticks(xi)+3)));
        xsum    = sum(xsliver,2);
        gg = get_groups(xsum==max(xsum));
        if visualize_data
            plot(xticks(xi), gg, 'y*');
        end;
        if length(gg)~=3, error('?');end;
        dens_pix(xi) = gg(2); % we subtracted off 4 pixels above
    end;

    % Post-plot results
    if visualize_data
        plot(xticks, dens_pix, 'r*')
    end;

    % Get final values
    ab_fig1_g04_dens_vals = (3 + 0.5*(yticks(end)-dens_pix)/mean(diff(yticks)))*1E5;
    ab_fig1_cc_regions = {'G1' 'G2' 'G3' 'B1' 'B2' 'B3' 'I' 'S1' 'S2' 'S3'};

    if validate_data
        % Plot old vs. new values
        figure; set(gcf,'Position',[104         516        1093         158]);
        subplot(1,2,1);
        imshow(orig_img);
        subplot(1,2,2);
        plot(1:length(ab_fig1_cc_regions), ab_fig1_g04_dens_vals, '-o');
        set(gca, 'xtick',1:length(ab_fig1_cc_regions), 'xticklabel', ab_fig1_cc_regions, 'xlim', [0.5 length(ab_fig1_cc_regions)+0.5], 'ylim', 1E5*[2.75 4.45]);
    end;


    %% Fig 1b #2
    [img,orig_img] = scrub_image(fullfile(AB_img_dirpath, 'Fig1b-2.png'), -0);

    % Find x axes
    [gx,gxs] = get_groups(sum(img,2) > size(img,2)*0.75, 'right'); xaxis_last = gx(2); xaxis_first = gx(2)-gxs(2);
    [gy,gys] = get_groups(sum(img,1) > size(img,1)*0.75, 'left');  yaxis_first = gy(2); yaxis_last = gy(2)+gys(2);

    data_pix = img(gx(1)+2:xaxis_first-2, gy(1)+gys(1)+1:yaxis_first-2);
    j = get_groups(sum(data_pix,2),'right');
    data_pix = data_pix(j(1)+1:end,:);

    % Get x ticks after that
    % find x and y ticks
    yticks = get_groups(data_pix(:,end-2));
    xticks = get_groups(data_pix(end, :));

    if visualize_data
        figure; imshow(data_pix); hold on;
        plot(xticks,(size(data_pix,1)-1)*ones(size(xticks)),'*g');
        plot(size(data_pix,2)-1, yticks, '*g');
    end;
    
    dens_pix = zeros(size(xticks));
    for xi = 1:length(xticks)
        xsliver = data_pix(:, max(1,floor(xticks(xi)-3)):min(size(data_pix,2),floor(xticks(xi)+3)));
        xsum    = sum(xsliver,2);
        [gg,ggs] = get_groups(xsum==max(xsum));
        if visualize_data
            plot(xticks(xi), gg, 'y*');
        end;
        if length(gg)~=3, xi, gg, error('?'); end;
        dens_pix(xi) = gg(2); % we subtracted off 4 pixels above
    end;
    if visualize_data
        plot(xticks, dens_pix, 'r*')
    end;
    
    % Get final values
    ab_fig1_g1_dens_vals = (8 + 1*(yticks(end)-dens_pix)/mean(diff(yticks)))*1E4;

    % Plot old vs. new values
    if validate_data
        figure; set(gcf,'Position',[104         516        1093         158]);
        subplot(1,2,1);
        imshow(orig_img);
        subplot(1,2,2);
        plot(1:length(ab_fig1_cc_regions), ab_fig1_g1_dens_vals, '-o');
        set(gca, 'xtick',1:length(ab_fig1_cc_regions), 'xticklabel', ab_fig1_cc_regions, 'xlim', [0.5
        length(ab_fig1_cc_regions)+0.5], 'ylim', 1E4*[7.5 9.5]);
    end;



    %% Fig 1b #3
    [img,orig_img] = scrub_image(fullfile(AB_img_dirpath, 'Fig1b-3.png'), -0);

    % Find x axes
    [gx,gxs] = get_groups(sum(img,2) > size(img,2)*0.75, 'right'); xaxis_last = gx(2); xaxis_first = gx(2)-gxs(2);
    [gy,gys] = get_groups(sum(img,1) > size(img,1)*0.75, 'left');  yaxis_first = gy(2); yaxis_last = gy(2)+gys(2);

    data_pix = img(gx(1)+2:xaxis_first-2, gy(1)+gys(1)+1:yaxis_first-2);
    j = get_groups(sum(data_pix,2),'right');
    data_pix = data_pix(j(1)+1:end,:);

    % Get x ticks after that
    % find x and y ticks
    yticks = get_groups(data_pix(:,end-2));
    xticks = get_groups(data_pix(end, :));

    if visualize_data
        figure; imshow(data_pix); hold on;
        plot(xticks,(size(data_pix,1)-1)*ones(size(xticks)),'*g');
        plot(size(data_pix,2)-1, yticks, '*g');
    end;
 
    dens_pix = zeros(size(xticks));
    for xi = 1:length(xticks)
        xsliver = data_pix(:, max(1,floor(xticks(xi)-3)):min(size(data_pix,2),floor(xticks(xi)+3)));
        xsum    = sum(xsliver,2);
        [gg,ggs] = get_groups(xsum==max(xsum));
        if visualize_data
            plot(xticks(xi), gg, 'y*');
        end;
        if length(gg) > 3, gg = gg(end-2:end); end; % picked up the text
        if length(gg)~=3, xi, gg, error('?'); end;
        dens_pix(xi) = gg(2); % we subtracted off 4 pixels above
    end;
    if visualize_data
        plot(xticks, dens_pix, 'r*')
    end;
    
    % Get final values
    ab_fig1_g3_dens_vals = (5 + 5*(yticks(end)-dens_pix)/mean(diff(yticks)))*1E2;

    % Plot old vs. new values
    if validate_data
        figure; set(gcf,'Position',[104         516        1093         158]);
        subplot(1,2,1);
        imshow(orig_img);
        subplot(1,2,2);
        plot(1:length(ab_fig1_cc_regions), ab_fig1_g3_dens_vals, '-o');
        set(gca, 'xtick',1:length(ab_fig1_cc_regions), 'xticklabel', ab_fig1_cc_regions, 'xlim', [0.5
        length(ab_fig1_cc_regions)+0.5], 'ylim', 1E2*[0 10.]);
    end;



    %% Fig 1b #4
    [img,orig_img] = scrub_image(fullfile(AB_img_dirpath, 'Fig1b-4.png'), -0);

    % Find x axes
    [gx,gxs] = get_groups(sum(img,2) > size(img,2)*0.75, 'right'); xaxis_last = gx(2); xaxis_first = gx(2)-gxs(2);
    [gy,gys] = get_groups(sum(img,1) > size(img,1)*0.75, 'left');  yaxis_first = gy(2); yaxis_last = gy(2)+gys(2);

    data_pix = img(gx(1)+2:xaxis_first-2, gy(1)+gys(1)+1:yaxis_first-2);
    j = get_groups(sum(data_pix,2),'right');
    data_pix = data_pix(j(1)+1:end,:);

    % Get x ticks after that
    % find x and y ticks
    yticks = get_groups(data_pix(:,end-2));
    xticks = get_groups(data_pix(end, :));

    if visualize_data
        figure; imshow(data_pix); hold on;
        plot(xticks,(size(data_pix,1)-1)*ones(size(xticks)),'*g');
        plot(size(data_pix,2)-1, yticks, '*g');
    end;

    dens_pix = zeros(size(xticks));
    for xi = 1:length(xticks)
        xsliver = data_pix(:, max(1,floor(xticks(xi)-3)):min(size(data_pix,2),floor(xticks(xi)+3)));
        xsum    = sum(xsliver,2);
        [gg,ggs] = get_groups(xsum==max(xsum));
        if visualize_data
            plot(xticks(xi), gg, 'y*');
        end;
        if length(gg) > 3, gg = gg(end-2:end); % picked up the text
        elseif length(gg) == 1, gg=[0 gg]; end;
        if visualize_data
            ens_pix(xi) = gg(2); % we subtracted off 4 pixels above
        end;
        dens_pix(xi) = gg(2); % we subtracted off 4 pixels above
    end;
    if visualize_data
        plot(xticks, dens_pix, 'r*')
    end;
    
    % Get final values
    ab_fig1_g5_dens_vals = (1 + 1*(yticks(end)-dens_pix)/mean(diff(yticks)))*1E2;

    % Plot old vs. new values
    if validate_data
        figure; set(gcf,'Position',[104         516        1093         158]);
        subplot(1,2,1);
        imshow(orig_img);
        subplot(1,2,2);
        plot(1:length(ab_fig1_cc_regions), ab_fig1_g5_dens_vals, '-o');
        set(gca, 'xtick',1:length(ab_fig1_cc_regions), 'xticklabel', ab_fig1_cc_regions, 'xlim', [0.5
        length(ab_fig1_cc_regions)+0.5], 'ylim', 1E2*[0 3.5]);
    end;


    % Determined this through manual image editing interventions
    ab_fig1_cc_areas = [27554 12007 12106 7952 7806 9128 8969 10858 15946 13889];
    ab_fig1_cc_rel_areas = ab_fig1_cc_areas./sum(ab_fig1_cc_areas);
    
    
    %% Construct outputs
    varnames = who('ab_*');
    varvals = cellfun(@eval, varnames, 'UniformOutput', false);
    vars = cell2struct(varvals, varnames);
