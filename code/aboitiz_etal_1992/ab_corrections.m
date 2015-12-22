% Aboitiz et al (1992) claim 20% of fibers were not imaged in their light
% microscope readings.  That kind of gross correction is OK for the whole
% callosum, but not helpful for individual regions, where thin fibers may
% be greater or lesser.
%
% Aboitiz et al. did provide histograms for five callosal regions, however.
%  Using those regions, and their approximate weightings based on figure
%  1a, we can try to reproportion things appropriately.

project_path = fullfile(fileparts(which(mfilename)), '..', '..');
addpath(genpath(fullfile(project_path, 'code')));
vars = ab_data(false, false);
close all;

% Save the variables by decomposing the struct,
%   assigning the vars locally, saving,
%   then cleaning up.
varnames = fieldnames(vars);
varvals = struct2cell(vars);
for vi=1:length(varnames)
    eval(sprintf('%s = varvals{%d};', varnames{vi}, vi));
end;

% Step 0: normalize histogram data
total_data = repmat(sum(ab_fig4_data, 2), [1 size(ab_fig4_data, 2)]);
ab_fig4_data_normd = ab_fig4_data ./ total_data;

% Choose histogram bars <= 0.4um (use something larger, for round-off
% error)
pct_unimaged = sum(ab_fig4_data_normd(:, ab_fig4_xbin_vals <= 0.401), 2);

% Make sure this corresponds to the reported 20% by weighting these values
% by the estimated area proportions, then making sure they add up to 20% of
% the total.

%Normalize by total % of fibers (across all samples)
%data.ab_factors = ab_factors / sum(ab_factors(:));
computed_pct = ab_fig4_cc_rel_areas * pct_unimaged;
ab_factors = pct_unimaged * (0.20 / computed_pct);

% Redistribute.  According to fig4, ditribute as:
% genu: G1,G2
% ant body: G3,B1
% mid body: B2
% post body: B3,I
% splenium: S1-S3
ab_fig1_all_dens_vals = ab_fig1_g04_dens_vals;
ab_fig1_all_dens_vals(1:2)  = ab_fig1_all_dens_vals(1:2)  * (1 + ab_factors(1));
ab_fig1_all_dens_vals(3:4)  = ab_fig1_all_dens_vals(3:4)  * (1 + ab_factors(2));
ab_fig1_all_dens_vals(5)    = ab_fig1_all_dens_vals(5)    * (1 + ab_factors(3));
ab_fig1_all_dens_vals(6:7)  = ab_fig1_all_dens_vals(6:7)  * (1 + ab_factors(4));
ab_fig1_all_dens_vals(8:10) = ab_fig1_all_dens_vals(8:10) * (1 + ab_factors(5));

save(fullfile(project_path, 'matfiles', 'aboitiz_etal_1992_with_imaging_correction.mat'));
