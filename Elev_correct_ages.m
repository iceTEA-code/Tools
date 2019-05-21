%% CALCULATE ELEVATION-CORRECTED EXPOSURE AGES 
%
% Calculates time-dependent scaling factors and exposure ages from sample 
% data in a spreadsheet (.xlsx, .csv) or text file (.txt) using selected 
% glacial isostatic adjustment (GIA) model or rate of elevation change, and
% a specified scaling model.
%
% Plots the elevation-corrected production through time, and the corrected 
% exposure ages as a kernel density estimate.
%
% Input data must include the following information:
% 1. Sample name
% 2. Latitude (decimal degrees)
% 3. Longitude (decimal degrees)
% 4. Elevation (m asl)
% 5. Pressure (hPa) (zero if not known)
% 6. Elevation uncertainty (m) (zero if not known)
% 7. Relative position (NaN as not relevant here)
% 8. Sample thickness (cm)
% 9. Bulk density (g/cm^3)
% 10. Shielding factor for terrain, snow, etc. (unitless)
% 11. Sample 10-Be concentration (atoms of 10-Be/g)
% 12. Sample 10-Be concentration 1 sigma uncertainty (atoms of 10-Be/g)
% 13. Sample 26-Al concentration (atoms of 26-Al/g)
% 14. Sample 26-Al concentration 1 sigma uncertainty (atoms of 26-Al/g)
% 15. Year the sample was collected (calendar year)
%
% For improved computation time, the age calculation determines age 
% uncertainties based only on the elevation (pressure) and measurement 
% uncertainties, in addition to the inherent production rate uncertainty.
% Note: the calculations take approximately three times longer for 
% time-dependent scaling schemes.
%
% Written by Richard Selwyn Jones, Durham University
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite, which is built on versions of 
% CRONUS-Earth and CRONUScalc code.
%
%
%%

clear % Start fresh
addpath(genpath(pwd))

data_name = 'elev_examples_input.xlsx';  % File name used for sample data
ages_name = 'examples_corr_i6g';  % Name of dataset to be used to save corrected data (e.g. _ages.mat)

% SET scaling model
scaling_model = 'LSD';      % 'DE','DU','LI','ST','LM','LSD','LSDn'

% SET elevation correction
correction_type = 'model';  % Select 'model' or 'rate'
GIA_model = 'W12';         % If 'model', set GIA model to use - 'I5G', 'I6G' or 'W12'(Antarctica only)
elev_rate = 2;             % If 'rate', set rate of elevation change (m/ka) - 
                           % A positive/negative value would correspond to lower/higher elevations in the past (uplift/subsidence)

% Load sample data
sample_data = get_data(data_name);


%% Calculate Elevation-corrected Exposure Ages

if strcmp(correction_type,'model')
    elev_input = GIA_model;
elseif strcmp(correction_type,'rate')
    elev_input = elev_rate;
end

% Perform correction
corrected = elev_correct(sample_data,scaling_model,correction_type,elev_input);

% Plot production through time
save_plot = 0;  % Save plot?  1 to save as .png and .eps, otherwise 0
y_lim = [];     % Optionally set y-axis limits (in atoms/g/year)
plot_prod_time(corrected.plot.uncorr,save_plot,ages_name,y_lim,corrected.plot.corr);

% Save corrected data to file
save_name = strcat(ages_name,'.mat');
save(save_name,'corrected');

% Export results table
format = 'xls'; % SET export format - 'xls' or 'txt'
export_elevcorr_results(sample_data,corrected,ages_name,format);


%% Plot Corrected Ages and Uncorrected Ages as Kernel Density Estimates

% Load uncorrected ages
% (Need to have previously used Calc_Plot_age.m or Import_Plot_age.m to obtain and save age data)
uncorrected_ages_name = 'elev_examples';    % Name used for saved dataset (_ages.mat)
load(strcat(uncorrected_ages_name,'_ages.mat'));

% Plot settings
mask = [];      % Select samples to plot (default is all)
x_lim = [];     % Optionally set x-axis limits (in ka)
save_plot = 0;  % Save plot?  1 to save as .png and .eps, otherwise 0

% Plot and Save figure
uncorr_ages_ka = ages_ka;
corr_ages_ka = corrected.plot.corr; corr_ages_ka.ages_name = ages_name;
plot_corr_kernel(uncorr_ages_ka,corr_ages_ka,save_plot,mask,x_lim);

