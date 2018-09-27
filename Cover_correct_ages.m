%% CALCULATE SURFACE COVER-CORRECTED EXPOSURE AGES 
%
% Calculates surface cover shielding factors and resulting exposure ages 
% from sample data in a spreadsheet (.xlsx, .csv) or text file (.txt) using
% selected surface cover type or manual cover density, and a specified 
% cover depth and scaling model.
%
% Plots the corrected exposure ages as a kernel density estimate.
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
% 10. Shielding factor for terrain (unitless)
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

data_name = 'example_moraine_input.xlsx';  % File name used for sample data
ages_name = 'example_cover';  % Name of dataset to be used to save corrected data (e.g. _ages.mat)

% SET scaling model
scaling_model = 'LSD';  % 'DE','DU','LI','ST','LM','LSD','LSDn'

% SET surface cover correction
cover_type = 'snow';  % Select 'snow', 'freshwater', 'seawater', 'loess', 'till', 'soil', 'ash', or 'manual'
cover_density = 1.6;  % If 'manual', set density of surface cover (g/cm^3)

cover_depth = 20;     % Select depth of surface cover (cm)


% Load sample data
sample_data = get_data(data_name);


%% Calculate Exposure Ages Corrected for Surface Cover

% Perform correction
corr_ages_ka = cov_correct(sample_data,scaling_model,cover_type,cover_depth,cover_density);

% Save corrected data to file
save_name = strcat(ages_name,'.mat');
save(save_name,'corr_ages_ka');


%% Plot Corrected Ages as Kernel Density Estimates

% Plot settings
feature = 1;    % Data from single feature?  yes = 1, no = 0
save_plot = 0;  % Save plot?  1 to save as .png and .eps, otherwise 0
mask = [];      % Select samples to plot (default is all)
time_lim = [];  % Optionally set x-axis limits (in ka)
weighted = [];  % Optionally select weighted (1) or unweighted (0) mean and standard deviation (default is weighted)

% Plot figure
corr_ages_ka.ages_name = ages_name;
plot_kernel(corr_ages_ka,feature,save_plot,mask,time_lim,weighted);

