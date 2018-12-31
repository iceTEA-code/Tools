%% CALCULATE AND PLOT EXPOSURE AGES FROM 10BE AND 26AL NUCLIDE CONCENTRATIONS
%
% Calculates exposure ages data from sample data in a spreadsheet (.xlsx, 
% .csv) or text file (.txt) using a specified scaling model. A new file is 
% then created with those ages.
%
% Optionally identifies and removes outliers if ages are from a single
% feature.
%
% Plots those exposure ages as a kernel density estimate, with statistics 
% computed if the ages are for a single feature, or as a (horizontal or 
% vertical) transect. The figures can be automatically saved using the
% given output name (ages_name), nuclide and scaling method.
%
% Input data must include the following information:
% 1. Sample name
% 2. Latitude (decimal degrees)
% 3. Longitude (decimal degrees)
% 4. Elevation (m asl)
% 5. Pressure (hPa) (zero if not known)
% 6. Elevation uncertainty (m) (zero if not known)
% 7. Relative position (distance from terminus, km; elevation above ice, m) (zero or NaN if not relevant or known)
% 8. Sample thickness (cm)
% 9. Bulk density (g/cm^3)
% 10. Shielding factor for terrain, snow, etc. (unitless)
% 11. Sample 10-Be concentration (atoms of 10-Be/g)
% 12. Sample 10-Be concentration 1 sigma uncertainty (atoms of 10-Be/g)
% 13. Sample 26-Al concentration (atoms of 26-Al/g)
% 14. Sample 26-Al concentration 1 sigma uncertainty (atoms of 26-Al/g)
% 15. Year the sample was collected (calendar year)
%
% Optional data should include the following information:
% 16. Sample 10-Be exposure age (mean; years)
% 17. Sample 10-Be exposure 1 sigma uncertainty (internal; years)
% 18. Sample 10-Be exposure 1 sigma uncertainty (external; years)
% 19. Sample 26-Al exposure age (mean; years)
% 20. Sample 26-Al exposure 1 sigma uncertainty (internal; years)
% 21. Sample 26-Al exposure 1 sigma uncertainty (external; years)
% 22. Scaling model used (i.e. 'DE','DU','LI','ST','LM','LSD'/'SF','LSDn'/'SA')
% NOTE: To plot previously calculated ages use Import_Plot_age.m.
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
ages_name = 'example_moraine';  % Name of dataset to be used to save calculated ages (_ages.mat)

% Load sample data
sample_data = get_data(data_name);


%% Calculate and Save Ages
% Run section only if ages need to be calculated

% SET scaling model
scaling_model = 'LSDn'; % 'DE','DU','LI','ST','LM','LSD','LSDn'

% Plot production rate through time?
plot_prod = 0; % yes = 1, no = 0


% Calculate
ages_ka = age_calc(sample_data,scaling_model,plot_prod);
ages_ka.ages_name = ages_name;

% Plot production through time for samples
save_plot = 0;  % Save plot?  1 to save as .png and .eps, otherwise 0
y_lim = [];     % Optionally set y-axis limits (in atoms/g/year)
if plot_prod == 1
    plot_prod_time(ages_ka,save_plot,ages_name,y_lim);
end

% Save ages to file
save_name = strcat(ages_name,'_ages.mat');
save(save_name,'ages_ka');

% Export results table
format = 'xls'; % SET export format - 'xls' or 'txt'
export_calcages(sample_data,ages_ka,ages_name,format);


%% Plot Ages as Kernel Density Estimates

load_name = strcat(ages_name,'_ages.mat');
load(load_name);

% Plot settings
feature = 1;    % Data from single feature?  yes = 1, no = 0
save_plot = 0;  % Save plot?  1 to save as .png and .eps, otherwise 0
mask = [];      % Select samples to plot (default is all)
time_lim = [];  % Optionally set x-axis limits (in ka)
weighted = [];  % Optionally select weighted (1) or unweighted (0) mean and standard deviation (default is weighted)

% Plot figure
plot_kernel(ages_ka,feature,save_plot,mask,time_lim,weighted);


%% Remove Outliers and Re-plot Ages

% Data from single feature?
feature = 1; % yes = 1, no = 0

% Find and remove outliers (generalised ESD test)
sig = []; % Optionally set significance level (default is 0.05)
new_ages_ka = find_outliers(ages_ka,feature,sig);

% Plot figure
plot_outliers = 1; % Also plot outliers?  yes = 1, no = 0
save_plot = 0;  % Save plot?  1 to save as .png and .eps, otherwise 0
mask = [];      % Select samples to plot (default is all)
time_lim = [];  % Optionally set x-axis limits (in ka)
weighted = [];  % Optionally select weighted (1) or unweighted (0) mean and standard deviation (default is weighted)
plot_outlier_kernel(ages_ka,new_ages_ka,plot_outliers,feature,save_plot,mask,time_lim,weighted);


%% Plot Ages as a Transect

load_name = strcat(ages_name,'_ages.mat');
load(load_name);

% Plot settings
transect_type = 'vert'; % SET as 'vert' or 'horiz'
save_plot = 0;  % Save plot?  1 to save as .png and .eps, otherwise 0
mask = [];      % Select samples to plot (default is all)
time_lim = [0 12];  % Optionally set limits of time axis (in ka)
pos_lim = [];   % Optionally set limits of relative position axis (in m or km)

% Plot figure
plot_transect(ages_ka,transect_type,save_plot,mask,time_lim,pos_lim);

