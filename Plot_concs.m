%% PLOT ISOTOPE CONCENTRATIONS FROM 10BE AND 26AL NUCLIDE CONCENTRATIONS
%
% Plots nuclide concentrations using sample data in a spreadsheet (.xlsx, 
% .csv) or text file (.txt). Saves figures as PNG and EPS files.
%
% Sample concentrations can be plotted on a two-isotope diagram (currently,
% 26Al/10Be vs. 10Be), or plotted against depth (designed for cores).
%
% Input data must include the following information:
% 1. Sample name
% 2. Latitude (decimal degrees)
% 3. Longitude (decimal degrees)
% 4. Elevation (m asl)
% 5. Pressure (hPa) (zero if not known)
% 6. Relative position (distance from terminus, km; elevation above ice, m)
% 7. Sample thickness (cm)
% 8. Bulk density (g/cm^3)
% 9. Shielding factor for terrain, snow, etc. (unitless)
% 10. Sample 10-Be concentration (atoms of 10-Be/g)
% 11. Sample 10-Be concentration 1 sigma uncertainty (atoms of 10-Be/g)
% 12. Sample 26-Al concentration (atoms of 26-Al/g)
% 13. Sample 26-Al concentration 1 sigma uncertainty (atoms of 26-Al/g)
% 14. Top depth of sample (cm)
% 15. Bottom depth of sample (cm)
% 16. Final mineral weight (g)
% 17. Year the sample was collected (calendar year)
%
% Written by Richard Selwyn Jones, Durham University
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite.
%
%
%%

clear % Start fresh
addpath(genpath(pwd))

% SET inputs
input_name = 'GISP2_input_complex.xlsx'; % File name used for sample data
cover_z = 0;                        % Cover depth above surface (g/cm^2) - default is zero
scaling_model = 'LSD';              % Scaling model - 'DE','DU','LI','ST','LM','LSD','LSDn'
concs_name = 'GISP2_data';          % Name of dataset to be used to save figures


% Load sample data
sample_data = get_data_complex(input_name); 
sample_data.cover.z = cover_z;

% Get parameters for samples
sample_data = get_pars(sample_data,scaling_model);


%% Plot Concentrations on a Two-Isotope Diagram

% Plot Settings
sigma = 2;      % Set sigma for plotted concentrations (1 or 2)
save_plot = 0;  % Save plot?  1 to save as .png and .eps, otherwise 0
expo_int = [];  % Optionally set exposure intervals to be shown (ka)
bur_int = [];   % Optionally set burial intervals to be shown (ka)
x_lim = [];     % Optionally set x-axis limits (norm. 10Be conc, atoms g)
y_lim = [];     % Optionally set y-axis limits (norm. 26Al/10Be ratio)

% Plot figure
plot_2iso_concs(sample_data,sigma,save_plot,concs_name,expo_int,bur_int,x_lim,y_lim);


%% Plot Concentrations with Depth

% Plot Settings
save_plot = 0;  % Save plot?  1 to save as .png and .eps, otherwise 0
x_lim = [];     % Optionally set x-axis limits (nuclide concentration, atoms g)
y_lim = [];     % Optionally set y-axis limits (depth, m)

% Plot figure
plot_concs_depth(sample_data,save_plot,concs_name,x_lim,y_lim);

