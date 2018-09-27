%% ANALYSE LINEAR RATES OF RETREAT OR THINNING
%
% Determines linear estimates of retreat/thinning rates for exposure age
% data in a horizontal or vertical transect. Least-squares regression is 
% applied randomly to normally-distributed exposure ages (2 sigma) through 
% a Monte Carlo simulation.
%
% Plots the probability distribution of computed rates, with estimates at
% 68% and 95% confidence bounds, and the models and bounds as a transect,
% with or without corresponding exposure ages.
%
% Exposure ages need to be in a correctly structured Matlab file 
% (_ages.mat). This can be done using Calc_Plot_age.m or Import_Plot_age.m,
% to calculate ages or import existing ages, respectively.
%
% Written by Richard Selwyn Jones, Durham University
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite.
%
%
%%

clear % Start fresh
addpath(genpath(pwd))

ages_name = 'example_retreattransect'; % Name of dataset containting calculated ages (_ages.mat)

load_name = strcat(ages_name,'_ages.mat');
load(load_name);


%% Determine Rates

transect_type = 'horiz'; % SET as 'vert' or 'horiz'
regress_type = 'weighted'; % SET least-square regression as 'unweighted' or 'weighted'
n_iter = 10000; % SET number of Monte Carlo iterations (default is 5000)
mask = []; % Select samples to analyse (leave empty for all samples; [])
save_plot = 0; % Save plot?  1 to save as .png and .eps, otherwise 0

regressed_rates = transect_regress_linear(ages_ka,transect_type,regress_type,save_plot,n_iter,mask);


%% Plot Rates as a Transect

% Plot settings
plot_ages = 1;  % Also plot exposure ages?  yes = 1, no = 0
save_plot = 0;  % Save plot?  1 to save as .png and .eps, otherwise 0
time_lim = [];  % Optionally set limits of time axis (in ka), otherwise leave empty
pos_lim = [];   % Optionally set limits of relative position axis (in m or km), otherwise leave empty

% Plot
plot_transect_linear_regress(regressed_rates,ages_ka,transect_type,plot_ages,save_plot,mask,time_lim,pos_lim);

