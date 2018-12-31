%% CALCULATE TOPOGRAPHIC SHIELDING FOR SAMPLE
%
% Calculates topographic shielding factor for cosmogenic nuclide production
% based on given horizon and dip values.
%
% Written by Richard Selwyn Jones, Durham University
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite, which is built on versions of 
% CRONUS-Earth and CRONUScalc code.
%
%
%%

clear % Start fresh
addpath(genpath(pwd));


% Specifiy shielding input
azimuths = [192 222 272 322 22 102 142]; % Measured azimuths for the horizon (0-360°)
horizon_angles = [3 3 6 8 12 8 0]; % Measured angle to the horizon for each azimuth (0-90°)

dip_direction = []; % Optional - dip direction of a rock surface (0-360°)
dip_angle = []; % Optional - dip angle of the surface (°)


% Calculate shielding factor
topo_out = get_topo(azimuths,horizon_angles,dip_direction,dip_angle);

