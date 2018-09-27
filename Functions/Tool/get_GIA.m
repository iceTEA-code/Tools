%
% elev_data = get_GIA(sample_data,GIA_model,time_arr)
%
% Gets elevation change through time, as derived from a glacio-isostatic 
% adjustment (GIA) model. The GIA model data was produced by P. Whitehouse,
% using a three-layer approximation of the VM2 Earth model for both the
% ICE-5G and ICE-6G global ice models (Peltier, 2004; Peltier et al., 
% 2015). These ice models have time steps from 122 ka and 26 ka before 
% present, respectively - the earliest model value is used prior to these 
% times. The computed vertical deformation is then calculated relative to 
% time-varying sea level using a sea-level model that includes global mean
% sea level, rotational feedbacks and gravitational effects. The relative 
% elevation is determined for each sample site by interpolating the model 
% data in space and time.
%
% sample_data is a struct containing necessary sample information produced 
% using get_data.m.
%
% GIA_model should be either "I5G" or "I6G", corresponding to the ice model
% used to calculate GIA (i.e. ICE-5G or ICE-6G).
%
% time_arr is an array of time (present to past), created in elev_correct.m.
%
% The output struct contains the time array, type of GIA model, and present
% elevation and array of elevation differences (relative to present) for 
% each sample. These are to be used with elev_correct.m, which determines 
% the absolute elevation of each sample through time.
%
% Written by Richard Selwyn Jones, Durham University
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite, which is built on versions of 
% CRONUS-Earth and CRONUScalc code.
%
%
%%

function elev_data = get_GIA(sample_data,GIA_model,time_arr)
  
  % Load GIA-modelled elevation data (values are difference from present, in metres)
  if strcmpi(GIA_model,'I5G')
      load('i5g.mat');
      GIA_data = i5g;
  elseif strcmpi(GIA_model,'I6G')
      load('i6g.mat');
      GIA_data = i6g;
  else
      error('the GIA_model for get_GIA should be "I5G" or "I6G"!');
  end
  lat_grid = GIA_data.lat_grid;
  lon_grid = GIA_data.lon_grid;
  GIA_time = GIA_data.time;
  elev_diff = GIA_data.elev_diff_relpastSLmodel;
  %elev_diff = GIA_data.elev_diff; % Uncomment to instead look at effects without accounting for past SL change (vertical deformation only)
  
  % Get elevations for each sample
  elevdiff_arr = cell(1,numel(sample_data.s));
  for a = 1:numel(sample_data.s)
      
      % Get sample details
      s_lat = sample_data.CC(a,1); % Latitude
      s_lon = sample_data.CC(a,2); % Longitude
      s_lon_360 = wrapTo360(s_lon); % Convert longitude from -180-180 (180W-180E) to 0-360
      s_elev = sample_data.CC(a,3); % Modern sample elevation
      
      % Get elevations for each timestep
      s_elev_diff = zeros(1,length(elev_diff));
      lat_arr = reshape(lat_grid,[],1); % Create latitude array
      lon_arr = reshape(lon_grid,[],1); % Create longitude array     
      
      for b = 1:length(elev_diff)
          this_elevdiff_arr = reshape(elev_diff{b},[],1); % Create elevation difference array
          interpolant = scatteredInterpolant(lat_arr,lon_arr,this_elevdiff_arr,'natural'); % Make interpolant
          this_elev_diff = interpolant(s_lat,s_lon_360); % Interpolate for sample site
          s_elev_diff(1,b) = this_elev_diff; % Record elevation relative to present
      end
      
      % Interpolate in time
      GIA_time_flip = fliplr(GIA_time); % Flip time - youngest to oldest
      s_elev_diff_flip = fliplr(s_elev_diff); % Flip elevation array
      s_elev_diff_full = interp1(GIA_time_flip,s_elev_diff_flip,time_arr); % Interpolate for full time array
      
      % Extrapolate oldest elevation value back in time
      elev_nan_log = isnan(s_elev_diff_full);
      oldest_ind = find(~elev_nan_log,1,'last');
      oldest_val = s_elev_diff_full(oldest_ind);
      s_elev_diff_full(elev_nan_log) = oldest_val;
      
      elevdiff_arr{a} = s_elev_diff_full;
      sample_elevs(a) = s_elev;
      
  end
  
  
  % Export
  elev_data.present_elev = sample_elevs;
  elev_data.GIA_model = GIA_model;
  elev_data.time_arr = time_arr;
  elev_data.elev_change = elevdiff_arr;
 
end
