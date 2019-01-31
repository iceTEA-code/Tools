%
% sample_data = get_data(input_name)
% 
% Loads input data (either from an .xlsx file, or headerless .txt or .csv), 
% extracts sample names and concentrations, sorts data for functions, and 
% exports to single output.
%
% For CRONUScalc calculations, the erosion rate and inheritance is set to 
% zero. Attenuation lengths are automatically determined using the Sato 
% model, based on the location information of each sample.
%
% Input data must include the following information:
% 1. Sample name
% 2. Latitude (decimal degrees)
% 3. Longitude (decimal degrees)
% 4. Elevation (m asl)
% 5. Pressure (hPa) (zero if not known)
% 6. Elevation uncertainty (m asl) (zero if not known)
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
%
% Written by Richard Selwyn Jones, Durham University
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite, which is built on versions of 
% CRONUS-Earth and CRONUScalc code.
%
%
%%

function sample_data = get_data(input_name)
  
  if ispc
      input_path = strcat(pwd,'\',input_name);
  else
      input_path = strcat(pwd,'/',input_name);
  end
  [~,~,ext] = fileparts(input_path);
  if isempty(ext)
      error('File name of sample data does not include the file type extension (e.g. .xlsx)')
  end

  % Load input data
  if strcmpi(ext,'.xlsx') || strcmpi(ext,'.xls')
      [in_data,in_txt,in_raw] = xlsread(input_name);
      for a = 1:length(in_raw(1,:))
          this_header = in_raw{1,a};
          strrow(a) = ~isnumeric(this_header);
      end
      n_raw_columns = numel(in_raw(1,:));
      if (n_raw_columns > 15) && all(~strrow(16:end)) % Remove any empty columns
          in_raw = in_raw(:,1:15);
      elseif (n_raw_columns > 22) && all(~strrow(23:end))
          in_raw = in_raw(:,1:22);
      end
      if all(strrow) % Remove header
          in_raw = in_raw(2:end,:); 
          in_txt = in_txt(2:end,:);
      end
  elseif strcmpi(ext,'.txt') || strcmpi(ext,'.csv')
      warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')
      in_raw = table2cell(readtable(input_name));
      in_txt = cell(size(in_raw));
      in_data = zeros(size(in_raw));
      for a = 1:length(in_raw(:,1))
          for b = 1:length(in_raw(1,:))
              this_raw = in_raw{a,b};
              this_num_log = isnumeric(this_raw);
              if this_num_log
                  in_data(a,b) = this_raw;
              else
                  in_txt{a,b} = this_raw;
              end
          end
      end
      in_data = in_data(:,2:end);
  end
  
  n_raw_columns = numel(in_raw(1,:));
  n_data_columns = numel(in_data(1,:));
  
  % Check inputs
  if (n_raw_columns > 22)
      error('sample input data has too many fields!');
  end
  if (n_raw_columns > 15 && n_raw_columns < 22)
      error('sample input data seems to have exposure ages but fields are missing!');
  end
  
  % Check sample names
  col1_nan_log = isnan(in_data(:,1)); % Check for NaNs in column 1
  if any(col1_nan_log) && (n_data_columns == 15 || n_data_columns == 22)
      miss_names = cellstr(num2str(in_data(~col1_nan_log,1))); % Get and convert numeric names from in_data
      in_txt_names = in_raw(:,1);
      in_txt_names(~col1_nan_log,1) = miss_names; % Add missed non-numeric names to data
      in_txt = cell(size(in_raw));
      in_txt(1,:) = in_raw(1,:); in_txt(:,1) = in_txt_names; % Combine for new text data
      in_data = in_data(:,2:end); % Adjust in_data (remove names column)
      warning('some sample names are numeric... fixed.');
  end
  
  if (n_raw_columns > 15)
      ages = in_data(:,15:20);
      scaling = in_txt(:,end);
  else
      ages = [];
      scaling = [];
  end
  
  % Set data columns
  name_txt_c = 1;
  lat_c = 1;
  lon_c = 2;
  elev_c = 3;
  pressure_c = 4;
  elev_err_c = 5;
  position_c = 6;
  thk_c = 7;
  density_c = 8;
  shield_c = 9;
  Be_c = 10;
  Be_sig_c = 11;
  Al_c = 12;
  Al_sig_c = 13;
  year_c = 14;
  Be_age_c = 15;
  Al_age_c = 18;
  
  
  % Get number of samples (and samples with 10-Be and 26-Al concs.)
  n_samples = numel(in_data(:,lat_c));
  
  % Get sample names
  names = in_txt(:,name_txt_c)';
  
  
  % If concentrations are NaN, make zero
  for a = 1:n_samples
      if (isnan(in_data(a,Be_c)))
          in_data(a,Be_c) = 0;
          in_data(a,Be_sig_c) = 0;
      end
      if (isnan(in_data(a,11)))
          in_data(a,Al_c) = 0;
          in_data(a,Al_sig_c) = 0;
      end
  end
 
  % Sort details for each sample and nuclide (currently only 10Be and 26Al)
  logical_10 = any(in_data(:,Be_c),2)';
  logical_26 = any(in_data(:,Al_c),2)';
  if ~any(logical_10) && ~any(logical_26)
      logical_10 = any(in_data(:,Be_age_c),2)';
      logical_26 = any(in_data(:,Al_age_c),2)';
  end
    
  NN = [any(logical_10) any(logical_26)];
  
  % Determine if sample details should be combined for each nuclide measurement, and sort if necessary
  if NN(1)
      in_data_10.data = in_data(logical_10,:);
      in_data_10.names = names(logical_10)';
  else
      in_data_10.data = [];
  end
  if NN(2)
      in_data_26.data = in_data(logical_26,:);
      in_data_26.names = names(logical_26)';
  else
      in_data_26.data = [];
  end
    
  logical_1026 = any(logical_10' & logical_26'); % Determine whether samples have both nuclide measurements
  
  % Get data
  if NN(1) && ~NN(2)
      comb_names = in_data_10.names;
      comb_data = in_data_10.data;
  elseif ~NN(1) && NN(2)
      comb_names = in_data_26.names;
      comb_data = in_data_26.data;
  elseif logical_1026
      comb_names = [in_data_10.names; in_data_26.names];
      comb_data = [in_data_10.data; in_data_26.data];
  end
  
  % Calculate number of samples
  n_comb = numel(comb_data(:,1));
  logical_comb_10 = any(comb_data(:,Be_c),2)';
  logical_comb_26 = any(comb_data(:,Al_c),2)';
  if ~any(logical_comb_10) && ~any(logical_comb_26)
      logical_comb_10 = any(comb_data(:,Be_age_c),2)';
      logical_comb_26 = any(comb_data(:,Al_age_c),2)';
  end
  
  
  % Collate sample details
  for c = 1:n_comb
      
      sample_data.s{c}.name = comb_names(c);
      
      if logical_comb_10(c)
          sample_data.s{c}.nuclide10 = 1;
          sample_data.s{c}.N10 = comb_data(c,Be_c);
          sample_data.s{c}.dN10 = comb_data(c,Be_sig_c);  
      else
          sample_data.s{c}.nuclide10 = 0;
      end
      if logical_comb_26(c)
          sample_data.s{c}.nuclide26 = 1;
          sample_data.s{c}.N26 = comb_data(c,Al_c);
          sample_data.s{c}.dN26 = comb_data(c,Al_sig_c);
      else
          sample_data.s{c}.nuclide26 = 0;
      end
           
  end
   
  
  % Assume certain sample details are unknown or zero
  e_rate = zeros(n_comb,1); % Erosion rate (mm/kyr)
  inh10 = zeros(n_comb,1);  % 10Be inheritance (atoms/g)
  inh26 = zeros(n_comb,1);  % 26Al inheritance (atoms/g)
  top_depth_gcm2 = zeros(n_comb,1); % Top depth (g/cm^2) - i.e. the surface
  
  % Initially set attenuation length to zero
  init_L = zeros(n_comb,1);
  
  % Combine necessary data for CronusCalc
  sample_data.CC = [comb_data(:,[lat_c:pressure_c,thk_c:shield_c]),e_rate,comb_data(:,[Be_c,Al_c]),inh10,inh26,init_L,top_depth_gcm2,comb_data(:,year_c)];
  
  % Determine atmospheric pressure (if unknown)
  sample_data.CC = atm_pressure(comb_data(:,lat_c),comb_data(:,lon_c),sample_data.CC);
  
%   % Determine attenuation lengths (g/cm^2) based on latitude (60 deg.)
%   if (sample_data.CC(:,1) > 60) | (sample_data.CC(:,1) < -60)
%       L = init_L + 140;  % Polar
%   else
%       L = init_L + 160;  % Non-polar
%   end

  % Determine attenuation lengths (g/cm^2) based on Sato model
  for d = 1:length(sample_data.CC(:,1))
      L(d,1) = attenuationlength(sample_data.CC(d,1),sample_data.CC(d,2),sample_data.CC(d,3),sample_data.CC(d,4));
  end
  sample_data.CC(:,13) = L; % Add attenuation lengths
  
  
  % Create uncertainties struct
  
  % Create a standard elevation uncertainty if not provided
  for eu = 1:length(comb_data(:,elev_err_c))
      if comb_data(eu,elev_err_c) == 0 || isnan(comb_data(eu,elev_err_c))
          elev_uncert(eu) = 5; % 5 metres
      else
          elev_uncert(eu) = comb_data(eu,elev_err_c);
      end
  end
  uncert = zeros(size(sample_data.CC)); % Assume zero uncertainty for now
  uncert(:,[9,10]) = comb_data(:,[Be_sig_c,Al_sig_c]); % Add nuclide concentration uncertainties
  uncert(:,3) = elev_uncert; % Add elevation uncertainties
  
  % Determine pressure uncertainties from elevation
  upelev_uncert = uncert;  upelev_uncert(:,3) = sample_data.CC(:,elev_c) + elev_uncert';
  upelev_uncert = atm_pressure(comb_data(:,lat_c),comb_data(:,lon_c),upelev_uncert);
  uppress_diff = sample_data.CC(:,4) - upelev_uncert(:,4);
  loelev_uncert = uncert;  loelev_uncert(:,3) = sample_data.CC(:,elev_c) - elev_uncert';
  loelev_uncert = atm_pressure(comb_data(:,lat_c),comb_data(:,lon_c),loelev_uncert);
  lopress_diff = sample_data.CC(:,4) - loelev_uncert(:,4);
  press_uncert = mean([abs(uppress_diff),abs(lopress_diff)],2); % Get average of upper and lower pressure uncertainties
  uncert(:,4) = press_uncert;
  
  sample_data.CC_uncert = uncert;
  
  
  % Create ages struct
  if ~isempty(ages)
      if NN(1)
          sample_data.ages.Be10 = ages(:,1:3);
      end
      if NN(2)
          sample_data.ages.Al26 = ages(:,4:6);
      end
      sample_data.ages.scaling_model = scaling{1}; % Assume the scaling used is the same for all samples
  end
  
  % Export logical of nuclides
  sample_data.logical_10 = logical_comb_10;
  sample_data.logical_26 = logical_comb_26;
  
  % Export sample position
  sample_data.position = comb_data(:,position_c)';
  
  % Create default cover depth
  sample_data.cover.z = 0;
  
  
end
