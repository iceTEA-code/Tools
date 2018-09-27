%
% corr_ages_ka = cov_correct(sample_data,scaling_model,cover_type,cover_depth)
% corr_ages_ka = cov_correct(sample_data,scaling_model,cover_type,cover_depth,cover_density)
% 
% Calculates exposure ages corrected for surface cover using the samples in 
% sample_data, dependent on the scaling_model set and cover inputs. Either 
% a preset cover type can be selected or a manual cover density can be 
% used. It uses calculations in CRONUScalc (Marrero et al., 2016).
%
% sample_data is a required struct, created using get_data.m.
%
% scaling_model should be one of 'DE','DU','LI','ST','LM','LSD'/'SF' or 
% 'LSDn'/'SA'.
%
% cover_type should be 'snow', 'freshwater', 'seawater', 'loess', 'till', 
% 'soil', 'ash', or 'manual'.
%
% cover_depth is a required value of the depth of surface cover (in cm).
%
% cover_density is an optional input that is to be used if the cover_type
% is set to "manual". A value should be in g/cm^3.
%
% Output is age for each nuclide, sample names, surface cover and total 
% shielding factors, scaling model used, logical of nuclides present, and 
% production rate with time if plot_prod is set.
%
% Written by Richard Selwyn Jones, Durham University
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite, which is built on versions of 
% CRONUS-Earth and CRONUScalc code.
%
%
%%

function corr_ages_ka = cov_correct(sample_data,scaling_model,cover_type,cover_depth,cover_density)
  
  % Check inputs
  if (nargin < 4 || nargin > 5)
      error('cov_correct has wrong number of inputs!');
  end
  if strcmpi(scaling_model,'LSD')
      scaling = 'SF';
  elseif strcmpi(scaling_model,'LSDn')
      scaling = 'SA';
  else
      scaling = scaling_model;
  end
  if nargin == 4 || ~strcmpi(cover_type,'manual')
      cover_density = [];
  end

  % Determine what nuclides are present
  NN = [any(sample_data.logical_10) any(sample_data.logical_26)];
  
  
  % Compute surface cover shielding factors
  for i = 1:length(sample_data.CC(:,1))
      cover_shield_facs(i,1) = get_cov_shield(cover_type,cover_depth,sample_data.CC(i,13),cover_density);
      sample_data.CC(i,7) = sample_data.CC(i,7) .* cover_shield_facs(i); % Combine with topographic shielding
  
      % Print shielding factors
      name = cellstr(sample_data.s{i}.name);
      disp(['Sample ' name{1} ':  Cover shielding factor ' sprintf('%0.4f',cover_shield_facs(i,1)) '  (total shielding ' sprintf('%0.4f',sample_data.CC(i,7)) ')']);
  end
  
  % Calculate exposure ages
  disp('Calculating cover-corrected exposure ages...')
  
  if NN(1)
      corr_ages_ka.Be10 = NaN(length(sample_data.CC(:,1)),3);
      for a = 1:nnz(sample_data.logical_10)
          [output,times,plotprod,~] = be10age_mod(sample_data.CC(a,:),sample_data.CC_uncert(a,:),scaling); % Full calculation with uncertainties
          age(a,:) = output([1,2,12],:)';
          prod_time{a} = [(0-times)',plotprod']; % Save production information (converting time to ka BP)
      end
      corr_ages_ka.Be10(sample_data.logical_10,:) = age(sample_data.logical_10,:);
      corr_ages_ka.prod_time_10 = prod_time;
      corr_ages_ka.logical_10 = sample_data.logical_10;
  end

  if NN(2)
      corr_ages_ka.Al26 = NaN(length(sample_data.CC(:,1)),3);
      for a = 1:nnz(sample_data.logical_26)
          [output,times,plotprod,~] = al26age_mod(sample_data.CC(a,:),sample_data.CC_uncert(a,:),scaling); % Full calculation with uncertainties
          age(a,:) = output([1,2,12],:)';
          prod_time{a} = [(0-times)',plotprod'];
      end
      corr_ages_ka.Al26(sample_data.logical_26,:) = age(sample_data.logical_26,:);
      corr_ages_ka.prod_time_26 = prod_time;
      corr_ages_ka.logical_26 = sample_data.logical_26;
  end
  
  
  % Get sample names
  for b = 1:length(sample_data.s)
      sample_names(b) = sample_data.s{b}.name;
  end
  
  % Append sample names, cover and total shielding factors, elevations and positions, as well as scaling model and logical of nuclides
  corr_ages_ka.sample_names = sample_names;
  corr_ages_ka.cover_shielding = cover_shield_facs;
  corr_ages_ka.total_shielding = sample_data.CC(:,7);
  corr_ages_ka.elevation = sample_data.CC(:,3)';
  corr_ages_ka.position = sample_data.position;
  corr_ages_ka.scaling_model = scaling_model;
  corr_ages_ka.NN = NN;
  
  disp('done.')
  
end
