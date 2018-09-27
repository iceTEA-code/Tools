%
% ages_ka = age_calc(sample_data,scaling_model)
% ages_ka = age_calc(sample_data,scaling_model,plot_prod)
% 
% Calculates exposure ages for the samples in sample_data, dependent on the
% scaling_model set. It uses calculations in CRONUScalc (Marrero et al., 
% 2016).
%
% sample_data is a required struct, created using get_data.m.
%
% scaling_model should be one of 'DE','DU','LI','ST','LM','LSD'/'SF' or 
% 'LSDn'/'SA'.
%
% plot_prod is an optional binary input (1 or 0). If set as 1, the 
% production rate through time used to calculate the age of each sample is 
% outputted. Default is 0.
%
% Output is age for each nuclide, sample names, scaling model used, logical
% of nuclides present, and production rate with time if plot_prod is set.
%
% Written by Richard Selwyn Jones, Durham University
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite, which is built on versions of 
% CRONUS-Earth and CRONUScalc code.
%
%
%%

function ages_ka = age_calc(sample_data,scaling_model,plot_prod)
  
  % Check inputs
  if (nargin < 2 || nargin > 3)
      error('age_calc has wrong number of inputs!');
  end
  if (nargin < 3 || isempty(plot_prod))
      plot_prod = 0;
  end
  if strcmpi(scaling_model,'LSD')
      scaling = 'SF';
  elseif strcmpi(scaling_model,'LSDn')
      scaling = 'SA';
  else
      scaling = scaling_model;
  end

  % Determine what nuclides are present
  NN = [any(sample_data.logical_10) any(sample_data.logical_26)];
  
  
  % Calculate exposure ages
  disp('Calculating exposure ages...')
  if NN(1)
      ages_ka.Be10 = NaN(length(sample_data.CC(:,1)),3);
      for a = 1:nnz(sample_data.logical_10)
          [output,times,plotprod,~] = be10age_mod(sample_data.CC(a,:),sample_data.CC_uncert(a,:),scaling); % Full calculation with uncertainties
          age(a,:) = output([1,2,12],:)';
          prod_time{a} = [(0-times)',plotprod']; % Save production information (converting time to ka BP)
      end
      ages_ka.Be10(sample_data.logical_10,:) = age(sample_data.logical_10,:);
      if plot_prod == 1
          ages_ka.prod_time_10 = prod_time;
      end
      ages_ka.logical_10 = sample_data.logical_10;
  end

  if NN(2)
      ages_ka.Al26 = NaN(length(sample_data.CC(:,1)),3);
      for a = 1:nnz(sample_data.logical_26)
          [output,times,plotprod,~] = al26age_mod(sample_data.CC(a,:),sample_data.CC_uncert(a,:),scaling); % Full calculation with uncertainties
          age(a,:) = output([1,2,12],:)';
          prod_time{a} = [(0-times)',plotprod'];
      end
      ages_ka.Al26(sample_data.logical_26,:) = age(sample_data.logical_26,:);
      if plot_prod == 1
          ages_ka.prod_time_26 = prod_time;
      end
      ages_ka.logical_26 = sample_data.logical_26;
  end
  
  
  % Get sample names
  for b = 1:length(sample_data.s)
      sample_names(b) = sample_data.s{b}.name;
  end
  
  % Append sample names, elevations and positions, as well as scaling model and logical of nuclides
  ages_ka.sample_names = sample_names;
  ages_ka.elevation = sample_data.CC(:,3)';
  ages_ka.elevation_err = sample_data.CC_uncert(:,3)';
  ages_ka.position = sample_data.position;
  ages_ka.scaling_model = scaling_model;
  ages_ka.NN = NN;
  
  disp('done.')
  
end
