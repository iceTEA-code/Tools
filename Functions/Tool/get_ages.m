%
% ages_ka = get_ages(sample_data,ages_name)
% 
% Gets exposure ages that have been calculated elsewhere, as entered for 
% samples in the input data, as well as associated data.
%
% sample_data is a required struct, created using get_data.m.
%
% ages_name is a character input for the name of the exposure age dataset.
%
% Output is age for each nuclide, sample names, scaling model used, logical
% of nuclides present, and elevations and relative positions of samples.
%
% Written by Richard Selwyn Jones, Durham University
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite.
%
%
%%

function ages_ka = get_ages(sample_data,ages_name)

  % Check inputs
  if (nargin ~= 2)
      error('get_ages has wrong number of inputs!');
  end

  % Determine what nuclides are present
  NN = [any(sample_data.logical_10) any(sample_data.logical_26)];
  
  
  try
  
      if NN(1)
          ages_ka.Be10 = NaN(length(sample_data.CC(:,1)),3);
          ages = sample_data.ages.Be10 / 1000; % Convert to ka
          ages_ka.Be10(sample_data.logical_10,:) = ages(sample_data.logical_10,:);
          ages_ka.logical_10 = sample_data.logical_10;
      end
      
      if NN(2)
          ages_ka.Al26 = NaN(length(sample_data.CC(:,1)),3);
          ages = sample_data.ages.Al26 / 1000; % Convert to ka
          ages_ka.Al26(sample_data.logical_26,:) = ages(sample_data.logical_26,:);
          ages_ka.logical_26 = sample_data.logical_26;
      end
      
      
      % Get sample names
      for b = 1:length(sample_data.s)
          sample_names(b) = sample_data.s{b}.name;
      end
      
      % Append sample names, elevations and positions, as well as scaling model, logical of nuclides and name of ages dataset
      ages_ka.sample_names = sample_names;
      ages_ka.elevation = sample_data.CC(:,3)';
      ages_ka.elevation_err = sample_data.CC_uncert(:,3)';
      ages_ka.position = sample_data.position;
      ages_ka.scaling_model = sample_data.ages.scaling_model;
      ages_ka.NN = NN;
      ages_ka.ages_name = ages_name;
  
  
  catch
      error('the input data (.xlsx) does not contain exposure ages!');
  end
    
  
end
