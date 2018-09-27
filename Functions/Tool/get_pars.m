%
% sample_data = get_pars(sample_data,scaling_model)  
%
% Gets parameters for each sample using CRONUScalc functions:
% pp=physical constants
% sp=sample parameters
% sf=scaling factors
% cp=computed parameters
% Exports to sample_data struct.
%
% Requires sample data to be loaded and sorted with get_data.m, and a
% scaling model to be set (i.e. 'DE','DU','LI','ST','LM','LSD'/'SF',
% 'LSDn'/'SA').
%
% Written by Richard Selwyn Jones, Durham University
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite, which is built on versions of 
% CRONUS-Earth and CRONUScalc code.
%
%
%%

function sample_data = get_pars(sample_data,scaling_model)

  if strcmpi(scaling_model,'LSD')
      scaling = 'SF';
  elseif strcmpi(scaling_model,'LSDn')
      scaling = 'SA';
  else
      scaling = scaling_model;
  end
  
  % Find number of samples
  n_samples = numel(sample_data.CC(:,1));

  % Get the physical parameters
  sample_data.pp = physpars(scaling);

  
  % Get sample-specific parameters
  for s = 1:n_samples
  
      % Extract the sample parameters from a sample_data vector
      sample_data.sp1026{s} = samppars1026(sample_data.CC(s,:));
      
      % Get the scale factors
      sample_data.sf1026{s} = scalefacs1026(sample_data.sp1026{s});
      
      % Computed parameters
      maxdepth = 2500; % Maxdepth specifies the maximum depth (in g/cm^2) for which production rates will be computed
      sample_data.cp1026{s} = comppars1026(sample_data.pp,sample_data.sp1026{s},sample_data.sf1026{s},maxdepth);
      
      % Go ahead and produce contemporary scaling factors
      sample_data.sf1026{s}.currentsf = getcurrentsf(sample_data.sf1026{s},0,scaling,'albe');
      
  end
  
  % Append scaling model
  sample_data.scaling_model = scaling;
      
end
