%
% sample_data = norm_concs(sample_data)
%
% Normalises nuclide concentrations of given samples.
%
% sample_data is a required struct containing necessary sample details, 
% initially created using get_data.m.
%
% Output is the input with the addition of production rates and normalised
% concentrations for each sample and nuclide.
%
% Written by Richard Selwyn Jones, Durham University.
% Modified from code by Greg Balco, Berkeley Geochronology Center, that was
% published in Schaefer et al., Nature, 2016.
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite.
%
%
%%

function sample_data = norm_concs(sample_data)

  % Check inputs
  if (nargin ~= 1)
      error('norm_concs has wrong number of inputs!');
  end
  if isempty(sample_data.cover.z)
      sample_data.cover.z = 0;
  end
  

  s = sample_data.s;
  
  for a = 1:length(s)
      
      % Compute depths
      this_top_z = s{a}.top_z_gcm2 + sample_data.cover.z;
      this_bottom_z = s{a}.bottom_z_gcm2 + sample_data.cover.z;
      
      this_weight = s{a}.weight;
      
      % Get sample parameters
      sf = sample_data.sf1026{a};
      cp = sample_data.cp1026{a};
      
      
      if s{a}.nuclide10 == 1
          
          % Define production rate function
          pr_func = @(z,t) PR_Z(z,sample_data.pp,sf,cp,10);
          
          % Integrate over depth for each sample
          Integrated_PR = zeros(1,length(this_top_z));
          for b = 1:length(this_top_z)
              % Integral of production rate in interval divided by interval thickness
              Integrated_PR(b) = integral(pr_func,this_top_z(b),this_bottom_z(b),'RelTol',1e-3,'AbsTol',1e-3) ./ (this_bottom_z(b) - this_top_z(b));
          end
          
          % Average by mineral weight
          sample_data.s{a}.PR_10 =  (sum(Integrated_PR.*this_weight)) ./ sum(this_weight);
          
          % Normalise nuclide concentration
          sample_data.s{a}.norm_N10 = s{a}.N10 ./ sample_data.s{a}.PR_10;
          sample_data.s{a}.norm_dN10 = s{a}.dN10 ./ sample_data.s{a}.PR_10;
      
      end
          
      if s{a}.nuclide26 == 1
          
          % Define production rate function
          pr_func = @(z,t) PR_Z(z,sample_data.pp,sf,cp,26);
          
          % Integrate over depth for each sample
          Integrated_PR = zeros(1,length(this_top_z));
          for b = 1:length(this_top_z)
              % Integral of production rate in interval divided by interval thickness.
              Integrated_PR(b) = integral(pr_func,this_top_z(b),this_bottom_z(b),'RelTol',1e-3,'AbsTol',1e-3) ./ (this_bottom_z(b) - this_top_z(b));
          end
          
          % Average by mineral weight
          sample_data.s{a}.PR_26 =  (sum(Integrated_PR.*this_weight)) ./ sum(this_weight);
          
          % Normalise nuclide concentration
          sample_data.s{a}.norm_N26 = s{a}.N26 ./ sample_data.s{a}.PR_26;
          sample_data.s{a}.norm_dN26 = s{a}.dN26 ./ sample_data.s{a}.PR_26;
          
      end
      
  end
  
end
