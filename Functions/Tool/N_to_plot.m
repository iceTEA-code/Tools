%
% N_data = N_to_plot(sample_data)
% N_data = N_to_plot(sample_data,predN)
% 
% Sorts measured (sample_data) or predicted (predN) nuclide concentrations 
% for plotting. If predN is excluded then concentrations are sorted for 
% measured data only, if predN is included then only the predicted 
% concentrations are normalised and sorted.
%
% sample_data is a required struct, created using get_data.m. It must
% include the normalised concentrations, calculated using norm_concs.m.
%
% predN is data struct of predicted nuclide concentrations. The summed 
% concentrations derived from the complex history model are used if
% present. Otherwise it uses whatever the predN input is, but must be a
% struct containing at least nuclide variables (e.g. N10).
%
% In cases where sample material has multiple concentration measurements
% for one nuclide and only a single measurement for another nuclide (e.g. 
% 2x 10Be measurements but 1x 26Al measurement from 2x core segments), the 
% multiple concentrations (i.e. for 10Be in above example) are mixed. This 
% is determined from which samples have common depths.
%
% Output is used with plot_2iso_concs.m. It is a struct of normalised 
% nuclide concentrations for all samples, for each nuclide.
%
% Written by Richard Selwyn Jones, Durham University
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite.
%
%
%%

function N_data = N_to_plot(sample_data,predN)

  % Check inputs
  if (nargin > 2)
      error('N_to_plot has wrong number of inputs!');
  end
  
  try
      pred = predN.sum; % Use summed predicted nuclide concentrations
  catch
      if (nargin < 2)
          pred = [];
      else
          pred = predN; % If summed concentrations are missing, use whatever the input is
      end
  end

  
  % Isolate sample data for each nuclide
  if any(sample_data.logical_10)
      data_N10 = sample_data.s(sample_data.logical_10);
      % Get top and bottom depths for each sample measurement
      for a = 1:length(data_N10)
          top_z_10(a) = min(data_N10{a}.top_z);
          bottom_z_10(a) = max(data_N10{a}.bottom_z);
      end
  end
  if any(sample_data.logical_26)
      data_N26 = sample_data.s(sample_data.logical_26);
      % Get top and bottom depths for each sample measurement
      for a = 1:length(data_N26)
          top_z_26(a) = min(data_N26{a}.top_z);
          bottom_z_26(a) = max(data_N26{a}.bottom_z);
      end
  end
  
  
  % Find common depths between nuclides
  if any(sample_data.logical_10) && any(sample_data.logical_26)
      common_top_z_1026 = intersect(top_z_10,top_z_26);
      common_bottom_z_1026 = intersect(bottom_z_10,bottom_z_26);
  end
  
  
  % Find samples to mix for each nuclide
  if any(sample_data.logical_10)
      % Create logical matrix (common depths x sample measurements)
      logical_mix_10 = true(length(common_top_z_1026),length(top_z_10));
      % Get sample depths
      for b = 1:length(common_top_z_1026)
          for c = 1:length(top_z_10)
              % Find if sample depth is within range of common depths, for each sample measurement and common depth;
              % Save to matrix
              logical_mix_10(b,c) = (top_z_10(c) >= common_top_z_1026(b) & top_z_10(c) <= common_bottom_z_1026(b));
          end
      end
  end
  if any(sample_data.logical_26)
      % Create logical matrix (common depths x sample measurements)
      logical_mix_26 = true(length(common_top_z_1026),length(top_z_26));
      % Get sample depths
      for b = 1:length(common_top_z_1026)
          for c = 1:length(top_z_26)
              % Find if sample depth is within range of common depths, for each sample measurement and common depth;
              % Save to matrix
              logical_mix_26(b,c) = (top_z_26(c) >= common_top_z_1026(b) & top_z_26(c) <= common_bottom_z_1026(b));
          end
      end
  end
  
  
  % Sort data for plotting
  
  if isempty(pred) % Do for sample data
      
      if any(~all(logical_mix_10,2)) || any(~all(logical_mix_26,2)) % Mix data
          if any(sample_data.logical_10)
              for d = 1:length(common_top_z_1026)
                  log_10 = logical_mix_10(d,:);
                  for e = 1:length(data_N10)
                      norm_N10_all(e) = data_N10{e}.norm_N10;
                      norm_dN10_all(e) = data_N10{e}.norm_dN10;
                      sum_weights_10_all(e) = sum(data_N10{e}.weight);
                  end
                  norm_N10 = norm_N10_all(log_10);
                  norm_dN10 = norm_dN10_all(log_10);
                  sum_weights_10 = sum_weights_10_all(log_10);
                  % Combine concentrations
                  for f = 1:length(norm_N10)
                      norm_N10_weight(f) = norm_N10(f) .* sum_weights_10(f);
                      norm_dN10_weight_sq(f) = ((norm_dN10(f) .* sum_weights_10(f)) ./ sum(sum_weights_10))^2;
                  end
                  mixed_norm_N10 = sum(norm_N10_weight) ./ sum(sum_weights_10);
                  mixed_norm_dN10 = sqrt(sum(norm_dN10_weight_sq));
                  % Export mixed concentrations
                  N_data.norm_N10(d) = mixed_norm_N10;
                  N_data.norm_dN10(d) = mixed_norm_dN10;
              end
          end
          if any(sample_data.logical_26)
              for d = 1:length(common_top_z_1026)
                  log_26 = logical_mix_26(d,:);
                  for e = 1:length(data_N26)
                      norm_N26_all(e) = data_N26{e}.norm_N26;
                      norm_dN26_all(e) = data_N26{e}.norm_dN26;
                      sum_weights_26_all(e) = sum(data_N26{e}.weight);
                  end
                  norm_N26 = norm_N26_all(log_26);
                  norm_dN26 = norm_dN26_all(log_26);
                  sum_weights_26 = sum_weights_26_all(log_26);
                  % Combine concentrations
                  for f = 1:length(norm_N26)
                      norm_N26_weight(f) = norm_N26(f) .* sum_weights_26(f);
                      norm_dN26_weight_sq(f) = ((norm_dN26(f) .* sum_weights_26(f)) ./ sum(sum_weights_26))^2;
                  end
                  mixed_norm_N26 = sum(norm_N26_weight) ./ sum(sum_weights_26);
                  mixed_norm_dN26 = sqrt(sum(norm_dN26_weight_sq));
                  % Export mixed concentrations
                  N_data.norm_N26(d) = mixed_norm_N26;
                  N_data.norm_dN26(d) = mixed_norm_dN26;
              end
          end
          
      else % Don't mix
          if any(sample_data.logical_10) && any(sample_data.logical_26)
              for d = 1:length(sample_data.s)
                  N_data.norm_N10(d) = sample_data.s{d}.norm_N10;
                  N_data.norm_dN10(d) = sample_data.s{d}.norm_dN10;
                  N_data.norm_N26(d) = sample_data.s{d}.norm_N26;
                  N_data.norm_dN26(d) = sample_data.s{d}.norm_dN26;
              end
          end
      end
  
      
  else % Do for predicted nuclide concentrations
      
      if any(~all(logical_mix_10,2)) || any(~all(logical_mix_26,2)) % Mix data
          if any(sample_data.logical_10)
              pred_N10 = pred.N10(sample_data.logical_10);
              for d = 1:length(pred_N10)
                  % Normalise nuclide concentrations
                  norm_N10_all(d) = pred_N10(d) ./ data_N10{d}.PR_10;
              end
              for e = 1:length(common_top_z_1026)
                  log_10 = logical_mix_10(e,:);
                  for f = 1:length(data_N10)
                      sum_weights_10_all(f) = sum(data_N10{f}.weight);
                  end
                  norm_N10 = norm_N10_all(log_10);
                  sum_weights_10 = sum_weights_10_all(log_10);
                  % Combine concentrations
                  for h = 1:length(norm_N10)
                      norm_N10_weight(h) = norm_N10(h) .* sum_weights_10(h);
                  end
                  mixed_norm_N10 = sum(norm_N10_weight) ./ sum(sum_weights_10);
                  % Export mixed concentrations
                  N_data.norm_N10(e) = mixed_norm_N10;
              end
          end
          if any(sample_data.logical_26)
              pred_N26 = pred.N26(sample_data.logical_26);
              for d = 1:length(pred_N26)
                  % Normalise nuclide concentrations
                  norm_N26_all(d) = pred_N26(d) ./ data_N26{d}.PR_26;
              end
              for e = 1:length(common_top_z_1026)
                  log_26 = logical_mix_26(e,:);
                  for g = 1:length(data_N26)
                      sum_weights_26_all(g) = sum(data_N26{g}.weight);
                  end
                  norm_N26 = norm_N26_all(log_26);
                  sum_weights_26 = sum_weights_26_all(log_26);
                  % Combine concentrations
                  for i = 1:length(norm_N26)
                      norm_N26_weight(i) = norm_N26(i) .* sum_weights_26(i);
                  end
                  mixed_norm_N26 = sum(norm_N26_weight) ./ sum(sum_weights_26);
                  % Export mixed concentrations
                  N_data.norm_N26(e) = mixed_norm_N26;
              end
          end
          
      else % Don't mix
          % Normalise nuclide concentrations
          if any(sample_data.logical_10)
              for d = 1:length(pred.N10)
                  N_data.norm_N10(d) = pred.N10(d) ./ sample_data.s{d}.PR_10;
              end
          end
          if any(sample_data.logical_26)
              for d = 1:length(pred.N26)
                  N_data.norm_N26(d) = pred.N26(d) ./ sample_data.s{d}.PR_26;
              end
          end
      end
  end
    
end
