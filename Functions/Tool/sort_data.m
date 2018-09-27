%
% sorted_data = sort_data(in_data,N,names)
% 
% Sorts in_data if multiple samples were combined for nuclide measurements.
% This is determined if nuclide measurements (mean and uncertainty) are
% identical for different samples.
%
% If such duplicate measurements are found then a new data matrix is
% produced with each new sample representing a unique nuclide measurement.
% Alternatively, if no duplicates are found, the new matrix is the same as
% the old matrix.
%
% Output is to be used with get_data_complex.m. It is a struct of the 
% sample matrix, sample names, logical of which samples had duplicate 
% measurements, and the other details for these combined samples.
%
% Written by Richard Selwyn Jones, Durham University
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite, which is built on versions of 
% CRONUS-Earth and CRONUScalc code.
%
%
%%

function sorted_data = sort_data(in_data,N,names)

  % Find duplicate nuclide measurements
  [~,N_ia,N_ic] = unique(N(:,1),'stable');
  [~,dup_N] = ismember(N(:,1), N(N_ia(accumarray(N_ic,1)>1),1));
  
  [~,N_ia,N_ic] = unique(N(:,2),'stable');
  [~,dup_dN] = ismember(N(:,2), N(N_ia(accumarray(N_ic,1)>1),2));
  
  % Check if any duplicates are consistent for both mean and 1 sig. measurements
  true_dup_N = dup_N & dup_dN;
  num_dup_N = dup_N(true_dup_N);
  
  
  % Sort only if duplicates are found
  if any(true_dup_N) 
      
      in_data_N = in_data(~true_dup_N,:); % Get non-duplicate data
      in_data_names = names(~true_dup_N)'; % Get non-duplicate names
      uni_dup = unique(num_dup_N); % Find number of duplicates
      sorted_names = cell(max(num_dup_N),1);
      sorted = zeros(max(num_dup_N),length(in_data(1,:)));
      
      
      for a = 1:max(num_dup_N)
          
          num = uni_dup(a);
          ind = find(dup_N == num);
          
          name = names(ind(1));
          
          sorted_N = in_data(ind(1),1:5);
          sorted_N = [sorted_N mean(in_data(ind,6:7))];
          sorted_N = [sorted_N in_data(ind(1),8:12)];
          sorted_N = [sorted_N mean(in_data(ind,13:15))];
          sorted_N = [sorted_N in_data(ind(1),16)];
          
          sorted_names(a) = name;
          sorted(a,:) = sorted_N(1,:);
          
          % Retain certain individual sample values
          sorted_data.dup_names{a} = names(ind);
          sorted_data.dup_top_z{a} = in_data(ind,13)';
          sorted_data.dup_bottom_z{a} = in_data(ind,14)';
          sorted_data.dup_weight{a} = in_data(ind,15)';
          
      end
      
      sorted_data.names = [in_data_names; sorted_names]; % Append sorted names
      sorted_data.data = [in_data_N; sorted]; % Append sorted data
      sorted_data.duplicates = [zeros(1,length(in_data_N(:,1))) ones(1,length(uni_dup))] > 0 ; % Identify duplicates
  
  % Otherwise just take in_data for nuclide
  else
      
      sorted_data.names = names(~true_dup_N)';
      sorted_data.data = in_data(~true_dup_N,:);
      sorted_data.duplicates = [];
      
  end
  
  
end
