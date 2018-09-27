%
% new_ages_ka = find_outliers(ages_ka,feature)
% new_ages_ka = find_outliers(ages_ka,feature,sig)
%
% Finds outliers in an exposure age dataset using a generalised extreme 
% Studentized deviate test, which iteratively tests one mean exposure age 
% at a time.
%
% ages_ka is a required struct, containing ages calculated using 
% age_calc.m, sample names and a logical of the nuclides that were 
% measured.
%
% feature should be a binary input of whether these ages come from a single
% feature (1) or not (0). Outliers are only determined if the ages come
% from a feature.
%
% sig is an optional input to specify the significance level for
% determining outliers. Default is 0.05.
%
% Output is a new struct of exposure age data, additionally identifying 
% any outliers and replacing exposure ages (means and uncertainties) with 
% NaNs for those outliers.
%
% Written by Richard Selwyn Jones, Durham University
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite.
%
%
%%

function new_ages_ka = find_outliers(ages_ka,feature,sig)

  % Check inputs
  if (nargin < 2 || nargin > 3)
      error('find_outliers has wrong number of inputs!');
  end
  if (nargin < 3 || isempty(sig))
      sig = 0.05;
  end
  
  
  if feature ~= 1
      error('Outliers can only be looked for if the dataset is from a feature!');
  else
      
      % Create output from input
      new_ages_ka = ages_ka;
      
      
      % Find and remove outliers for each nuclide
      
      if ages_ka.NN(1)
          sample_names = ages_ka.sample_names(ages_ka.logical_10);
          ages = ages_ka.Be10;
          ages_m = ages(:,1)';
          
          max_outliers = numel(ages_m)-1; % Maximum number of outliers returned from the gESD test
          idx = gESDtest(ages_m,sig,max_outliers);
          outliers = ages_m(idx);
          
          ages(idx,:) = NaN; % Remove outliers from existing exposure age data
          
          new_ages_ka.Be10_outliers = outliers; % Export outliers
          new_ages_ka.Be10 = ages; % Export ages without outliers
          
          if ~isempty(outliers)
              disp('Outliers (Be-10):');
              for a = 1:length(outliers)
                  name = cellstr(sample_names(a));
                  disp(['sample ' name{1} ' (mean of ' sprintf('%0.2f',outliers(a)) ' ka)']);
              end
          else
              disp('Generalised ESD test found no outliers in Be-10 data');
          end
      end
      
      
      if ages_ka.NN(2)
          sample_names = ages_ka.sample_names(ages_ka.logical_26);
          ages = ages_ka.Al26;
          ages_m = ages(:,1)';
          
          max_outliers = numel(ages_m)-1; % Maximum number of outliers returned from the ESD test
          idx = gESDtest(ages_m,sig,max_outliers);
          outliers = ages_m(idx);

          ages(idx,:) = NaN; % Remove outliers from existing exposure age data

          new_ages_ka.Al26_outliers = outliers; % Export outliers
          new_ages_ka.Al26 = ages; % Export ages without outliers
                    
          if ~isempty(outliers)
              disp('Outliers (Al-26):');
              for a = 1:length(outliers)
                  name = cellstr(sample_names(a));
                  disp(['sample ' name{1} ' (mean of ' sprintf('%0.2f',outliers(a)) ' ka)']);
              end
          else
              disp('Generalised ESD test found no outliers in Al-26 data');
          end
      end
      
  end
  
end


%%%%%%%% Perform the generalised extreme Studentized deviate test %%%%%%%%

function out = gESDtest(ages_m,sig,max_outliers)
  % This function adopts the framework in MATLAB's isoutlier, implementing 
  % the gESD test to identify potential outliers.
  
  % Check inputs
  if size(ages_m,2) > 1
      ages_m = ages_m';
  end
  
  ages_size = size(ages_m);
  ncols = prod(ages_size(2:end));
  lowerbound = NaN([1 ages_size(2:end)]);
  upperbound = lowerbound;
  center = lowerbound;
  agesflat = ages_m(:,:);
  
  out = false(ages_size);
  for j=1:ncols
      indvec = (j-1)*size(agesflat,1)+1:j*size(agesflat,1); % linear indices
      agestemp = agesflat(:,j);
      indvec(isnan(agestemp)) = [];
      agestemp(isnan(agestemp)) = [];
      n = numel(agestemp);
      if n > 0
          ages_mean = zeros(max_outliers,1);
          ages_std = zeros(max_outliers,1);
          lambda = zeros(max_outliers,1);
          R = zeros(max_outliers,1);
          Rloc = zeros(max_outliers,1);
          
          for i = 1:max_outliers
              ages_mean(i) = mean(agestemp, 'omitnan');
              ages_std(i) = std(agestemp, 'omitnan');
              [ages_max,loc] = max(abs(agestemp - ages_mean(i)));
              R(i) = ages_max/ages_std(i);
              agestemp(loc) = [];
              Rloc(i) = indvec(loc);
              indvec(loc) = [];
              
              % compute lambda
              pp = 1 - sig / (2*(n-i+1));
              t = tinv(pp,n-i-1);
              lambda(i) = (n-i)*t/sqrt((n-i-1+t.^2)*(n-i+1));
          end
          
          lastindex = find(R > lambda, 1, 'last');
          
          out(Rloc(1:lastindex)) = true;
          if isempty(lastindex)
              tindex = max_outliers;
          else
              tindex = min(lastindex+1,max_outliers);
          end
          center(j) = ages_mean(tindex);
          lowerbound(j) = ages_mean(tindex) - ages_std(tindex)*lambda(tindex);
          upperbound(j) = ages_mean(tindex) + ages_std(tindex)*lambda(tindex);
      end
  end
  
  out = (ages_m > upperbound | ages_m < lowerbound);

end
