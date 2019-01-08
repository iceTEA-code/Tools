%
% out = transect_regress_fourier(ages_ka,transect_type)
% out = transect_regress_fourier(ages_ka,transect_type,n_terms)
% out = transect_regress_fourier(ages_ka,transect_type,n_terms,mask)
%
% Performs fourier-based regression analysis for exposure age data, in a 
% horizontal or vertical transect, to derive a profile of retreat/thinning 
% and corresponding rates through time. The analysis uses only mean values;
% neither the exposure age or sample position uncertainties are used. For a
% more robust analysis use transect_regress_spline.m, which accounts for 
% both age and position uncertainties within a Bayesian framework, and
% assumes that exposure ages represent either stability or continuous 
% retreat/thinning through time (i.e. no advance/thickening occurs).
%
% The MATLAB Curve Fitting Toolbox needs to be installed in order to run
% the analysis.
% 
% ages_ka is a required struct, containing ages and uncertainties 
% calculated using age_calc.m, elevations and positions, and a logical of 
% the nuclides that were measured.
%
% transect_type should be either "vert" or "horiz", which determines
% whether the transect is vertical (e.g. elevation above the modern ice) or
% horizontal (e.g. distance from modern ice terminus). The units are m/yr
% and km/yr respectively.
%
% n_terms is an optional input to specify the number of terms in the
% Fourier Series analysis, which should be a value between 1 and 8 (default
% is 3). The higher the number of terms, the more sinusoidal the fit.
%
% mask option can also be included to specify which samples to plot. This
% should be based on the original input sample data. Default is to plot all
% samples.
%
% Outputs various model parameters, the rates, the quantile estimates.
% plot_transect_spline_rates.m can subsequently be used to plot the
% results.
%
% Written by Richard Selwyn Jones, Durham University.
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite.
%
%
%%

function out = transect_regress_fourier(ages_ka,transect_type,n_terms,mask)
  
  % Check inputs
  if (nargin < 2 || nargin > 4)
      error('transect_regress_fourier has wrong number of inputs!');
  end
  
  if (~strcmp(transect_type,'vert') && ~strcmp(transect_type,'horiz'))
      error('transect_type should be "vert" or "horiz"!');
  end
  if (nargin < 3 || isempty(n_terms))
      n_terms = 3;
  end
  if (n_terms < 1 || n_terms > 8)
      error('n_terms should be a value between 1 and 8!');
  end
  if (nargin < 4 || isempty(mask))
      mask = 1:length(ages_ka.sample_names);
  end
  
  if ~license('test', 'Curve_Fitting_Toolbox')
      error('This analysis requires the Curve Fitting Toolbox to be installed.');
  end
  
  
  % Get relative position of samples (use elevation if input is NaN)
  pos = zeros(length(ages_ka.position),1);
  for a = 1:length(ages_ka.position)
      if isnan(ages_ka.position(a))
          pos(a) = ages_ka.elevation(a);
      else
          pos(a) = ages_ka.position(a);
      end
  end
  pos = pos(mask);
  
  
  % Get indices of nuclides measured
  NN_ind = find(ages_ka.NN);
  
 
  % Regress for each nuclide
  disp('Estimating rates from fourier series analysis...')
  
  for n = 1:numel(NN_ind)
      
      N_ind = NN_ind(n);
      
      % Get data for nuclide; Mask; Convert from ka to yrs
      if N_ind == 1
          age_m = ages_ka.Be10(mask,1) * 1000;
          age_err = ages_ka.Be10(mask,3) * 1000;
          
      elseif N_ind == 2
          age_m = ages_ka.Al26(mask,1) * 1000;
          age_err = ages_ka.Al26(mask,3) * 1000;
      end
      
      % Collate data
      pos_err = NaN(size(pos));
      %age_err = NaN(size(age_m));
      data = [pos,pos_err,age_m,age_err];
      
      % Perform regression
      n_sims = 10000;
 
      if n_terms == 1
          fourier_type = 'fourier1';
      elseif n_terms == 2
          fourier_type = 'fourier2';
      elseif n_terms == 3
          fourier_type = 'fourier3';
      elseif n_terms == 4
          fourier_type = 'fourier4';
      elseif n_terms == 5
          fourier_type = 'fourier5';
      elseif n_terms == 6
          fourier_type = 'fourier6';
      elseif n_terms == 7
          fourier_type = 'fourier7';
      elseif n_terms == 8
          fourier_type = 'fourier8';
      end
      
      y = pos;
      X = age_m;
      [fou_fit,gof,output] = fit(X,y,fourier_type);
      X_arr = round(min(X),-1):10:round(max(X),-1);
      fit_Y = fou_fit(X_arr);
      intopt = 'functional'; % 'observation','functional'
      simopt = 'off';
%       level = 0.68;
%       Y_bounds68 = predint(fou_fit,X_arr,level,intopt,simopt);
      level = 0.95;
      Y_bounds95 = predint(fou_fit,X_arr,level,intopt,simopt);
      
      % Smooth bounds
%       [fou_bounds,~,~] = fit(X_arr',Y_bounds68(:,1),fourier_type);
%       Y_bounds68(:,1) = fou_bounds(X_arr);
%       [fou_bounds,~,~] = fit(X_arr',Y_bounds68(:,2),fourier_type);
%       Y_bounds68(:,2) = fou_bounds(X_arr);
      [fou_bounds,~,~] = fit(X_arr',Y_bounds95(:,1),fourier_type);
      Y_bounds95(:,1) = fou_bounds(X_arr);
      [fou_bounds,~,~] = fit(X_arr',Y_bounds95(:,2),fourier_type);
      Y_bounds95(:,2) = fou_bounds(X_arr);
      
%       % Interpolate for yearly array
%       yr_arr = min(X):1:max(X);
%       fit_Y_yr = interp1(X_arr,fit_Y,yr_arr);
%       Y_bounds_yr(:,1) = interp1(X_arr,Y_bounds(:,1),yr_arr);
%       Y_bounds_yr(:,2) = interp1(X_arr,Y_bounds(:,2),yr_arr);
      
      reg.med_pos = fit_Y'; %fit_Y_yr;
%       reg.u68_pos = Y_bounds68(:,2)';
%       reg.l68_pos = Y_bounds68(:,1)';
      reg.u95_pos = Y_bounds95(:,2)'; %Y_bounds_yr(:,2)';
      reg.l95_pos = Y_bounds95(:,1)'; %Y_bounds_yr(:,1)';    
      
      % Sort analyis output     
      xbck = fliplr(X_arr); %yr_arr
%       l68bck = fliplr(reg.l68_pos);
      l95bck = fliplr(reg.l95_pos);
      plotting.reg_time = [X_arr,xbck]; %yr_arr
%       plotting.reg_pos68 = [reg.u68_pos,l68bck];
      plotting.reg_pos95 = [reg.u95_pos,l95bck];
      
      % Calculate rate
      rate = diff(reg.med_pos)./diff(X_arr);
      med_rate = [NaN,rate];
      rate_scale = med_rate./nanstd(med_rate);
%       u68_diff_scale = ((reg.u68_pos-reg.med_pos)./(std(reg.u68_pos-reg.med_pos)));
%       l68_diff_scale = ((reg.med_pos-reg.l68_pos)./(std(reg.med_pos-reg.l68_pos)));
%       rate_u68_scale = rate_scale + u68_diff_scale;
%       rate_l68_scale = rate_scale - l68_diff_scale;
%       rate_u68 = rate_u68_scale*nanstd(med_rate);
%       rate_l68 = rate_l68_scale*nanstd(med_rate);
      u95_diff_scale = ((reg.u95_pos-reg.med_pos)./(std(reg.u95_pos-reg.med_pos)));
      l95_diff_scale = ((reg.med_pos-reg.l95_pos)./(std(reg.med_pos-reg.l95_pos)));
      rate_u95_scale = rate_scale + u95_diff_scale;
      rate_l95_scale = rate_scale - l95_diff_scale;
      rate_u95 = rate_u95_scale*nanstd(med_rate);
      rate_l95 = rate_l95_scale*nanstd(med_rate);
      
      reg.med_rate = med_rate;
%       reg.u68_rate = rate_u68;
%       reg.l68_rate = rate_l68;
      reg.u95_rate = rate_u95;
      reg.l95_rate = rate_l95;
%       plotting.reg_rate68 = [reg.u68_rate,fliplr(reg.l68_rate)];
      plotting.reg_rate95 = [reg.u95_rate,fliplr(reg.l95_rate)];
      reg.yearsBP = X_arr; %yr_arr
      
      % Calculate max and min rates
      [max_rate,max_rate_idx] = max(reg.med_rate);
      max_rate_time = reg.yearsBP(max_rate_idx);
      [min_rate,min_rate_idx] = min(reg.med_rate);
      min_rate_time = reg.yearsBP(min_rate_idx);
      if min_rate < 0
          min_rate = 0;
      end
      
      if (strcmp(transect_type,'vert'))
          rate_factor = 100;
          disp(['Maximum median rate: ' sprintf('%0.2f',(max_rate*rate_factor)) ' cm/yr at ' sprintf('%0.2f',(max_rate_time/1000)) ' ka']) 
          disp(['Minimum median rate: ' sprintf('%0.2f',(min_rate*rate_factor)) ' cm/yr at ' sprintf('%0.2f',(min_rate_time/1000)) ' ka'])
          
      else
          rate_factor = 1000;
          disp(['Maximum median rate: ' sprintf('%0.2f',(max_rate*rate_factor)) ' m/yr at ' sprintf('%0.2f',(max_rate_time/1000)) ' ka'])
          disp(['Minimum median rate: ' sprintf('%0.2f',(min_rate*rate_factor)) ' m/yr at ' sprintf('%0.2f',(min_rate_time/1000)) ' ka'])
      end
      
      
      % Export outputs
      if N_ind == 1
          out.N10.elev_ages = data;
          out.N10 = reg;
          out.N10.plot_data = plotting;
      elseif N_ind == 2
          out.N26.elev_ages = data;
          out.N26 = reg;
          out.N26.plot_data = plotting;
      end
      
  end
  
  out.regress_type = 'fourier';
  disp('done.')

end
