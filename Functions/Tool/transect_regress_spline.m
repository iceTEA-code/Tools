%
% out = transect_regress_spline(ages_ka,transect_type)
% out = transect_regress_spline(ages_ka,transect_type,n_iter)
% out = transect_regress_spline(ages_ka,transect_type,n_iter,mask)
%
% Performs penalised spline regression for exposure age data, in a 
% horizontal or vertical transect, to derive a profile of retreat/thinning 
% and corresponding rates through time. Both the normally-distributed 
% exposure ages (2 sigma) and sample elevation uncertainties are used 
% within a Bayesian framework. It is assumed that exposure ages represent 
% either stability or continuous retreat/thinning through time (i.e. no 
% advance/thickening occurs).
%
% Just Another Gibbs Sampler (JAGS) is used to efficiently perform the
% analysis. If not found, then the program is downloaded and installed.
% MATJAGs is used as the Matlab interface for JAGS.
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
% n_iter is an optional input to specify the number of Monte Carlo
% iterations. Default is 20000.
%
% mask option can also be included to specify which samples to plot. This
% should be based on the original input sample data. Default is to plot all
% samples.
%
% Outputs various model parameters, the rates, the quantile estimates.
% plot_transect_spline_rates.m can subsequently be used to plot the
% results.
%
% Written by Niamh Cahill, University College Dublin, and
% Richard Selwyn Jones, Durham University.
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite.
%
%
%%

function out = transect_regress_spline(ages_ka,transect_type,n_iter,mask)
  
  % Check that the system can find JAGS
  if ispc
      [status,~] = system('where "jags-terminal.exe"');
  elseif ismac
      [status,~] = system('which jags');
      if status == 1
          status = 0;
      end
  end
  % Install if it cannot be found
  if status ~= 0
      jags_id = get_jags();
      % Re-check that the system can find JAGS
      if jags_id == 1 && ispc              
          os = computer('arch');
          if strcmp(os,'win32')
              jagsbin = 'x32\bin;';
          elseif strcmp(os,'win64')
              jagsbin = 'x64\bin;';
          end
          jagspath = strcat(num2str(pwd),'\Functions\JAGS\',jagsbin);
          setenv('PATH', [jagspath getenv('PATH')]); % Assign to system path
          [status,~] = system('where "jags-terminal.exe"'); % Check
      elseif jags_id == 1 && ismac
          [status,~] = system('which jags');
          if status == 1
              status = 0;
          end
      end
      if jags_id == 0 && status ~= 0 % If JAGS wasn't installed and can't be found, don't continue
          warning('JAGS was not installed and can''t be found on the system, so will not continue. Install manually and try to run the tool again.')
          return
      end
  end

  
  % Check inputs
  if (nargin < 2 || nargin > 4)
      error('transect_regress_spline has wrong number of inputs!');
  end
  
  if (~strcmp(transect_type,'vert') && ~strcmp(transect_type,'horiz'))
      error('transect_type should be "vert" or "horiz"!');
  end
  if (nargin < 3 || isempty(n_iter))
      n_iter = 20000;
  end
  if (nargin < 4 || isempty(mask))
      mask = 1:length(ages_ka.sample_names);
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
  if strcmp(transect_type,'vert')
      pos_err = ages_ka.elevation_err';
  else
      pos_err = zeros(length(pos),1); % Have to assume that uncertainty on the position is zero
  end
  pos = pos(mask);
  pos_err = pos_err(mask);
  
  
  % Get indices of nuclides measured
  NN_ind = find(ages_ka.NN);
  
 
  % Regress for each nuclide
  disp('Estimating rates from Bayesian penalised spline analysis...')
  
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
      data = [pos,pos_err,age_m,age_err];
      
      % Sort and set up the data
      standardise = 1;
      data_sorted = sortrows(data,3);
      
      Position = data_sorted(:,1);
      PositionError = data_sorted(:,2);
      Age = data_sorted(:,3);
      AgeError = data_sorted(:,4);
      
      if logical(standardise)
          year = (-Age-mean(-Age))/std(-Age);
          N = length(year);
          position = (Position-mean(Position))/std(Position);
          sigma_y = PositionError/std(Position);
          sigma_x = AgeError/std(-Age);
      end
      if ~logical(standardise)
          year = -Age/1000;
          N = length(year);
          position = Position;
          sigma_y = PositionError;
          sigma_x = AgeError/1000;
      end
      
      
      % Create required JAGS data
      
      % Because we have x error we need to create the spline basis functions within JAGS
      % Use function getknots to get the knots for the splines
      nseg = 20;
      res = getknots(year,nseg,[]);
      
      % JAGS data needed for basis functions
      n_knots = length(res.knots);
      D = diff(diag(ones(n_knots,1)),res.deg + 1) / (gamma(res.deg + 1) * res.dx ^ res.deg);
      K = length(D(:,1));
      
      % Set up prediction data
      x_spacing = 100;
      x_star = linspace(min(year),max(year),x_spacing);
      n_grid = length(x_star);
      
      mydata = struct('N',N, ...
          'y',position, ...
          'x',year, ...
          'xstar',x_star, ...
          'sigmay',sigma_y, ...
          'sigmax',sigma_x, ...
          'nknots',n_knots, ...
          'deg',res.deg, ...
          'D',D, ...
          'knots',res.knots, ...
          'ngrid',n_grid, ...
          'K',K);
      
      % Parameters to look at/save
      mypars = {'betak','sigmadelta','muypred'};
      
      % MCMC settings
      n_chains = 4; % Number of MCMC chains
      n_burnin = 4000; % Number of burnin steps
      n_iterations = n_iter; % Number of samples (iterations) after burnin
      n_thin = 4; % Keep every n'th step
      n_samples = n_iterations/(n_thin+1); % Determine number of samples per chain following burnin and thinning out
      
      % Set initial values
      for i = 1:n_chains
          in_val.nuy = 0.5;
          in_val.sigmadelta = 1;
          %in_val.betak(1) = 1;
          init_vals(i) = in_val;
      end
      
      
      % Run the model
      if ispc
          spline_mod_path = '\Functions\Tool\SplineModel.txt';
          jagstmp_path = '\Functions\jagstmp';
      else
          spline_mod_path = '/Functions/Tool/SplineModel.txt';
          jagstmp_path = '/Functions/jagstmp';
      end
      [~, ~, mcmcmodel_jags] = matjags( ...
          mydata, ...
          fullfile(pwd,spline_mod_path), ...
          init_vals, ...
          'workingDir', strcat(num2str(pwd),jagstmp_path), ...
          'nchains', n_chains, ...
          'nburnin', n_burnin, ...
          'nsamples', n_samples, ...
          'thin', n_thin, ...
          'monitorparams', mypars );
      
      
      % Sort outputs
      muypred_struct = cat(3,mcmcmodel_jags(:).muypred);
      muypred_struct_sort = permute(muypred_struct,[1,3,2]);
      mu_ypred = reshape(muypred_struct_sort,n_samples*n_chains,x_spacing);
      if logical(standardise)
          ypred = ((mu_ypred)*std(Position))+mean(Position);
          x_star = -((x_star*std(-Age))+mean(-Age));
      end
      if ~logical(standardise)
          ypred = mu_ypred;
          x_star = -x_star*1000;
      end
      
      reg.med_y = median(ypred);
      reg.u68_pos = quantile(ypred,0.84);
      reg.l68_pos = quantile(ypred,0.16);
      l68bck = fliplr(reg.l68_pos);
      reg.u95_pos = quantile(ypred,0.975);
      reg.l95_pos = quantile(ypred,0.025);
      l95bck = fliplr(reg.l95_pos);
      xbck = fliplr(x_star);
      plotting.reg_time = [x_star,xbck];
      plotting.reg_pos68 = [reg.u68_pos,l68bck];
      plotting.reg_pos95 = [reg.u95_pos,l95bck];
      
      
      % Calculate rate
      ypred_diff = diff(ypred')';
      x_star_diff = diff(x_star);
      rate = bsxfun(@rdivide,ypred_diff,x_star_diff);
      reg.med_rate = [0 median(rate)];
      reg.u68_rate = [0 quantile(rate,0.84)];
      reg.l68_rate = [0 quantile(rate,0.16)];
      reg.u95_rate = [0 quantile(rate,0.975)];
      reg.l95_rate = [0 quantile(rate,0.025)];
      plotting.reg_rate68 = [reg.u68_rate,fliplr(reg.l68_rate)];
      plotting.reg_rate95 = [reg.u95_rate,fliplr(reg.l95_rate)];
      reg.yearsBP = x_star;
      
      
      % Calculate max and min rates
      if reg.med_rate(1) == 0
          med_rates_mod = reg.med_rate(2:end);
          med_time_mod = reg.yearsBP(2:end);
      end
      [max_rate,max_rate_idx] = max(med_rates_mod);
      max_rate_time = med_time_mod(max_rate_idx);
      [min_rate,min_rate_idx] = min(med_rates_mod);
      min_rate_time = med_time_mod(min_rate_idx);
      
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
          out.N10.elev_ages = data_sorted;
          out.N10 = reg;
          out.N10.plot_data = plotting;
      elseif N_ind == 2
          out.N26.elev_ages = data_sorted;
          out.N26 = reg;
          out.N26.plot_data = plotting;
      end
      
  end
  
  disp('done.')

end


%%%%%%%%%%%%% Spline function %%%%%%%%%%%%%

function out = getknots(x, nseg, deg)
  
  % Defaults
  if isempty(nseg)
      nseg = 30;
  end
  if isempty(deg)
      deg = 3;
  end
  
  xl = min(x);
  xr = max(x);
  
  % Construct B-spline basis
  dx = (xr - xl) / nseg;
  knots = xl - deg * dx : dx : xr + deg * dx;
  
  out.knots = knots;
  out.dx = dx;
  out.deg = 3;
  
end
