%
% out = transect_regress_linear(ages_ka,transect_type,regress_type,save_plot)
% out = transect_regress_linear(ages_ka,transect_type,regress_type,save_plot,n_iter)
% out = transect_regress_linear(ages_ka,transect_type,regress_type,save_plot,n_iter,mask)
%
% Performs least-squares linear regression randomly to normally-distributed
% exposure ages (2 sigma) through a Monte Carlo simulation, for each 
% nuclide separately. Only regression slopes that are consistent with
% thinning/retreat are kept. Plots the probability distribution of modelled
% rates, and generates estimates at 68% and 95% confidence.
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
% regress_type corresponds to the type of least-squares regressions, and
% should be either "unweighted" or "weighted".
%
% save_plot is a binary input to save the figure(s) in .png and .eps 
% formats [1], or not [0].
%
% n_iter is an optional input to specify the number of Monte Carlo
% iterations.
%
% mask option can also be included to specify which samples to plot. This
% should be based on the original input sample data. Default is to plot all
% samples.
%
% Outputs various model parameters, the rates, the quantile estimates, and
% the figure (fig) and axes (ax) handles of the probability plot.
%
% Written by Richard Selwyn Jones, Durham University
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite.
%
%
%%

function out = transect_regress_linear(ages_ka,transect_type,regress_type,save_plot,n_iter,mask)

  % Check inputs
  if (nargin < 3 || nargin > 6)
      error('transect_regress_linear has wrong number of inputs!');
  end
  
  if (~strcmp(transect_type,'vert') && ~strcmp(transect_type,'horiz'))
      error('transect_type should be "vert" or "horiz"!');
  end
  if (~strcmp(regress_type,'unweighted') && ~strcmp(regress_type,'weighted'))
      error('regress_type should be "unweighted" or "weighted"!');
  end
  if (nargin < 4 || isempty(save_plot))
      save_plot = 0;
  end  
  if (nargin < 5 || isempty(n_iter))
      n_iter = 5000;
  end
  if (nargin < 6 || isempty(mask))
      mask = 1:length(ages_ka.sample_names);
  end
  
  % Suppress warnings
  warning('off','stats:regress:NoConst');
  
  
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
  for n = 1:numel(NN_ind)
      
      N_ind = NN_ind(n);

      % Get data for nuclide; Mask; Convert from ka to yrs
      if N_ind == 1
          N = 10;
          age_m = ages_ka.Be10(mask,1) * 1000;
          age_err = ages_ka.Be10(mask,3) * 1000;
          
      elseif N_ind == 2
          N = 26;
          age_m = ages_ka.Al26(mask,1) * 1000;
          age_err = ages_ka.Al26(mask,3) * 1000;
      end
  
  
      n_pos = length(pos);
      midx = mean(pos);
      X0=[ones(n_pos,1) pos-midx*ones(n_pos,1)]; % Axis changed to make midpoint the origin, for numerical stability
     
      
      W = ones(n_pos,1) ./ age_err;

      switch regress_type
          
          case 'unweighted'
              X = X0;
              W = ones(n_pos,1);
             
          case 'weighted'
              W = ones(n_pos,1) ./ age_err;
              ncol = 2;
              for b = 1:ncol
                  WX(:,b) = W .* X0(:,b);
              end
              X = WX;
      end
      
      y = W .* age_m;
      
      [B,B_ints] = LS_reg(y,X);
      modl = X0 * B;
      
      %DOF = n_pos-2;
      %t_crit = tinv(0.975,DOF);
      xm = [1.1*min(X0(:,2)):1:1.1*max(X0(:,2))]';
      n_xm = length(xm);
      
      %xmodl = B(1) * ones(n_xm,1)+B(2) * xm;
      %stdx = t_crit * sig_err * sqrt((1+1/n_pos) * ones(n_xm,1) + (1/((n_pos-1)*x_sd^2)) * (xm) .^2); % Formula for conventional 95% model CI
      Xm = [ones(n_xm,1), xm];

      B_rand_s = randn(n_iter,2);
      sig_B = (0.5/1.96) * (B_ints(:,2)-B_ints(:,1));
      B_rand_r = (B*ones(1,n_iter))' + [B_rand_s(:,1)*sig_B(1,1) B_rand_s(:,2)*sig_B(2,1)];
      [~, B_idx] = sort(B_rand_r(:,2));
      B_rand_s = B_rand_r(B_idx,1:2);
      bad_slope = B_rand_s(:,2) < 0;
      
      n_bad = sum(bad_slope); % Number of negative slopes
      B_rand = B_rand_s(n_bad+1:n_iter,1:2);
      
      rate = sort(ones(n_iter-n_bad,1)./B_rand(:,2));
      
      
      % Get rates at 68% and 95% quantiles
      quant_68 = [quantile(rate,.159,1); quantile(rate,.842,1)];
      quant_95 = [quantile(rate,.025,1); quantile(rate,.975,1)];
      
      % Get random models at 95% quantile
      rand_mod = Xm * B_rand';
      qu = quantile(rand_mod,.975,2);
      ql = quantile(rand_mod,.025,2);
      
      xm = xm + midx*ones(n_xm,1);
      d_xm = X0(:,2) + midx*ones(n_pos,1); % Readjust x-axis for plotting
      
      
      % Assimilate data
      reg.regress_type = regress_type;
      reg.quant_68 = quant_68;
      reg.quant_95 = quant_95;
      reg.rate = rate;
      reg.n_iter = n_iter;
      reg.n_bad = n_bad;
      

      % Print result
      if (strcmp(transect_type,'vert'))
          disp(['Thinning rate from ' regress_type ' least squares regression:'])
          disp([sprintf('%0.2f',quant_68(1)) '-' sprintf('%0.2f',quant_68(2)) ' m/yr  (68%)']);
          disp([sprintf('%0.2f',quant_95(1)) '-' sprintf('%0.2f',quant_95(2)) ' m/yr  (95%)']);
      else
          disp(['Retreat rate from ' regress_type ' least squares regression:'])
          disp([sprintf('%0.2f',quant_68(1)) '-' sprintf('%0.2f',quant_68(2)) ' km/yr  (68%)']);
          disp([sprintf('%0.2f',quant_95(1)) '-' sprintf('%0.2f',quant_95(2)) ' km/yr  (95%)']);
      end
      
      
      % Plot histogram of rates
      hist_plot = plot_rates_hist(ages_ka,reg,transect_type,N,save_plot);
      

      % Export
      plot.modl = modl;
      plot.d_xm = d_xm;
      plot.rand_mod = rand_mod;
      plot.xm = xm;
      plot.qu = qu;
      plot.ql = ql;
            
      if N_ind == 1
          out.N10 = reg;
          out.N10.plot_data = plot;
          out.N10.hist_plot_h = hist_plot;
      elseif N_ind == 2
          out.N26 = reg;
          out.N26.plot_data = plot;
          out.N26.hist_plot_h = hist_plot;
      end
  
  end
  
  
end


%%%%%%%%%%%% Perform multiple least squares linear regression %%%%%%%%%%%%%

function [coeff,coeff_ints] = LS_reg(y,X)
  
  alpha = 0.05;  % Use 95% confidence intervals (100-alpha)
  
  [n,ncolX] = size(X);
  
  % Remove any missing values
  wasnan = (isnan(y) | any(isnan(X),2));
  havenans = any(wasnan);
  if havenans
     y(wasnan) = [];
     X(wasnan,:) = [];
     n = length(y);
  end

  % Remove dependent columns of X using orthogonal-triangular decomposition
  [Q,R,perm] = qr(X,0);
  if isempty(R)
      p = 0;
  elseif isvector(R)
      p = double(abs(R(1))>0);
  else
      p = sum(abs(diag(R)) > max(n,ncolX)*eps(R(1)));
  end
  if p < ncolX
      R = R(1:p,1:p);
      Q = Q(:,1:p);
      perm = perm(1:p);
  end
  
  
  % Compute the LS coefficients
  coeff = zeros(ncolX,1);
  coeff(perm) = R \ (Q'*y);
  
  % Find confidence intervals for x
  RI = R\eye(p);
  nu = max(0,n-p);
  yhat = X*coeff;
  r = y-yhat;
  normr = norm(r);
  if nu ~= 0
      rmse = normr/sqrt(nu);
      tval = tinv((1-alpha/2),nu);
  else
      rmse = NaN;
      tval = 0;
  end
  se = zeros(ncolX,1);
  se(perm,:) = rmse*sqrt(sum(abs(RI).^2,2));
  coeff_ints = [coeff-tval*se, coeff+tval*se];

end
