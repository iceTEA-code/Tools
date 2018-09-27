%
% out = elev_correct(sample_data,scaling_model,correction_type,elev_input)
%
% Computes a time-dependent production rate and exposure age for each
% sample and nuclide, corrected for elevation change through time.
% Calculations are based on a modified version of the CRONUScalc code
% (Marrero et al., 2016).
%
% sample_data is a struct containing necessary sample information produced 
% using get_data.m.
%
% scaling_model should be one of 'DE','DU','LI','ST','LM','LSD'/'SF' or 
% 'LSDn'/'SA'.
%
% correction_type corresponds to the method of elevation change correction,
% and should be either "model" or "rate" for GIA-derived or linear rate of 
% change, respectively.
%
% elev_input should be either "I5G" or "I6G" if correction_type is "model",
% or a value (m/ka) if correction_type is "rate".
%
% Outputs the uncorrected and corrected production rate, age and plotting 
% data for each nuclide, as well as the scaling scheme used.
%
% Written by Richard Selwyn Jones, Durham University
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite, which is built on versions of 
% CRONUS-Earth and CRONUScalc code.
%
%
%%

function out = elev_correct(sample_data,scaling_model,correction_type,elev_input)
  
  % Check inputs
  if (nargin ~= 4)
      error('elev_correct has wrong number of inputs!');
  end
  if (~strcmpi(correction_type,'rate') && ~strcmp(correction_type,'model'))
      error('correction_type should be "rate" or "model"!');
  end
  if (strcmpi(correction_type,'rate') && length(elev_input) ~= 1)
      error('the elev_input for correction_type "rate" should be a single value (m/ka)!');
  end
  if strcmpi(correction_type,'model') && (~strcmpi(elev_input,'I5G') && ~strcmpi(elev_input,'I6G'))
      error('the elev_input for correction_type "model" should be "I5G" or "I6G"!');
  end
  if strcmpi(scaling_model,'LSD')
      scaling = 'SF';
  elseif strcmpi(scaling_model,'LSDn')
      scaling = 'SA';
  else
      scaling = scaling_model;
  end
  
  
  % Create array for elevation change
  max_age_ka = 8160; % Max age (6* 10Be half-life = 8160 ka)
  t_interval = 100;  % Time interval (years)
  time_arr = 0:t_interval:max_age_ka*1000;
  
  if strcmpi(correction_type,'rate')
      elev_data.time_arr = time_arr;
      elev_change_rev = -elev_input;
      if elev_input == 0
          elev_data.elev_change = zeros(1,length(time_arr));
          warning('Rate of elevation change is set to zero.');
      else
          elev_data.elev_per_yr = elev_input/1000;
          elev_data.elev_change = 0:t_interval*elev_change_rev:(max_age_ka*1000)*elev_change_rev;
      end
      
  elseif strcmpi(correction_type,'model')
      % Get GIA-modelled elevations for samples
      disp('Getting GIA-modelled elevations for samples...')
      elev_data = get_GIA(sample_data,elev_input,time_arr);
  end
  
  elev_data.type = correction_type; % Append correction type
  
  
  % Get sample names
  for a = 1:length(sample_data.s)
      sample_names(a) = sample_data.s{a}.name;
  end
  
  
  % Determine what nuclides are present
  NN = [any(sample_data.logical_10) any(sample_data.logical_26)];
  
  
  % Setup outputs
  out.plot.uncorr.scaling_model = scaling_model;
  out.plot.uncorr.NN = NN;
  out.plot.corr = out.plot.uncorr;
  
  
  % Perform for each nuclide
  disp('Calculating elevation-corrected exposure ages...')
  
  if NN(1) % Be10
      nuclide = 10;
      logical_N10 = sample_data.logical_10;
      names_N10 = sample_names(logical_N10);
      data_N10 = sample_data.CC(logical_N10,:);
      data_N10_uncert = sample_data.CC_uncert(logical_N10,:);
      out.N10 = compute_correction(data_N10,data_N10_uncert,names_N10,nuclide,logical_N10,elev_data,scaling);
      out.plot.uncorr.prod_time_10 = out.N10.uncorr_prod_time;
      out.plot.corr.prod_time_10 = out.N10.corr_prod_time;
      out.plot.corr.Be10 = [out.N10.corr_mean_age/1000;out.N10.corr_age_int/1000;out.N10.corr_age_ext/1000]';
      out.plot.corr.sample_names = names_N10;
  end
  
  if NN(2) % Al26
      nuclide = 26;
      logical_N26 = sample_data.logical_26;
      names_N26 = sample_names(logical_N26);
      data_N26 = sample_data.CC(logical_N26,:);
      data_N26_uncert = sample_data.CC_uncert(logical_N26,:);
      out.N26 = compute_correction(data_N26,data_N26_uncert,names_N26,nuclide,logical_N26,elev_data,scaling);
      out.plot.uncorr.prod_time_26 = out.N26.uncorr_prod_time;
      out.plot.corr.prod_time_26 = out.N26.corr_prod_time;
      out.plot.corr.Al26 = [out.N26.corr_mean_age/1000;out.N26.corr_age_int/1000;out.N26.corr_age_ext/1000]';
      out.plot.corr.sample_names = names_N26;
  end
  
  out.scaling_model = scaling_model;
  out.correction_type = correction_type;
  out.elev_input = elev_input;
  
  disp('done.')
      
end



  %%%%%%%% Perform elevation corrections for the specified nuclide %%%%%%%%

  function corr_out = compute_correction(data_N,data_N_uncert,sample_names,nuclide,logical_N,elev_data,scaling_model)

  disp(' ');
  if nuclide == 10
      disp('Elevation-corrected ages (Be-10):');
  elseif nuclide == 26
      disp('Elevation-corrected ages (Al-26):');
  end
  
  % Do for each sample
  for a = 1:nnz(logical_N)
      
      % Get data for sample
      s_data = data_N(a,:);
      s_data_uncert = data_N_uncert(a,:);
      
      % Get physical constants
      pp = physpars(scaling_model); 
  
            
      % Compute for means
      
      % Create array of sample elevations
      if strcmpi(elev_data.type,'model')
          elev_arr = s_data(3) + elev_data.elev_change{a};
      else
          elev_arr(1,:) = s_data(3) + elev_data.elev_change;
          
          elev_arr_NaN = elev_arr < 0;
          elev_arr(elev_arr_NaN) = NaN; % Make elevations below zero NaNs
      end
      
      % Get parameters, scaling factors and ages
      [uncorr_mean_SF,corr_mean_SF,uncorr_mean_age,corr_mean_age,uncorr_times,uncorr_prod,corr_times,corr_prod] = get_fac_par_age(s_data,elev_data.time_arr,elev_arr,scaling_model,pp,nuclide,1);
      
      
      % Compute for uncertainties
      
      % Work through pressure and measurement uncertainties
      % Do only for elevation (pressure and measurement uncertainties)
      if nuclide == 10
          uncert_mask = [4,9];
      elseif nuclide == 26
          uncert_mask = [4,10];
      end
      derivs = zeros(length(uncert_mask),1);
      corr_age_uncert = 0.0;
      for b = 1:length(uncert_mask)
          mask_idx = uncert_mask(b);
          if (s_data_uncert(mask_idx) ~= 0.0)
              if (s_data(mask_idx) ~= 0.0)
                  this_delta = 0.01*abs(s_data(mask_idx));
              else
                  this_delta = 0.01;
              end
              delta_s_data = s_data;
              delta_s_data(mask_idx) = delta_s_data(mask_idx)+this_delta;
              
              % Create array of elevations if not present for uncertainties 
              if strcmpi(elev_data.type,'model')
                  delta_elev_arr = delta_s_data(3) + elev_data.elev_change{a};
              else
                  delta_elev_arr(1,:) = delta_s_data(3) + elev_data.elev_change;
              end
              
              % Get parameters, scaling factors and ages for uncertainties
              [uncorr_delta_SF,corr_delta_SF,uncorr_delta_age,corr_delta_age,~,~,~,~] = get_fac_par_age(delta_s_data,elev_data.time_arr,delta_elev_arr,scaling_model,pp,nuclide,0);
              d_aged = (corr_delta_age-corr_mean_age) / (this_delta);
              derivs(mask_idx) = d_aged;
              if (~isnan(s_data_uncert(mask_idx)))
                  corr_age_uncert = corr_age_uncert + (d_aged^2*s_data_uncert(mask_idx)^2);
              end
          end
      end
      
      % Add in terms for the uncertainty in production rates
      delta_pp = pp;
      if nuclide == 10
          delta_pp.PsBe = pp.PsBe+0.01*abs(pp.PsBe);
          [uncorr_deltapp_SF,corr_deltapp_SF,uncorr_deltapp_age,corr_deltapp_age,~,~,~,~] = get_fac_par_age(s_data,elev_data.time_arr,elev_arr,scaling_model,delta_pp,nuclide,0);
          d_aged_ps = (corr_deltapp_age-corr_mean_age) / (0.01*abs(pp.PsBe));
          corr_age_uncert_ext = corr_age_uncert + (d_aged_ps^2*pp.sigmaPsBe^2);
      elseif nuclide == 26
          delta_pp.PsAl = pp.PsAl+0.01*abs(pp.PsAl);
          [uncorr_deltapp_SF,corr_deltapp_SF,uncorr_deltapp_age,corr_deltapp_age,~,~,~,~] = get_fac_par_age(s_data,elev_data.time_arr,elev_arr,scaling_model,delta_pp,nuclide,0);
          d_aged_ps = (corr_deltapp_age-corr_mean_age) / (0.01*abs(pp.PsAl));
          corr_age_uncert_ext = corr_age_uncert + (d_aged_ps^2*pp.sigmaPsAl^2);
      end
      
      % Take the square root of uncertainty to get a standard deviation
      corr_age_uncert_int = sqrt(corr_age_uncert); % Internal uncertainty (analytical only)
      corr_age_uncert_ext = sqrt(corr_age_uncert_ext); % External/total uncertainty
      
      
      % Calculate the mean age difference
      age_diff = corr_mean_age - uncorr_mean_age;
      % Calculate the percentage of the age difference
      per_age_diff = (age_diff/uncorr_mean_age) * 100;
      
      
      % Calculate the average elevation since exposure
      [~,exp_ind] = min(abs(elev_data.time_arr-(corr_mean_age*1000))); % Find index of mean corrected exposure age (within 100 yrs)
      mean_elev = mean(elev_arr(1:exp_ind)); % Mean of elevations since then
      
      % Save elevation information
      uncorr_elev = ones(length(uncorr_prod),1) * s_data(3);
      uncorr_elev_time = [(0-uncorr_times)',uncorr_elev];
      corr_elev_time = [elev_data.time_arr(1:exp_ind)'/1000,elev_arr(1:exp_ind)'];
      
      % Save production information
      uncorr_prod_time = [(0-uncorr_times)',uncorr_prod'];
      corr_prod_time = [(0-corr_times)',corr_prod'];
            
      
      % Export
      corr_out.mean_elev(a) = mean_elev;
      corr_out.corr_mean_SF(a,:) = corr_mean_SF;
      corr_out.uncorr_mean_SF(a,:) = uncorr_mean_SF;
      corr_out.uncorr_mean_age(a) = uncorr_mean_age .* 1000;
      corr_out.corr_mean_age(a) = corr_mean_age .* 1000;
      corr_out.mean_age_diff(a) = age_diff .* 1000;
      corr_out.per_age_diff(a) = per_age_diff;
      corr_out.corr_age_int(a) = corr_age_uncert_int .* 1000;
      corr_out.corr_age_ext(a) = corr_age_uncert_ext .* 1000;
      corr_out.uncorr_elev_time{a} = uncorr_elev_time;
      corr_out.corr_elev_time{a} = corr_elev_time;
      corr_out.uncorr_prod_time{a} = uncorr_prod_time;
      corr_out.corr_prod_time{a} = corr_prod_time;
      
      
      % Print result
      name = cellstr(sample_names(a));
      disp(['Sample ' name{1} '  ' sprintf('%0.2f',corr_mean_age) ' +/- ' sprintf('%0.1f',corr_age_uncert_ext) '(' sprintf('%0.1f',corr_age_uncert_int) ') ka']);
      if age_diff < 0
          disp(['  ' sprintf('%0.1f',abs(age_diff)) ' ka (' sprintf('%0.0f',abs(per_age_diff)) '%) younger than the mean uncorrected age.']);
      elseif age_diff > 0
          disp(['  ' sprintf('%0.1f',abs(age_diff)) ' ka (' sprintf('%0.0f',abs(per_age_diff)) '%) older than the mean uncorrected age.']);
      end
      
  end

  end

  
  %%%%%%%%%%%%%% Compute scaling factors, parameters and age %%%%%%%%%%%%%%
  
  function [uncorr_SF,corr_SF,uncorr_age,corr_age,uncorr_times,uncorr_prod,corr_times,corr_prod] = get_fac_par_age(s_data,time_arr,elev_arr,scaling_model,pp,nuclide,mean_or_not)
  
  % Compute parameters
  
  if (nuclide == 10 || nuclide == 26)
      sp = samppars1026(s_data);
      sf = scalefacs1026(sp,scaling_model);
      max_age_ka = 8160;   % Max age (6* 10Be half-life = 8160 ka)
      max_depth = sp.depthtotop+max_age_ka*sp.epsilon+sp.ls*sp.rb+1000;
      cp = comppars1026(pp,sp,sf,max_depth);
  end
  
  tdsfsample.lat = sp.latitude;
  tdsfsample.long = sp.longitude;
  tdsfsample.elevation = sp.elevation;
  tdsfsample.pressure = sp.P;
  tdsfsample.scaling = scaling_model;
  load pmag_consts;
  
  
  % Compute nuclide production without correction
  
  tdsf_uncorr = get_tdsf(tdsfsample,pmag_consts);
  
  if strcmpi(scaling_model,'st')
      uncorr_SF = tdsf_uncorr.SF_St;
  elseif strcmpi(scaling_model,'du')
      uncorr_SF = tdsf_uncorr.SF_Du;
  elseif strcmpi(scaling_model,'de')
      uncorr_SF = tdsf_uncorr.SF_De;
  elseif strcmpi(scaling_model,'li')
      uncorr_SF = tdsf_uncorr.SF_Li;
  elseif strcmpi(scaling_model,'lm')
      uncorr_SF = tdsf_uncorr.SF_Lm;
  elseif strcmpi(scaling_model,'sf') || strcmpi(scaling_model,'LSD')
      uncorr_SF = tdsf_uncorr.SF_Sf;
  elseif strcmpi(scaling_model,'sa') || strcmpi(scaling_model,'LSDn')
      uncorr_SF = tdsf_uncorr.SF_Sa10;
  end
  
  
  % Compute nuclide production with correction
  
  tdsf_corr = get_tdsf_elev(tdsfsample,pmag_consts,time_arr,elev_arr);
  
  if strcmpi(scaling_model,'st')
      corr_SF = tdsf_corr.SF_St;
  elseif strcmpi(scaling_model,'du')
      corr_SF = tdsf_corr.SF_Du;
  elseif strcmpi(scaling_model,'de')
      corr_SF = tdsf_corr.SF_De;
  elseif strcmpi(scaling_model,'li')
      corr_SF = tdsf_corr.SF_Li;
  elseif strcmpi(scaling_model,'lm')
      corr_SF = tdsf_corr.SF_Lm;
  elseif strcmpi(scaling_model,'sf') || strcmpi(scaling_model,'LSD')
      corr_SF = tdsf_corr.SF_Sf;
  elseif strcmpi(scaling_model,'sa') || strcmpi(scaling_model,'LSDn')
      if nuclide == 10
          corr_SF = tdsf_corr.SF_Sa10;
      elseif nuclide == 26
          corr_SF = tdsf_corr.SF_Sa26;
      elseif nuclide == 3
          corr_SF = tdsf_corr.SF_Sa3;
      elseif nuclide == 14
          corr_SF = tdsf_corr.SF_Sa14;
      else
          corr_SF.Ca = tdsf_corr.SF_Sa36Ca;
	      corr_SF.K = tdsf_corr.SF_Sa36K;
	      corr_SF.Ti = tdsf_corr.SF_Sa36Ti;
	      corr_SF.Fe = tdsf_corr.SF_Sa36Fe;
      end
  end
  
  
  % Compute ages
  if nuclide == 10
      if mean_or_not == 1 % Do only for mean values
          uncorr_raw = computeage10(pp,sp,sf,cp,scaling_model); % Calculate for uncorrected
          uncorr_age = uncorr_raw(1);
          sf.currentsf = getcurrentsf(sf,0,scaling_model,'be'); % Get contemporary surface production rates (atoms/g)
          [~,~,ProdsBe,~,~,~] = prodz1026(0,pp,sf,cp);
          [uncorr_times,uncorr_prod] = getPlotData(uncorr_age,sf,ProdsBe,'be',scaling_model); % Get production data
      else
          uncorr_age = NaN;
          uncorr_times = NaN;
          uncorr_prod = NaN;
      end
      corr_sf = sf; corr_sf.tdsf = tdsf_corr; corr_sf.P = mean(tdsf_corr.pressure,'omitnan'); % Add new time-dependent scaling factor and pressure
      corr_cp = comppars1026(pp,sp,corr_sf,max_depth); % Compute parameters (including muon scaling)
      try
          corr_raw = computeage10(pp,sp,corr_sf,corr_cp,scaling_model); % Calculate corrected
      catch
          error('The sample elevation is below sea level at age of sample. Use a smaller rate of elevation change.');
      end
      corr_age = corr_raw(1);
      sf.currentsf = getcurrentsf(sf,0,scaling_model,'be'); % Get contemporary surface production rates (atoms/g)
      [~,~,ProdsBe,~,~,~] = prodz1026(0,pp,sf,cp);
      [corr_times,corr_prod] = getPlotData(corr_age,corr_sf,ProdsBe,'be',scaling_model); % Get production data

  elseif nuclide == 26
      if mean_or_not == 1 % Do only for mean values
          uncorr_raw = computeage26(pp,sp,sf,cp,scaling_model); % Calculate for uncorrected
          uncorr_age = uncorr_raw(1);
          sf.currentsf = getcurrentsf(sf,0,scaling_model,'al'); % Get contemporary surface production rates (atoms/g)
          [~,~,~,~,ProdsAl,~] = prodz1026(0,pp,sf,cp);
          [uncorr_times,uncorr_prod] = getPlotData(uncorr_age,sf,ProdsAl,'al',scaling_model); % Get production data          
      else
          uncorr_age = NaN;
          uncorr_times = NaN;
          uncorr_prod = NaN;
      end
      corr_sf = sf; corr_sf.tdsf = tdsf_corr; corr_sf.P = mean(tdsf_corr.pressure,'omitnan'); % Add new time-dependent scaling factor and pressure
      corr_cp = comppars1026(pp,sp,corr_sf,max_depth); % Compute parameters (including muon scaling)
      try
          corr_raw = computeage26(pp,sp,corr_sf,corr_cp,scaling_model); % Calculate corrected
      catch
          error('the sample elevation is below sea level at age of sample. Use a smaller rate of elevation change.');
      end
      corr_age = corr_raw(1);
      sf.currentsf = getcurrentsf(sf,0,scaling_model,'al'); % Get contemporary surface production rates (atoms/g)
      [~,~,~,~,ProdsAl,~] = prodz1026(0,pp,sf,cp);
      [corr_times,corr_prod] = getPlotData(corr_age,corr_sf,ProdsAl,'al',scaling_model); % Get production data
  end
  
  end
