%
% out = plot_kernel(ages_ka,feature,save_plot)
% out = plot_kernel(ages_ka,feature,save_plot,mask)
% out = plot_kernel(ages_ka,feature,save_plot,mask,x_lim,weight_type)
%
% Plots kernel density estimates for exposure ages, separately for each 
% nuclide. The mode of the distribution, weighted mean and standard 
% deviation, and reduced chi-squared with critical value are computed if 
% the ages are from the same feature.
%
% ages_ka is a required struct, containing ages and uncertainties 
% calculated using age_calc.m or imported using get_ages.m, and a logical 
% of the nuclides that were measured.
%
% feature should be a binary input of whether these ages come from a single
% feature [1] or not [0].
%
% save_plot is a binary input to save the figure(s) in .png and .eps 
% formats [1], or not [0].
%
% mask option can also be included to specify which samples to plot. This
% should be based on the original input sample data. Default is to plot all
% samples.
%
% x_lim is an optional input to specify the x-axis limits (in unit ka), in 
% the form [lower, upper].
%
% weighted is an optional binary input of whether the weighted mean and 
% standard deviation should be calculated [1] or not [0]. If 0, then the 
% arithmetric mean and standard deviation (unweighted) is used. The default
% is weighted. Applies only if feature = 1.
%
% Outputs the figure (fig) and axes (ax) handles.
%
% Written by Richard Selwyn Jones, Durham University
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite.
%
%
%%

function out = plot_kernel(ages_ka,feature,save_plot,mask,x_lim,weighted)

  % Check inputs
  if (nargin < 2 || nargin > 6)
      error('plot_kernel has wrong number of inputs!');
  end
  if (nargin < 3 || isempty(save_plot))
      save_plot = 0;
  end  
  if (nargin < 4 || isempty(mask))
      mask = 1:length(ages_ka.sample_names);
  end
  if (nargin < 6 || isempty(weighted))
      weighted = 1;
  end
  
  if strcmpi(ages_ka.scaling_model,'SF')
      scaling = 'LSD';
  elseif strcmpi(ages_ka.scaling_model,'SA')
      scaling = 'LSDn';
  else
      scaling = ages_ka.scaling_model;
  end
  
  % Suppress warnings
  warning('off','stats:chi2gof:LowCounts');

  % Get indices of nuclides measured
  NN_ind = find(ages_ka.NN);
  
  % Define colours
  d_colours = [0.86,0.098,0.106; 0.208,0.475,0.694; 0.291,0.66,0.279; 0.574,0.294,0.615; 0.962,0.479,0; 0.996,0.996,0.199]; % Dark
  l_colours = [0.98,0.703,0.68; 0.699,0.801,0.887; 0.797,0.918,0.77; 0.867,0.793,0.891; 0.992,0.848,0.648; 0.996,0.996,0.598]; % Light

  
  % Plot for each nuclide
  for a = 1:numel(NN_ind)
      
      N_ind = NN_ind(a);

      if N_ind == 1
          N_name = 'Be-10';
          col = [d_colours(1,:); l_colours(1,:)];
          ages_errs = ages_ka.Be10(mask,:);
          
      elseif N_ind == 2
          N_name = 'Al-26';
          col = [d_colours(2,:); l_colours(2,:)];
          ages_errs = ages_ka.Al26(mask,:);
      end

      % Remove any NaNs (outliers)
      ages_errs = ages_errs(all(~isnan(ages_errs),2),:);
      
      if feature == 1
          ages = ages_errs(:,1);
          errs = ages_errs(:,2); % Use interior uncertainties

          n_samples = numel(ages);

          if weighted == 1
              % Calculate weighted mean and standard deviation
              for b = 1:n_samples
                  weights(b,1) = (1/errs(b)) / sum(1./errs);
              end
              age_mean = sum(ages.*weights);
              age_stdev = sqrt(var(ages,weights));
              calc_txt = 'Wtd.';
          
          else
              % Calculate arithmetic mean and standard deviation
              age_mean = mean(ages);
              age_stdev = std(ages);
              calc_txt = 'Arith.';
              
          end
          
          % Perform reduced chi-squared test
          chi2_df = n_samples-1;
          for c = 1:n_samples
              xsig2(c) = ((ages(c) - age_mean)^2) / (errs(c)^2);
          end
          chi2_r = (1/chi2_df) * sum(xsig2);
          chi2_crit = 1+2 * sqrt(2/chi2_df);
          %p_val = chi2pval(chi2_r, chi2_df);
          
          if chi2_r <= chi2_crit
              chi_result = '(\chi^2_R < \kappa)';
          else
              chi_result = '(\chi^2_R > \kappa)';
          end
          

      else
          ages = ages_errs(:,1);
          errs = ages_errs(:,3); % Use exterior uncertainties
      end
      
      
      % Create figure
      fig = figure;
      hold on;
      h = fancy_pants_camelplot(ages,errs);
      set(h(1:end-1),'Color',col(2,:),'Linewidth',1.5,'Linestyle','-');
      set(h(end),'Color',col(1,:),'Linewidth',2,'Linestyle','-');
      ax=gca;
      YLim_max = ax.YLim(2);
      XLim_min = ax.XLim(1);
      XLim_max = ax.XLim(2);
      if feature == 1
          mode_val = max(h(end).YData); % Get mode of summed distribution
          mod_ind = h(end).YData == mode_val; % Get mode index
          mode_age = h(end).XData(mod_ind); % Find age of mode
          plot([mode_age mode_age],[0 YLim_max],'--k','linewidth',2);
          plot([age_mean age_mean],[0 YLim_max],'-k','linewidth',2);
          plot([age_mean-age_stdev age_mean-age_stdev],[0 YLim_max],':k','linewidth',1.5);
          plot([age_mean+age_stdev age_mean+age_stdev],[0 YLim_max],':k','linewidth',1.5);
          summary_text = {strcat('{Nuclide: }',N_name),strcat('{Scaling model: }',num2str(scaling)),strcat(calc_txt,'{ mean: }',num2str(round(age_mean,2)),'{ ka}'),strcat(calc_txt,'{ SD (1\sigma): }',num2str(round(age_stdev,2)),'{ ka}'),strcat('{Mode: }',num2str(round(mode_age,2)),'{ ka}'),strcat('{\chi^2_R: }',num2str(round(chi2_r,2)),'{  (}',num2str(round(chi2_df,1)),' d.f.)'),strcat('{\kappa: }',sprintf('%.2f',chi2_crit),'{  }',chi_result)};
      else
          summary_text = {strcat('{Nuclide: }',N_name),strcat('{Scaling model: }',num2str(scaling))};
      end
      hold off;
      ylim([0 YLim_max]);
      if (nargin == 5 || ~isempty(x_lim))
          xlim(x_lim);
      else
          if XLim_min < 0
              xlim([0 XLim_max]);
          else
              xlim([XLim_min XLim_max]);
          end
      end
      ax.XDir = 'reverse';
      box_dim = [0.15 0.70 0.19 0.2];
      annotation('textbox',box_dim,'String',summary_text,'FitBoxToText','on','Color',col(1,:),'FontSize',10,'FontWeight','bold');
      ylabel('Density');
      xlabel('Exposure age (ka)');
      box on;
      set(gca,'FontSize',10);
      
      
      % Save figure
      if save_plot == 1
          fig_name = strcat(ages_ka.ages_name,'_ExposureAge_KDE_',N_name,'_',num2str(scaling));
          export_fig(fig_name,'-png','-r300','-transparent'); % Save as a PNG (raster) file
          saveas(gcf,fig_name,'epsc'); % Save as a EPS (vector) file
      end

      % Export figure handles
      out{a} = [fig ax];
  
  end
  
end
