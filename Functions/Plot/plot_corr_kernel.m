%
% out = plot_corr_kernel(uncorr_ages_ka,corr_ages_ka,save_plot)
% out = plot_corr_kernel(uncorr_ages_ka,corr_ages_ka,save_plot,mask)
% out = plot_corr_kernel(uncorr_ages_ka,corr_ages_ka,save_plot,mask,x_lim)
%
% Plots kernel density estimates for both uncorrected and corrected
% exposure ages, separately for each nuclide.
%
% uncorr_ages_ka is a required struct, containing uncorrected ages and 
% uncertainties calculated using age_calc.m, and a logical of the nuclides
% that were measured.
%
% corr_ages_ka is a required struct, containing corrected ages and 
% uncertainties calculated using elev_correct.m within the plot substruct
% (e.g. corrected.plot.corr).
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
% Outputs the figure (fig) and axes (ax) handles.
%
% Written by Richard Selwyn Jones, Durham University
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite, which is built on versions of 
% CRONUS-Earth and CRONUScalc code.
%
%
%%

function out = plot_corr_kernel(uncorr_ages_ka,corr_ages_ka,save_plot,mask,x_lim)

  % Check inputs
  if (nargin < 2 || nargin > 5)
      error('plot_corr_kernel has wrong number of inputs!');
  end
  if (nargin < 3 || isempty(save_plot))
      save_plot = 0;
  end  
  if (nargin < 4 || isempty(mask))
      mask = 1:length(uncorr_ages_ka.sample_names);
  end
  
  if strcmpi(corr_ages_ka.scaling_model,'SF')
      corr_scaling = 'LSD';
  elseif strcmpi(corr_ages_ka.scaling_model,'SA')
      corr_scaling = 'LSDn';
  else
      corr_scaling = corr_ages_ka.scaling_model;
  end
  if strcmpi(uncorr_ages_ka.scaling_model,'SF')
      uncorr_scaling = 'LSD';
  elseif strcmpi(uncorr_ages_ka.scaling_model,'SA')
      uncorr_scaling = 'LSDn';
  else
      uncorr_scaling = uncorr_ages_ka.scaling_model;
  end
  if ~strcmp(corr_scaling,uncorr_scaling)
      warning(['Different scaling models were used to calculate ages (uncorrected, ' uncorr_scaling '; corrected, ' corr_scaling ')!']);
  end

  % Get indices of nuclides measured
  NN_ind = find(uncorr_ages_ka.NN);
  
  % Define colours
  d_colours = [0.86,0.098,0.106; 0.208,0.475,0.694; 0.291,0.66,0.279; 0.574,0.294,0.615; 0.962,0.479,0; 0.996,0.996,0.199]; % Dark
  l_colours = [0.98,0.703,0.68; 0.699,0.801,0.887; 0.797,0.918,0.77; 0.867,0.793,0.891; 0.992,0.848,0.648; 0.996,0.996,0.598]; % Light

  
  % Plot for each nuclide
  for a = 1:numel(NN_ind)
      
      N_ind = NN_ind(a);
      
      if N_ind == 1
          N_name = 'Be-10';
          col = [d_colours(1,:); l_colours(1,:)];
          uncorr_ages_errs = uncorr_ages_ka.Be10(mask,:);
          corr_ages_errs = corr_ages_ka.Be10(mask,:);
          
      elseif N_ind == 2
          N_name = 'Al-26';
          col = [d_colours(2,:); l_colours(2,:)];
          uncorr_ages_errs = uncorr_ages_ka.Al26(mask,:);
          corr_ages_errs = corr_ages_ka.Al26(mask,:);
      end

      % Remove any NaNs (outliers)
      uncorr_ages_errs = uncorr_ages_errs(all(~isnan(uncorr_ages_errs),2),:);
      corr_ages_errs = corr_ages_errs(all(~isnan(corr_ages_errs),2),:);
      
      
      uncorr_ages = uncorr_ages_errs(:,1);
      uncorr_err = uncorr_ages_errs(:,2); % Use exterior uncertainties
      corr_ages = corr_ages_errs(:,1);
      corr_err = corr_ages_errs(:,2); % Use exterior uncertainties
      
      
      % Create figure
      fig = figure;
      hold on;
      h1 = fancy_pants_camelplot(uncorr_ages,uncorr_err);
      set(h1(1:end-1),'Color',col(2,:),'Linewidth',1.5,'Linestyle','-');
      delete(h1(end)); % Remove summed probability line
      h2 = fancy_pants_camelplot(corr_ages,corr_err);
      set(h2(1:end-1),'Color',col(1,:),'Linewidth',1.5,'Linestyle','-');
      hold off;
      delete(h2(end)); % Remove summed probability line
      
      leg = legend([h1(1),h2(1)],'uncorrected','corrected','Location','NorthWest'); %legend('boxoff');
      ax = gca;

      summary_text = {strcat('{Nuclide: }',N_name),strcat('{Scaling model: }',num2str(corr_scaling))};
      title(leg,summary_text,'Color',col(1,:),'FontSize',10,'FontWeight','bold');
      ylim(ax.YLim);
      if (nargin == 5 && ~isempty(x_lim))
          xlim(x_lim);
      else
          xlim(ax.XLim);
      end
      ax.XDir = 'reverse';
      ylabel('Density');
      xlabel('Exposure age (ka)');
      box on;
      set(gca,'FontSize',10);
      
      
      % Save figure
      if save_plot == 1
          fig_name = strcat(corr_ages_ka.ages_name,'_ExposureAge_KDE_',N_name,'_',num2str(corr_scaling));
          export_fig(fig_name,'-png','-r300'); % Save as a PNG (raster) file
          saveas(gcf,fig_name,'epsc'); % Save as a EPS (vector) file
      end

      % Export figure handles
      out{a} = [fig ax];
  
  end
  
end
