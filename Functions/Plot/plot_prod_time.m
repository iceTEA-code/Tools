%
% out = plot_prod_time(prod_input_1)
% out = plot_prod_time(prod_input_1,save_plot)
% out = plot_prod_time(prod_input_1,save_plot,ages_name)
% out = plot_prod_time(prod_input_1,save_plot,ages_name,y_lim)
% out = plot_prod_time(prod_input_1,save_plot,ages_name,y_lim,prod_input_2)
%
% Plots sample-specific production rates through time, separately
% for each nuclide. These rates were used to calculate the exposure ages of
% the samples.
% 
% prod_input_1 is a required struct, containing production rates 
% calculated using age_calc.m or elev_correct.m, a logical of the nuclides that were 
% measured, and the scaling model used. If prod_input_2 is also inluded,
% then prod_input_1 is the uncorrected output of elev_correct.m.
%
% save_plot is a binary input to save the figure(s) in .png and .eps 
% formats [1], or not [0].
%
% ages_name is string input for the name of the data, which will be used
% to save the figure. It is required if save_plot is 1.
%
% y_lim is an optional input to specify the y-axis limits (in unit 
% atoms/g/year), in the form [lower, upper].
%
% prod_input_2 is an optional struct, containing the same production rate 
% information as in prod_input_1, but elevation-corrected using 
% elev_correct.m.
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

function out = plot_prod_time(prod_input_1,save_plot,ages_name,y_lim,prod_input_2)

  % Check inputs
  if (nargin < 1 || nargin > 5)
      error('plot_prod_time has wrong number of inputs!');
  end
  if (nargin < 2 || isempty(save_plot))
      save_plot = 0;
  end
  if save_plot == 1 && (nargin < 3 || isempty(ages_name))
      error('plot_prod_time requires ages_name to save the figure');
  end
  matver = version('-date'); matver = str2num(matver(end-3:end)); % Check MATLAB version
  

  % Get indices of nuclides measured
  NN_ind = find(prod_input_1.NN);
  
  % Define colours
  colours = [0.86,0.098,0.106; 0.208,0.475,0.694; 0.291,0.66,0.279; 0.574,0.294,0.615; 0.962,0.479,0];

  % Plot for each nuclide
  for a = 1:numel(NN_ind)
      
      N_ind = NN_ind(a);

      if N_ind == 1
          N_name = 'Be-10';
          col = colours(1,:);
          prod_time_1 = prod_input_1.prod_time_10;
          if nargin > 4
              prod_time_2 = prod_input_2.prod_time_10;
          end
      elseif N_ind == 2
          N_name = 'Al-26';
          col = colours(2,:);
          prod_time_1 = prod_input_1.prod_time_26;
          if nargin > 4
              prod_time_2 = prod_input_2.prod_time_26;
          end
      end

      
      % Create figure
      fig = figure;
      l_width = 1.5;
      hold on;
      for b = 1:length(prod_time_1)
          if nargin > 4
              plot(prod_time_1{b}(2:end,1),prod_time_1{b}(2:end,2),'--','Color',col,'LineWidth',l_width);
              plot(prod_time_2{b}(2:end,1),prod_time_2{b}(2:end,2),'-','Color',col,'LineWidth',l_width);
              leg = legend('uncorrected','corrected','Location','Best'); %'NorthEast'
          else
              plot(prod_time_1{b}(2:end,1),prod_time_1{b}(2:end,2),'-','Color',col,'LineWidth',l_width);
          end
      end
      hold off;
      ax = gca;
      summary_text = {strcat('{Nuclide: }',N_name),strcat('{Scaling model: }',num2str(prod_input_1.scaling_model))};
      if nargin > 4
          if matver > 2015
              title(leg,summary_text,'Color',col,'FontSize',10,'FontWeight','bold');
          end
      else
          annotation('textbox',[0.61 0.69 0.2 0.2],'String',summary_text,'FitBoxToText','on','Color',col,'FontSize',10,'FontWeight','bold');
      end
      hold off;
      ylabel({'Sample-specific production rate','(atoms g^{-1} yr^{-1})'});
      xlabel('Time exposed (ka)');
      if (nargin > 3 && ~isempty(y_lim))
          xlim(y_lim);
      else
          ylim(ax.YLim);
      end
      xlim([0 ax.XLim(2)]);
      ax.XDir = 'reverse';
      grid on;
      box on;
      set(gca,'FontSize',10);
      
      
      % Save figure
      if save_plot == 1
          fig_name = strcat(ages_name,'_ProdVsTime_',N_name,'_',num2str(prod_input_1.scaling_model));
          export_fig(fig_name,'-png','-r400'); % Save as a PNG (raster) file
          saveas(gcf,fig_name,'epsc'); % Save as a EPS (vector) file
      end

      % Export figure handles
      out{a} = [fig ax];
      
  end
      
end