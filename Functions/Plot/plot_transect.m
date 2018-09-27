%
% out = plot_transect(ages_ka,transect_type,save_plot)
% out = plot_transect(ages_ka,transect_type,save_plot,mask)
% out = plot_transect(ages_ka,transect_type,save_plot,mask,time_lim,pos_lim)
%
% Plots exposure ages as either a horizontal or vertical transect, 
% separately for each nuclide.
%
% ages_ka is a required struct, containing ages and uncertainties 
% calculated using age_calc.m, elevations and positions, and a logical of 
% the nuclides that were measured.
%
% transect_type should be either "vert" or "horiz", which determines
% whether the transect is vertical (e.g. elevation above the modern ice) or
% horizontal (e.g. distance from modern ice terminus). The units are m and
% km, respectively. The inputted relative position is used for each sample,
% unless it is NaN, in which case the elevation (m asl) is used.
%
% save_plot is a binary input to save the figure(s) in .png and .eps 
% formats [1], or not [0].
%
% mask option can also be included to specify which samples to plot. This
% should be based on the original input sample data. Default is to plot all
% samples.
%
% Optionally set the limit for the age axis with time_lim ([lower upper]), 
% and limit for the relative position axis with pos_lim ([lower upper]).
%
% Outputs the figure (fig) and axes (ax) handles.
%
% Written by Richard Selwyn Jones, Durham University
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite.
%
%
%%

function out = plot_transect(ages_ka,transect_type,save_plot,mask,time_lim,pos_lim)

  % Check inputs
  if (nargin < 2 || nargin > 6)
      error('plot_transect has wrong number of inputs!');
  end
  if (~strcmp(transect_type,'vert') && ~strcmp(transect_type,'horiz'))
      error('transect_type should be "vert" or "horiz"!');
  end
  if (nargin < 3 || isempty(save_plot))
      save_plot = 0;
  end  
  if (nargin < 4 || isempty(mask))
      mask = 1:length(ages_ka.sample_names);
  end
  
  if strcmpi(ages_ka.scaling_model,'SF')
      scaling = 'LSD';
  elseif strcmpi(ages_ka.scaling_model,'SA')
      scaling = 'LSDn';
  else
      scaling = ages_ka.scaling_model;
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
  

  % Get indices of nuclides measured
  NN_ind = find(ages_ka.NN);
  
  % Define colours
  d_colours = [0.86,0.098,0.106; 0.208,0.475,0.694; 0.291,0.66,0.279; 0.574,0.294,0.615; 0.962,0.479,0];
  l_colours = [0.98,0.703,0.68; 0.699,0.801,0.887; 0.797,0.918,0.77; 0.867,0.793,0.891; 0.992,0.848,0.648];

  
  % Plot for each nuclide
  for b = 1:numel(NN_ind)
      
      N_ind = NN_ind(b);

      if N_ind == 1
          N_name = 'Be-10';
          col = [d_colours(1,:); l_colours(1,:)];
          ages_errs = ages_ka.Be10(mask,:);
          N_pos = pos(ages_ka.logical_10);
          N_pos = N_pos(mask);
      elseif N_ind == 2
          N_name = 'Al-26';
          col = [d_colours(2,:); l_colours(2,:)];
          ages_errs = ages_ka.Al26(mask,:);
          N_pos = pos(ages_ka.logical_26);
          N_pos = N_pos(mask);
      end
      
      
      % Create figure
      fig = figure;      
      hold on;
      
      if (strcmp(transect_type,'vert'))
          h = errbar(ages_errs(:,1),N_pos,ages_errs(:,2),'-','horiz');
          set(h,'Color',col(1,:),'LineWidth',1.5);
          plot(ages_errs(:,1),N_pos,'o','MarkerSize',6,'MarkerEdgeColor',col(1,:),'MarkerFaceColor',col(2,:),'LineWidth',1);
          x_lab = 'Exposure age (ka)';
          y_lab = 'Relative elevation (m)';
          ax=gca;
          if (nargin > 4 && ~isempty(time_lim))
              xlim(time_lim);
          else
              xlim(ax.XLim);
          end
          if (nargin > 4 && ~isempty(pos_lim))
              ylim(pos_lim);
          else
              ylim(ax.YLim);
          end
          ax.XDir = 'reverse';
          box_dim = [0.61 0.7 0.2 0.2];
                
      elseif (strcmp(transect_type,'horiz'))
          if nnz(N_pos) == 0
             warning('All relative positions for the horizontal transect are zero.') 
          end          
          h = errbar(N_pos,ages_errs(:,1),ages_errs(:,2),'-');
          set(h,'Color',col(1,:),'LineWidth',2);
          plot(N_pos,ages_errs(:,1),'o','MarkerSize',6,'MarkerEdgeColor',col(1,:),'MarkerFaceColor',col(2,:),'LineWidth',1);
          x_lab = 'Relative distance (km)';
          y_lab = 'Exposure age (ka)';
          ax=gca;
          if (nargin > 3 && ~isempty(time_lim))
              ylim(time_lim);
          else
              ylim(ax.YLim);
          end
          if (nargin > 3 && ~isempty(pos_lim))
              xlim(pos_lim);
          else
              xlim(ax.XLim);
          end
          box_dim = [0.15 0.69 0.1 0.2];
      end
      hold off;
      ylim(ax.YLim);
      summary_text = {strcat('{Nuclide: }',N_name),strcat('{Scaling model: }',scaling)};
      annotation('textbox',box_dim,'String',summary_text,'FitBoxToText','on','Color',col(1,:),'FontSize',10,'FontWeight','bold');
      xlabel(x_lab);
      ylabel(y_lab);
      grid on;
      box on;
      set(gca,'FontSize',10);
      
      
      % Save figure
      if save_plot == 1
          fig_name = strcat(ages_ka.ages_name,'_ExposureAge_transect_',N_name,'_',scaling);
          export_fig(fig_name,'-png','-r300','-transparent'); % Save as a PNG (raster) file
          saveas(gcf,fig_name,'epsc'); % Save as a EPS (vector) file
      end

      % Export figure handles
      out{b} = [fig ax];
  
  end
  
end
