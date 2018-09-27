%
% depth_plot = plot_concs_depth(sample_data)
% depth_plot = plot_concs_depth(sample_data,save_plot)
% depth_plot = plot_concs_depth(sample_data,save_plot,concs_name)
% depth_plot = plot_concs_depth(sample_data,save_plot,concs_name,x_lim,y_lim)
%
% Plots the nuclide concentrations vs. depth.
%
% sample_data is a required struct, created using get_data.m.
%
% save_plot is a binary input to save the figure(s) in .png and .eps 
% formats [1], or not [0].
%
% concs_name is string input for the name of the data, which will be used
% to save the figure. It is required if save_plot is 1.
%
% x_lim and y_lim set the limits for the x-axis and y-axis ([lower upper]).
% If y_lim is empty, then the uppermost depth is set to zero (the surface).
%
% Outputs are handles for the figure and subplots.
%
% Written by Richard Selwyn Jones, Durham University
% Modified from code by Greg Balco, Berkeley Geochronology Center, that was
% published in Schaefer et al., Nature, 2016.
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite.
%
%
%%

function depth_plot = plot_concs_depth(sample_data,save_plot,concs_name,x_lim,y_lim)

  % Check inputs
  if (nargin < 1 || nargin > 5)
      error('plot_concs_depth has wrong number of inputs!');
  end
  if (nargin < 2 || isempty(save_plot))
      save_plot = 0;
  end
  if save_plot == 1 && (nargin < 3 || isempty(concs_name))
      error('plot_concs_depth requires concs_name to save the figure');
  end
  
  
  % Define colours
  d_colours = [0.86,0.098,0.106; 0.208,0.475,0.694; 0.291,0.66,0.279; 0.574,0.294,0.615; 0.962,0.479,0; 0.996,0.996,0.199]; % Dark
  l_colours = [0.98,0.703,0.68; 0.699,0.801,0.887; 0.797,0.918,0.77; 0.867,0.793,0.891; 0.992,0.848,0.648; 0.996,0.996,0.598]; % Light
  
  fig = figure;
  a1 = subplot(1,1,1);
  a1.Position(2) = 0.13;
  
  
  % Plot sample concentrations with depth
  axes(a1);
  for a = 1:length(sample_data.s)
      this_s = sample_data.s{a};
      
      if this_s.nuclide10 == 1
          d_col = d_colours(1,:);
          l_col = l_colours(1,:);
          
          % Plot mean
          xx = [this_s.N10 this_s.N10];
          yy = [max(this_s.bottom_z)/100 min(this_s.top_z)/100];
          plot(xx,yy,'Color',l_col,'LineWidth',1);
          hold on;
          
          % Plot 1-sigma box
          N10_low = this_s.N10 - this_s.dN10;
          N10_upp = this_s.N10 + this_s.dN10;
          xx = [N10_low N10_low N10_upp N10_upp N10_low];
          for b = 1:length(this_s.top_z)
              yy = [this_s.top_z(b)/100 this_s.bottom_z(b)/100 this_s.bottom_z(b)/100 this_s.top_z(b)/100 this_s.top_z(b)/100];
              plot(xx,yy,'Color',d_col,'LineWidth',1);
          end
      end
      
      if this_s.nuclide26 == 1
          d_col = d_colours(2,:);
          l_col = l_colours(2,:);
          
          % Plot mean
          xx = [this_s.N26 this_s.N26];
          yy = [max(this_s.bottom_z)/100 min(this_s.top_z)/100];
          plot(xx,yy,'Color',l_col,'LineWidth',1);
          hold on;
          
          % Plot 1-sigma box
          N26_low = this_s.N26 - this_s.dN26;
          N26_upp = this_s.N26 + this_s.dN26;
          xx = [N26_low N26_low N26_upp N26_upp N26_low];
          for b = 1:length(this_s.top_z)
              yy = [this_s.top_z(b)/100 this_s.bottom_z(b)/100 this_s.bottom_z(b)/100 this_s.top_z(b)/100 this_s.top_z(b)/100];
              plot(xx,yy,'Color',d_col,'LineWidth',1);
          end
      end
  end
  
  
  % Adjust axes
  if (nargin == 5 && ~isempty(y_lim))
      a1.YLim = y_lim;
  else
      a1.YLim(1) = 0; % Default - set top depth to zero (i.e. the surface)
  end
  if (nargin == 4 && ~isempty(x_lim))
      a1.XLim = x_lim;
  end
  a1.YDir = 'reverse';
  xlabel('Nuclide concentration (atoms g^{-1})');
  ylabel('Depth below surface (m)');
  grid on;
  a1.GridLineStyle = ':';
  
  
  % Save figure
  if save_plot == 1
      scaling_model = sample_data.scaling_model;
      if strcmpi(scaling_model,'SF')
          scaling_model = 'LSD';
      elseif strcmpi(scaling_model,'SA')
          scaling_model = 'LSDn';
      end
      fig_name = strcat(concs_name,'_conc-vs-depth_',scaling_model);
      export_fig(fig_name,'-png','-r300','-transparent'); % Save as a PNG (raster) file
      saveas(gcf,fig_name,'epsc'); % Save as a EPS (vector) file
  end
  
  % Export handles
  depth_plot = [fig a1];  

end
