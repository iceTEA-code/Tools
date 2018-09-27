%
% core_plot = plot_core_concs(sample_data,plot_ts)
% core_plot = plot_core_concs(sample_data,plot_ts,x_lim,y_lim)
%
% Plots the nuclide concentrations vs. depth down a core.
%
% sample_data is a required struct, created using get_data.m.
%
% plot_ts is a binary input for the plot type - '0' to plot just nuclide
% concentrations, or '1' to additionally plot the time-series data and
% corresponding periods of exposure/burial.
%
% x_lim and y_lim set the limits for the x-axis and y-axis ([lower upper]).
% If y_lim is empty, then the uppermost depth is set to zero (the surface).
%
% Can be used with plot_time_series.m and plot_ExpBur.m to subsequently
% plot the time-series data corresponding periods of exposure and burial as
% part of the same figure.
%
% Outputs are handles for the figure and subplots.
%
%
%%

function core_plot = plot_core_concs(sample_data,plot_ts,x_lim,y_lim)

  % Check inputs
  if (nargin < 2 || nargin > 4)
      error('plot_core_concs has wrong number of inputs!');
  end
  if plot_ts ~= 1 && plot_ts ~= 0
      error('plot_ts should be binary (1 or 0)!');
  end
  
  
  % Define colours
  d_colours = [0.86,0.098,0.106; 0.208,0.475,0.694; 0.291,0.66,0.279; 0.574,0.294,0.615; 0.962,0.479,0; 0.996,0.996,0.199]; % Dark
  l_colours = [0.98,0.703,0.68; 0.699,0.801,0.887; 0.797,0.918,0.77; 0.867,0.793,0.891; 0.992,0.848,0.648; 0.996,0.996,0.598]; % Light
  
  fig = figure;
  if plot_ts == 1
      a1 = subplot(3,1,1); a2 = subplot(3,1,2); a3 = subplot(3,1,3);
      set(a1,'pos',[0.14 0.40 0.78 0.54]);
      set(a2,'pos',[0.14 0.15 0.78 0.12]);
      set(a3,'pos',[0.14 0.11 0.78 0.03]);
  else
      a1 = subplot(1,1,1);
      a1.Position(2) = 0.13;
  end
  
  
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
          plot(xx,yy,'Color',l_col);
          hold on;
          
          % Plot 1-sigma box
          N10_low = this_s.N10 - this_s.dN10;
          N10_upp = this_s.N10 + this_s.dN10;
          xx = [N10_low N10_low N10_upp N10_upp N10_low];
          for b = 1:length(this_s.top_z)
              yy = [this_s.top_z(b)/100 this_s.bottom_z(b)/100 this_s.bottom_z(b)/100 this_s.top_z(b)/100 this_s.top_z(b)/100];
              plot(xx,yy,'Color',d_col);
          end
      end
      
      if this_s.nuclide26 == 1
          d_col = d_colours(2,:);
          l_col = l_colours(2,:);
          
          % Plot mean
          xx = [this_s.N26 this_s.N26];
          yy = [max(this_s.bottom_z)/100 min(this_s.top_z)/100];
          plot(xx,yy,'Color',l_col);
          hold on;
          
          % Plot 1-sigma box
          N26_low = this_s.N26 - this_s.dN26;
          N26_upp = this_s.N26 + this_s.dN26;
          xx = [N26_low N26_low N26_upp N26_upp N26_low];
          for b = 1:length(this_s.top_z)
              yy = [this_s.top_z(b)/100 this_s.bottom_z(b)/100 this_s.bottom_z(b)/100 this_s.top_z(b)/100 this_s.top_z(b)/100];
              plot(xx,yy,'Color',d_col);
          end
      end
  end
  
  
  % Adjust axes
  if (nargin == 4 && ~isempty(y_lim))
      a1.YLim = y_lim;
  else
      a1.YLim(1) = 0; % Default - set top depth to zero (i.e. the surface)
  end
  if (nargin == 3 && ~isempty(x_lim))
      a1.XLim = x_lim;
  end
  a1.YDir = 'reverse';
  xlabel('Nuclide concentration (atoms g^{-1})');
  ylabel('Core depth (m)');
  grid on;
  a1.GridLineStyle = ':';
  
  
  % Export handles
  if plot_ts == 1
      core_plot = [fig a1 a2 a3];
  else
      core_plot = [fig a1];
  end
  

end
