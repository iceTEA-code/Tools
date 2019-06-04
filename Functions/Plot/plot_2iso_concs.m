%
% iso_plot = plot_2iso_concs(sample_data,sigma,add_names,save_plot)
% iso_plot = plot_2iso_concs(sample_data,sigma,add_names,save_plot,concs_name)
% iso_plot = plot_2iso_concs(sample_data,sigma,add_names,save_plot,concs_name,...,expo_intervals,bur_intervals)
% iso_plot = plot_2iso_concs(sample_data,sigma,add_names,save_plot,concs_name,...,x_lim,y_lim)
%
% Plots normalised nuclide concentrations on a two-isotope plot (26Al/10Be
% vs 10Be).
%
% The steady-state erosion island is plotted as a black line, while 
% exposure and burial isochrons are plotted as grey dashed and dot-dashed
% lines, respectively.
%
% sample_data is a required struct, created using get_data.m.
%
% sigma should be either 1 or 2. If '1' then only 1 sigma concentration
% uncertainty is plotted, and if '2' then both 1 and 2 sigma uncertainties
% are shown.
%
% add_names is a binary input to show sample names on the plot [1], or not
% [0]. Names can only be added where samples have both Be-10 and Al-26 
% measurements
%
% save_plot is a binary input to save the figure(s) in .png and .eps 
% formats [1], or not [0].
%
% concs_name is string input for the name of the data, which will be used
% to save the figure. It is required if save_plot is 1.
%
% expo_intervals and bur_intervals can be vectors of values in unit ka,
% which corresponds to the plotted and labelled isochrons. If empty, then
% defaults are used.
%
% x_lim and y_lim set the limits for the x-axis and y-axis ([lower upper]).
% If empty, then defaults are used.
%
% The figure (fig) and axes (ax) handles are exported as an output.
%
% Written by Richard Selwyn Jones, Durham University
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite.
%
%
%%

function [iso_plot,sample_data] = plot_2iso_concs(sample_data,sigma,add_names,save_plot,concs_name,expo_intervals,bur_intervals,x_lim,y_lim)
  
  % Check inputs
  if (nargin < 2 || nargin > 9)
      error('plot_2iso_concs has wrong number of inputs!');
  end
  
  if (~sigma == 1 || ~sigma == 2)
      error('sigma should be 1 or 2!');
  end
  if (nargin < 3 || isempty(add_names))
      add_names = 0;
  end
  if (nargin < 4 || isempty(save_plot))
      save_plot = 0;
  end
  if save_plot == 1 && (nargin < 5 || isempty(concs_name))
      error('plot_2iso_concs requires concs_name to save the figure');
  end
  
  if (nargin > 5  && ~isempty(expo_intervals))
      expo_int = expo_intervals;
  else
      expo_int = [5,10,25,50,100,250,500,1000,2000,10000]; % Default exposure intervals to be shown (ka)
  end
  if (nargin > 5 && ~isempty(bur_intervals))
      bur_int = bur_intervals;
  else
      bur_int = 300:300:1500; % Default burial intervals to be shown (ka)
  end
  
  
  % Compute production rates for samples and normalise concentrations
  sample_data = norm_concs(sample_data);
  
   
  % Get average parameters for samples
  for b = 1:length(sample_data.sp1026)
      s_rho(b) = sample_data.sp1026{b}.rb;
      s_L(b) = sample_data.sp1026{b}.Lambdafe;
  end
  rho = mean(s_rho);
  L = mean(s_L);
  
   
  % Continually exposed (non-eroding)
  A_10 = sample_data.pp.lambda10Be + rho / L * 0; % Zero erosion
  A_26 = sample_data.pp.lambda26Al + rho / L * 0; % Zero erosion
  t_expos = linspace(0,10e6,10000+1); % Time exposed (a) - 1a to 10Ma
  
  be_conc_expo0ero = (1./A_10) .* (1-exp(-1*A_10*t_expos));
  al_conc_expo0ero = (1./A_26) .* (1-exp(-1*A_26*t_expos));
  
  ratio_expo0ero = al_conc_expo0ero ./ be_conc_expo0ero;
  
  % Continually exposed (steady-state erosion line)
  e_subaerial = linspace(0.1,0.0000001,20000); % Sub-aerial erosion (cm/a)
  A_10 = sample_data.pp.lambda10Be + rho / L * e_subaerial;
  A_26 = sample_data.pp.lambda26Al + rho / L * e_subaerial;
  t_expos = 5e7; % Constant time exposed (50 Ma)
  
  be_conc_expoSSero = (1./A_10) .* (1-exp(-1*A_10*t_expos));
  al_conc_expoSSero = (1./A_26) .* (1-exp(-1*A_26*t_expos));
  
  ratio_expoSSero = al_conc_expoSSero ./ be_conc_expoSSero;
  
  % Burial isochrons
  A_10 = sample_data.pp.lambda10Be + rho / L * 0; % Zero erosion
  A_26 = sample_data.pp.lambda26Al + rho / L * 0; % Zero erosion

  t_expos = linspace(0,10e6,10000+1); % Time exposed (a) - 1a to 10Ma
  t_bur = linspace(0,1e7,200+1); % Burial time (a) - 1a to 100Ma
  
  for k=1:length(t_expos) % Time consuming part
      be_conc(k) = (1./A_10) .* (1-exp(-1*A_10*t_expos(k)));
      al_conc(k) = (1./A_26) .* (1-exp(-1*A_26*t_expos(k)));
      for i=1:length(t_bur)
          be_conc_burial(k,i) = be_conc(k)*exp(-1*A_10*t_bur(i));
          al_conc_burial(k,i) = al_conc(k)*exp(-1*A_26*t_bur(i));
      end
  end
  
  ratio_burial = al_conc_burial ./ be_conc_burial;
  
  
  % Determine isochron grid
  labels_expo0ero = strcat(num2str(expo_int'),' ka');
  expo0ero_col = zeros(length(expo_int),1);
  for c = 1:length(expo_int)
      expo0ero_col(c) = find(t_expos==(expo_int(c)*1000));
  end
  labels_t_bur = strcat(num2str(bur_int'),' ka');
  burial_col = zeros(length(bur_int),1);
  for c = 1:length(bur_int)
      burial_col(c) = find(t_bur==(bur_int(c)*1000));
  end
  
  
  % Sort nuclide data to plot; Mix if necessary
  plot_N = N_to_plot(sample_data);
  
  
  % PLOT
  fig = figure;
  a1 = subplot(1,1,1);
  a1.Position(2) = 0.13;

  % Make colours
  grey = [.5 .5 .5];
  %colour = [0.86,0.098,0.106; 0.98,0.703,0.68]; % Sample ellipses = Red
  %colour = [0.48,0.20,0.58; 0.76,0.64,0.81]; % Sample ellipses = Purple
  colour = [0,0.53,0.216; 0.65,0.858,0.628]; % Sample ellipses = Green

  axes(a1);
  semilogx(be_conc_burial(:,burial_col),ratio_burial(:,burial_col),'-.','Color',grey,'Linewidth',1.2);
  hold on;
  semilogx(be_conc_burial(expo0ero_col,:)',ratio_burial(expo0ero_col,:)','--','Color',grey,'Linewidth',1.2);
  semilogx(be_conc_expo0ero,ratio_expo0ero,'-k','Linewidth',1.4);
  semilogx(be_conc_expoSSero,ratio_expoSSero,'-','Color','k','Linewidth',1.4);
  text(be_conc_expo0ero(1,expo0ero_col),ratio_expo0ero(1,expo0ero_col),labels_expo0ero,'VerticalAlignment','bottom','HorizontalAlignment','left','Color',grey)
  text(be_conc_burial(expo0ero_col(2)+1,burial_col),ratio_burial(expo0ero_col(2)+1,burial_col),labels_t_bur,'VerticalAlignment','bottom','HorizontalAlignment','left','Color',grey)
  
  % Plot sample ellipses
  for el = 1:length(plot_N.norm_N10)
      N10 = plot_N.norm_N10(el);
      dN10 = plot_N.norm_dN10(el);
      N26 = plot_N.norm_N26(el);
      dN26 = plot_N.norm_dN26(el);
      if sigma == 1
          h_e = ellipse(N10,dN10,N26,dN26,1);
          set(h_e,'Color',colour(1,:),'LineWidth',1,'LineStyle','-');
          plot(N10,(N26./N10),'.','Color',colour(1,:));
      else
          [h_e1,h_e2] = ellipse(N10,dN10,N26,dN26,2);
          set(h_e1,'Color',colour(1,:),'LineWidth',1,'LineStyle','-');
          set(h_e2,'Color',colour(2,:),'LineWidth',1,'LineStyle','-');
          plot(N10,(N26./N10),'.','Color',colour(1,:));
      end
  end
 
  hold off;
  
  
  % Adjust axes
  if ((nargin == 7 || nargin == 9) && ~isempty(y_lim))
      ylim(y_lim);
  else
      if a1.YLim(1) < 0.4 && a1.YLim(2) > 1
          ylim(a1.YLim);
      elseif a1.YLim(1) > 0.4 && a1.YLim(2) > 1
          ylim([0.4,a1.YLim(2)]);
      elseif a1.YLim(1) < 0.4 && a1.YLim(2) < 1
          ylim([a1.YLim(1),1.05]);
      else
          ylim([0.4,1.05]);
      end
  end
  if ((nargin == 7 || nargin == 9) && ~isempty(x_lim))
      xlim(x_lim);
  else
      xlim([2e3,3e6]);
  end
  set(a1,'xscale','log')
  xlabel('^{10}Be ^* (years)')
  ylabel('^{26}Al ^* / ^{10}Be ^*')
  box on;
  
  
  % Add sample names to plot
  if add_names == 1
     add_names_2iso(sample_data); 
  end
  
  
  % Save figure
  if save_plot == 1
      scaling_model = sample_data.scaling_model;
      if strcmpi(scaling_model,'SF')
          scaling_model = 'LSD';
      elseif strcmpi(scaling_model,'SA')
          scaling_model = 'LSDn';
      end
      fig_name = strcat(concs_name,'_2isotope_',scaling_model);
      export_fig(fig_name,'-png','-r300','-transparent'); % Save as a PNG (raster) file
      saveas(gcf,fig_name,'epsc'); % Save as a EPS (vector) file
  end
    
  % Export figure handles
  iso_plot = [fig a1];
   
  
end
