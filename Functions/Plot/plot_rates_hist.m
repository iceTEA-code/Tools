%
% out = plot_rates_hist(ages_ka,reg,transect_type,N,save_plot)
%
% Plots histograms of the non-negative thinning/retreat rates, and 
% identifies those rates at 68% and 95% confidence. It is designed to be 
% used with transect_regress.m.
%
% ages_ka is a required struct, containing ages and uncertainties 
% calculated using age_calc.m, elevations and positions, and a logical of 
% the nuclides that were measured.
%
% reg is a struct of necessary information computed in transect_regress.m.
%
% transect_type should be either "vert" or "horiz", which determines
% whether the transect is vertical (e.g. elevation above the modern ice) or
% horizontal (e.g. distance from modern ice terminus). The units are m/yr
% and km/yr respectively.
%
% N is the nuclide of the exposure age data (e.g. '10' for Be-10).
%
% save_plot is a binary input to save the figure(s) in .png and .eps 
% formats [1], or not [0].
%
% Outputs the figure (fig) and axes (ax) handles.
%
% Written by Richard Selwyn Jones, Durham University
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite.
%
%
%%

function out = plot_rates_hist(ages_ka,reg,transect_type,N,save_plot)

  % Check inputs
  if (nargin ~= 5)
      error('plot_rates_prob has wrong number of inputs!');
  end
  if (~strcmp(transect_type,'vert') && ~strcmp(transect_type,'horiz'))
      error('transect_type should be "vert" or "horiz"!');
  end
  if (N ~= 10 && N ~= 26)
      error('N should be "10" or "26"!');
  end
  
  % Suppress warning
  warning('off','stats:ksdensity:NoConvergence');
  
  % Define colours
  d_colours = [0.86,0.098,0.106; 0.208,0.475,0.694; 0.291,0.66,0.279; 0.574,0.294,0.615; 0.962,0.479,0];
  l_colours = [0.98,0.703,0.68; 0.699,0.801,0.887; 0.797,0.918,0.77; 0.867,0.793,0.891; 0.992,0.848,0.648];

  
  if N == 10
      N_name = 'Be-10';
      col = [d_colours(1,:); l_colours(1,:)];
  elseif N == 26
      N_name = 'Al-26';
      col = [d_colours(2,:); l_colours(2,:)];
  end
  
  % Generate bins
  max_rate = max(reg.rate); min_rate = min(reg.rate); rate_range = max_rate-min_rate;
  rate_95_range = reg.quant_95(2)-reg.quant_95(1);
  if rate_95_range > 5
      binlength = 0.25; % Bin width (m/yr)
  elseif rate_95_range < 1
      binlength = 0.01; % Bin width (m/yr)
  else
      binlength = 0.1; % Bin width (m/yr)
  end
  rate_bins = rate_range/binlength;
  nbins = round(rate_bins);

  rate_med = median(reg.rate);
  
  fig = figure;
  
  semilogx(0,0,'.k');
  hold on;
  
  h = histfit(reg.rate,nbins,'kernel');
  delete(h(2));
  set(h(1),'EdgeColor',col(1,:),'FaceColor',col(2,:));
  
  ax = gca;
  
  plot([rate_med rate_med],ax.YLim,'--','Color',col(1,:));
  
  hold off;
  
  mode_val = max(h(1).YData); % Get mode of distribution
  mod_ind = h(1).YData == mode_val; % Get mode index
  rate_mode = h(1).XData(mod_ind); % Find corresponding rate
  
  if (strcmp(transect_type,'vert'))
      %title('Distribution of (non-negative) random thinning rates');
      x_lab = 'Thinning rate (m yr^{-1})';
      summary_text = {strcat('{Nuclide: }',N_name),strcat('{Scaling model: }',ages_ka.scaling_model),strcat('{LS regression: }',reg.regress_type),strcat('{MC iterations: }',num2str(reg.n_iter)),strcat('{Mode: }',num2str(round(rate_mode,2)),'{ m yr^{-1}}'),strcat('{Median: }',num2str(round(rate_med,2)),'{ m yr^{-1}}'),strcat('{68%: }',num2str(round(reg.quant_68(1),2)),'-',num2str(round(reg.quant_68(2),2)),'{ m yr^{-1}}'),strcat('{95%: }',num2str(round(reg.quant_95(1),2)),'-',num2str(round(reg.quant_95(2),2)),'{ m yr^{-1}}')};
  else
      %title('Distribution of (non-negative) random retreat rates');
      x_lab = 'Retreat rate (km yr^{-1})';
      summary_text = {strcat('{Nuclide: }',N_name),strcat('{Scaling model: }',ages_ka.scaling_model),strcat('{LS regression: }',reg.regress_type),strcat('{MC iterations: }',num2str(reg.n_iter)),strcat('{Mode: }',num2str(round(rate_mode,2)),'{ km yr^{-1}}'),strcat('{Median: }',num2str(round(rate_med,2)),'{ km yr^{-1}}'),strcat('{68%: }',num2str(round(reg.quant_68(1),2)),'-',num2str(round(reg.quant_68(2),2)),'{ km yr^{-1}}'),strcat('{95%: }',num2str(round(reg.quant_95(1),2)),'-',num2str(round(reg.quant_95(2),2)),'{ km yr^{-1}}')};
  end
  
  box_dim = [0.55 0.42 0.2 0.2];
  annotation('textbox',box_dim,'String',summary_text,'VerticalAlignment','middle','FitBoxToText','on','Color',col(1,:),'FontSize',10,'FontWeight','bold');
  
  box on;
  ax.XScale = 'linear';
  if max_rate > 1
      xlim([0 round(reg.quant_95(2))]);
  end
  ylabel('Frequency');
  xlabel(x_lab);
  
  
  % Save figure
  if save_plot == 1
      fig_name = strcat(ages_ka.ages_name,'_TransRegress_hist_',N_name,'_',ages_ka.scaling_model);
      export_fig(fig_name,'-png','-r300','-transparent'); % Save as a PNG (raster) file
      saveas(gcf,fig_name,'epsc'); % Save as a EPS (vector) file
  end
  
  out = [fig ax];
  
end
