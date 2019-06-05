%
% add_names_2iso(norm_sample_data)
%
% Adds sample names to a two-isotope plot (26Al/10Be vs 10Be).
%
% norm_sample_data is a required struct, initially created using get_data.m
% with concentrations normalised in plot_2iso_concs.m.
%
% Written by Richard Selwyn Jones, Durham University
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite.
%
%
%%

function add_names_2iso(norm_sample_data)
  
  if nnz(norm_sample_data.logical_10) ~= nnz(norm_sample_data.logical_26)
      
      warning('Names can only be added where samples have both Be-10 and Al-26 measurements.')
  
  else
          
      % Get data for each sample
      names = cell(length(norm_sample_data.s),1);
      N10 = zeros(length(norm_sample_data.s),1);
      ratio = zeros(length(norm_sample_data.s),1);
      for a = 1:length(norm_sample_data.s)
          this_sample = norm_sample_data.s{a};
          names{a} = this_sample.name{1};
          N10(a) = this_sample.norm_N10;
          ratio(a) = this_sample.norm_N26/N10(a);
      end
      N10_5percent = (N10/100)*5;
      
      % Add text to plot
      h = text(N10+N10_5percent,ratio,char(names));
      set(h,'FontSize',10);
  
  end
  
end
