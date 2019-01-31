%
% export_calcages(sample_data,ages_ka,out_name,format)
%
% Exports the calculated exposure age data to either an Excel spreadsheet 
% or a tab-deliminated .txt file.
%
% sample_data is a struct containing necessary sample information produced 
% using get_data.m.
%
% ages_ka is a struct containing calculated exposure ages produced using 
% age_calc.m.
%
% out_name corresponds to the name of dataset that will be used to save the
% exposure age data (e.g. *_results.xlsx).
%
% format should be either "xls" or "txt" to export a .xlsx or .txt file.
%
% Written by Richard Selwyn Jones, Durham University
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite, which is built on versions of 
% CRONUS-Earth and CRONUScalc code.
%
%
%%

function export_calcages(sample_data,ages_ka,out_name,format)
  
  % Check inputs
  if (nargin ~= 4)
      error('export_calcages has wrong number of inputs!');
  end
    
  % Find nuclides
  logical_10 = sample_data.logical_10;
  logical_26 = sample_data.logical_26;
  
  % Get sample details
%   latitudes = sample_data.CC(:,1);
%   longitudes = sample_data.CC(:,2);
%   elevations = sample_data.CC(:,3);
  sample_names = ages_ka.sample_names;
  
  % Get calculation info
  scaling_model = ages_ka.scaling_model;
  
  % Create table header
  xls_header = ["Sample name","Nuclide","Exposure age (mean; yrs)","Measurement uncertainty (1 sig.; yrs)","Total uncertainty (1 sig.; yrs)","Scaling model"];
  txt_header = ["Sample_name","Nuclide","Exposure_age_mean_yrs","Measurement_uncertainty_1sig_yrs","Total_uncertainty_1sig_yrs","Scaling_model"];

    
  % Sort data for each nuclide
  if any(logical_10)
      Be_data = string(zeros(nnz(logical_10),6));
      nuclide = 'Be-10';
      for a = 1:nnz(logical_10)
          Be_data(a,2) = nuclide;
      end
      Be_data(:,1) = sample_names(logical_10);
      Be_data(:,3) = round(ages_ka.Be10(:,1)*1000,-1);
      Be_data(:,4) = round(ages_ka.Be10(:,2)*1000,-1);
      Be_data(:,5) = round(ages_ka.Be10(:,3)*1000,-1);
  end
  if any(logical_26)
      Al_data = string(zeros(nnz(logical_26),6));
      nuclide = 'Al-26';
      for a = 1:nnz(logical_26)
          Al_data(a,2) = nuclide;
      end
      Al_data(:,1) = sample_names(logical_26);
      Al_data(:,3) = round(ages_ka.Al26(:,1)*1000,-1);
      Al_data(:,4) = round(ages_ka.Al26(:,2)*1000,-1);
      Al_data(:,5) = round(ages_ka.Al26(:,3)*1000,-1);
  end
  
  if any(logical_10) && ~any(logical_26)
      out_data = Be_data;
  elseif ~any(logical_10) && any(logical_26)
      out_data = Al_data;
  else
      out_data = [Be_data; Al_data];
  end
  
  % Add universal info
  for c = 1:length(out_data(:,1))
      out_data(c,6) = scaling_model;
  end
  
  
  % Export
  out_table = table(out_data(:,1),out_data(:,2),out_data(:,3),out_data(:,4),out_data(:,5),out_data(:,6));
  out_table.Properties.VariableNames = cellstr(txt_header);
  
  if strcmpi(format,'xls') || strcmpi(format,'xlsx') || strcmpi(format,'xsl') || strcmpi(format,'xslx')
      
      % Export table to .xlsx file
      save_name = strcat(out_name,'_results.xlsx');
      if ispc
          xlswrite(save_name,[xls_header; out_data]);
      else
          writetable(out_table,save_name);
      end
      
  elseif strcmpi(format,'txt') || strcmpi(format,'text')
      
      % Export table to .txt file
      save_name = strcat(out_name,'_results.txt');
      writetable(out_table,save_name,'Delimiter','tab');
      
  else
      error('format should be "xls" (Excel spreadsheet) or "txt" (tab-deliminated .txt file)!');
  end
  
end
