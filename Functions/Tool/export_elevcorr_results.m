%
% export_elevcorr_results(sample_data,corrected_data,out_name,format)
%
% Exports the elevation-corrected data to either an Excel spreadsheet or a 
% tab-deliminated .txt file.
%
% sample_data is a struct containing necessary sample information produced 
% using get_data.m.
%
% corrected_data is a struct containing elevation-correction information 
% produced using elev_correct.m.
%
% out_name corresponds to the name of dataset that will be used to save the
% corrected data (e.g. *_results.xlsx).
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

function export_elevcorr_results(sample_data,corrected_data,out_name,format)
  
  % Check inputs
  if (nargin ~= 4)
      error('export_elevcorr_results has wrong number of inputs!');
  end
    
  % Find nuclides
  logical_10 = sample_data.logical_10;
  logical_26 = sample_data.logical_26;
  
  % Get sample details
  latitudes = sample_data.CC(:,1);
  longitudes = sample_data.CC(:,2);
  mod_elevations = sample_data.CC(:,3);
  for i = 1:length(sample_data.s)
      sample_names{i} = sample_data.s{1,i}.name;
  end
  
  % Get correction info
  scaling_model = corrected_data.scaling_model;
  correction_type = corrected_data.correction_type;
  elev_input = corrected_data.elev_input;
  
  % Create table header
  if strcmp(correction_type,'rate')
      correction_type = 'Rate';
      xls_header = ["Sample name","Latitude (DD)","Longitude (DD)","Correction type","Rate of elevation change (m/ka)","Modern elevation (m asl)","Average elevation since exposed (m above present sea level))","Nuclide","Uncorrected age (mean; yrs)","Corrected age (mean; yrs)","Corrected age (1 sig.; yrs)","Mean age difference (yrs)","Mean age difference (%)","Scaling model"];
      txt_header = ["Sample_name","Latitude_DD","Longitude_DD","Correction_type","Rate_of_elevation_change_mka","Modern_elevation_masl","Average_elevation_since_exposed_masl","Nuclide","Uncorrected_age_mean_yrs","Corrected_age_mean_yrs","Corrected_age_1sig_yrs","Mean_age_difference_yrs","Mean_age_difference_percent","Scaling_model"];
  elseif strcmp(correction_type,'model')
      correction_type = 'GIA model';
      xls_header = ["Sample name","Latitude (DD)","Longitude (DD)","Correction type","Ice model","Modern elevation (m asl)","Average elevation since exposed (m above present sea level)","Nuclide","Uncorrected age (mean; yrs)","Corrected age (mean; yrs)","Corrected age (1 sig.; yrs)","Mean age difference (yrs)","Mean age difference (%)","Scaling model"];
      txt_header = ["Sample_name","Latitude_DD","Longitude_DD","Correction_type","Ice_model","Modern_elevation_masl","Average_elevation_since_exposed_masl","Nuclide","Uncorrected_age_mean_yrs","Corrected_age_mean_yrs","Corrected_age_1sig_yrs","Mean_age_difference_yrs","Mean_age_difference_percent","Scaling_model"];
  end
    
  % Sort data for each nuclide
  if any(logical_10)
      Be_data = string(zeros(nnz(logical_10),14));
      nuclide = 'Be-10';
      for a = 1:nnz(logical_10)
          Be_data(a,8) = nuclide;
      end
      Be_data(:,1) = sample_names(logical_10);
      Be_data(:,2) = round(latitudes(logical_10),4);
      Be_data(:,3) = round(longitudes(logical_10),4);
      Be_data(:,6) = round(mod_elevations(logical_10),0);
      Be_data(:,7) = round(corrected_data.N10.mean_elev,0);
      Be_data(:,9) = round(corrected_data.N10.uncorr_mean_age,-1);
      Be_data(:,10) = round(corrected_data.N10.corr_mean_age,-1);
      Be_data(:,11) = round(corrected_data.N10.corr_age_ext,-1);
      Be_data(:,12) = round(corrected_data.N10.mean_age_diff,-1);
      Be_data(:,13) = round(corrected_data.N10.per_age_diff,0);
  end
  if any(logical_26)
      Al_data = string(zeros(nnz(logical_26),14));
      nuclide = 'Al-26';
      for a = 1:nnz(logical_26)
          Al_data(a,8) = nuclide;
      end
      Al_data(:,1) = sample_names(logical_26);
      Al_data(:,2) = round(latitudes(logical_26),4);
      Al_data(:,3) = round(longitudes(logical_26),4);
      Al_data(:,6) = round(mod_elevations(logical_26),0);
      Al_data(:,7) = round(corrected_data.N26.mean_elev,0);
      Al_data(:,9) = round(corrected_data.N26.uncorr_mean_age,-1);
      Al_data(:,10) = round(corrected_data.N26.corr_mean_age,-1);
      Al_data(:,11) = round(corrected_data.N26.corr_age_ext,-1);
      Al_data(:,12) = round(corrected_data.N26.mean_age_diff,-1);
      Al_data(:,13) = round(corrected_data.N26.per_age_diff,0);
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
      out_data(c,4) = correction_type;
      out_data(c,5) = elev_input;
      out_data(c,14) = scaling_model;
  end
  
  
  % Export
  out_table = table(out_data(:,1),out_data(:,2),out_data(:,3),out_data(:,4),out_data(:,5),out_data(:,6),out_data(:,7),out_data(:,8),out_data(:,9),out_data(:,10),out_data(:,11),out_data(:,12),out_data(:,13),out_data(:,14));
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
