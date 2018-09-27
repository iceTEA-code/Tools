%
% get_version()
% 
% Returns the version of iceTEA being used.
% 
% No input is required.
%
% Written by Richard Selwyn Jones, Durham University
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite.
%
%
%%

function get_version()
  
  vers_file = importdata('iceTEA_version.txt');
  version = vers_file.data;
  
  disp(['iceTEA version = v' num2str(version)])

end
