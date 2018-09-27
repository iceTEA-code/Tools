%
% sample_data = atm_pressure(sample_data)
% 
% Computes atmospheric pressure if not known (i.e. if set to zero in
% sampledata).
%
% It uses ERA-40 or, if in Antarctica, an elevation-pressure relationship.
%
% Written by Richard Selwyn Jones, Durham University
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite, which is built on versions of 
% CRONUS-Earth and CRONUScalc code.
%
%
%%

function sample_data = atm_pressure(lats,lons,sample_data)

for a = 1:length(lats)
    
    if sample_data(a,4) == 0
        
        % Computed from elevation-pressure relationship of Antarctic stations
        if lats(a) < -60
            sample_data(a,4) = antatm(sample_data(a,3));
            
            % Computed from ERA-40 reanalysis data
        else
            sample_data(a,4) = ERA40atm(lats(a),lons(a),sample_data(a,3));
        end
        
    end
    
end

end
