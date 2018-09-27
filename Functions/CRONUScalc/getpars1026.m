%
%  [pp,sp1026,sf1026,cp1026]=getpars1026(sampledata1026,scaling_model,maxdepth);  
%
%
% The sampledata1026 vector contains the following information:
%
%1. Latitude (decimal degrees)
%2. Longitude (decimal degrees)
%3. Elevation (meters)
%4. Pressure (hPa)       (Exactly one of 3 or 4 must be NaN.)
%5. sample thickness (cm)
%6. bulk density (g/cm^3)
%7. Shielding factor for terrain, snow, etc. (unitless)
%8. Erosion-rate epsilon (mm/kyr)
%9. Sample 10-Be concentration (atoms of 10-Be/g of target)
%10. Sample 26-Al concentration (atoms of 26-Al/g of target)
%11. Inheritance for Be (atoms 10-Be/g of target)
%12. Inheritance for Al (atoms 26-Al/g of target)
%13. Lambdafe Effective neutron attenuation length (g/cm^2)
%14. Depth to top of sample (g/cm^2)
%
% The maxdepth parameter is optional if it is given, it specifies
% the maximum depth (in g/cm^2) for which production rates will be
% computed.  If not given, the default value (currently 2500
% g/cm^2) is used instead. 
%
function [pp,sp1026,sf1026,cp1026]=getpars1026(sampledata1026,scaling_model,maxdepth)
%
% If maxdepth is not specified, use a default.
%
if (nargin < 3)
  maxdepth=2500;
end
%
% Get the physical parameters.
%
pp=physpars();
%
% Extract the sample parameters from the sampledatavector.
%
sp1026=samppars1026(sampledata1026);
%
% Get the scale factors.
%
sf1026=scalefacs1026(sp1026);
%
% Computed parameters.
%
cp1026=comppars1026(pp,sp1026,sf1026,maxdepth);
%
% Go ahead and produce contemporary scaling factors.
%
sf1026.currentsf=getcurrentsf(sf1026, 0, scaling_model, 'albe');


