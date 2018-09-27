function l=rawattenuationlength(pressure,rigiditycutoff)
%
% l=attenuationlength(pressure,rigiditycutoff)
%
% Interpolates an attenuation length for a given atmospheric depth
% (in g/cm^2) and rigidity cutoff.  The attenuation length returned
% is the "flat horizon" Lambdafe.  Further correction may be
% required to account for topography.  
%

%
% Take care of extreme rigidity cutoffs that should really be 0.
%
if (rigiditycutoff < 0)
  rigiditycutoff=0;
end
%
% Convert pressure (hPa) to atmospheric depth (g/cm^2).
%
atmdepth=pressure*1.019716;
%
% Sato Table.
%
T=[151.4 151.8 152.8 153.7 155.2 157.8 162.5 171.8;...
   152.1 152.4 154.1 156.1 159.1 163.7 171.2 185.0;...
   155.5 158.6 160.6 163.7 168.1 174.6 184.9 203.6;...
   162.0 164.8 167.7 171.2 176.3 183.7 195.7 217.9;...
   167.8 170.0 172.4 177.2 182.4 190.6 203.6 228.2;...
   168.8 170.9 174.5 178.6 184.0 192.3 205.8 230.9];
depths=[1100 1000 900 800 700 600 500 400];
cutoffs=[0; 4; 8; 12; 16; 20];
%
% interpolate from the table.
%
l=interp2(depths,cutoffs,T,atmdepth,rigiditycutoff);
%
% Check for NaN.
%
if (isnan(l))
  %warning('Attenuation length returned NaN. Using a value of 160.');
  l = 160;
end