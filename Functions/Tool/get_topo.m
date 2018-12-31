%
% out = get_topo(azimuths,horizon_angles,dip_direction,dip_angle)
% 
% Calculates topographic shielding using a function from CRONUScalc 
% (Marrero et al., 2016), which is a modified version of the Balco et al. 
% (2008) skyline.m function.
%
% azimuths is a required array of measured azimuths for the horizon 
% (0-360°).
%
% horizon_angles is a required array of the measured angle to the horizon 
% for each azimuth (0-90°).
%
% dip_direction is an optional input to specify the dip direction of a rock
% surface (0-360°).
%
% dip_angle is an optional input to specify the dip angle of the surface 
% (0-90°).
%
% Outputs shielding factors for topography and surface geometry, production
% change due to effective attenuation length, and total shielding.
%
% Written by Richard Selwyn Jones, Durham University
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite, which is built on versions of 
% CRONUS-Earth and CRONUScalc code.
%
%
%%

function out = get_topo(azimuths,horizon_angles,dip_direction,dip_angle)
  
  % Check inputs
  if (nargin < 2 || nargin > 4)
      error('get_topo has wrong number of inputs!');
  end
  if (nargin < 3) || isempty(dip_direction)
      dip_direction = 0;
      dip_angle = 0;
  end
  
  theta_horiz = [azimuths',horizon_angles'];
  topo_out = topofactor(theta_horiz,dip_direction,dip_angle);
  
  disp(['Shielding factor: ' num2str(topo_out.ST_calc)])
  
  % Export
  out.TopoShielding = topo_out.Stopog;
  out.AttenShielding = topo_out.Slambda;
  out.TotalShielding = topo_out.ST_calc;  
  
end
