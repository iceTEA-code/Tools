%
% cover_shield_fac = get_cov_shield(cover_type,cover_depth,attenuation,cover_density)
%
% Computes a surface cover shielding factor for correcting exposure ages.
% Either a preset cover type can be selected or a manual cover density
% can be used.
%
% cover_type should be 'snow', 'freshwater', 'seawater', 'loess', 'till', 
% 'soil', 'ash', or 'manual'.
%
% cover_depth is a required value of the depth of surface cover (in cm).
%
% atten_L is a required value for the apparent attenuation length 
% (in g/cm^2).
%
% cover_density is an optional input that is to be used if the cover_type
% is set to "manual". A value should be in g/cm^3.
%
% Written by Richard Selwyn Jones, Durham University
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite.
%
%
%%

function cover_shield_fac = get_cov_shield(cover_type,cover_depth,atten_L,cover_density)
  
  % Check inputs
  if (nargin < 3 || nargin > 4)
      error('get_cov_shield has wrong number of inputs!');
  end
  if (nargin == 3) && strcmpi(cover_type,'manual')
      error('cover_type is ''manual'', but density input is missing.');
  elseif ~isempty(cover_density) && ~strcmpi(cover_type,'manual')
      error('Density input seems to be included, but cover_type is not ''manual''.');
  end

  % Get density of surface cover (g/cm^3)
  if strcmpi(cover_type,'snow')
      density = 0.27;
  
  elseif strcmpi(cover_type,'freshwater') || strcmpi(cover_type,'fresh water')
      density = 0.999; % Fresh water near surface (1.1 bars) at 10°C.
      
  elseif strcmpi(cover_type,'seawater') || strcmpi(cover_type,'sea water')
      density = 1.027; % Seawater near surface (1.1 bars) at 10°C with salinity of 35 g/kg.
  
  elseif strcmpi(cover_type,'soil')
      density = 1.3; % Average of dry mineral soil (~1.0-1.6 g/cm^3)

  elseif strcmpi(cover_type,'loess')
      density = 1.6;
      
  elseif strcmpi(cover_type,'till')
      density = 1.8;
      
  elseif strcmpi(cover_type,'ash')
      density = 0.7;

  elseif strcmpi(cover_type,'manual')
      density = cover_density;
      
  else
      error('cover_type is not one of ''snow'', ''freshwater'', ''seawater'', ''loess'', ''till'', ''soil'', ''ash'', or ''manual''!')

  end
  
  
  % Calculate depth of cover (g/cm^2)
  depth_gcm2 = cover_depth*density;
  
  
  % Compute shielding factor
  cover_shield_fac = exp( - depth_gcm2 / atten_L );
  
end
