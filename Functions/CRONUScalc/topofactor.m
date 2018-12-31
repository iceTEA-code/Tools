function out=topofactor(input,DD,dip_angle,varargin)

% function out=topo(input,DD,dip_angle,options...)
%
% Duplicate of topography calculation in CHLOE
%
% Inputs:
%   input    - n by 2 array of Theta & Horizon values. e.g. [1 1;2 2;3 3]
%              If no Theta and Horizon measurements (e.g., dipping surface
%              only), just put one row of [0,0]
%   DD       - Direction of dip
%   dip_angle- dip angle
%   options  - Optional values for constants Ssnow Lf_t and Lf_a.
%              If given, they will override the default values.
%              Example usage: topo(input,40,25,'Ssnow',1,'Lf_a',180) 
%
% Returns:   A struct with fields:
%           .ST_calc : Shielding factor incorporating all effects
%           .Slambda : Relative production change due to effective 
%                      attenuation length
%           .Lf_calc : Effective fast neutron attenuation length
%           .Stopog  : Calculated topographic/surface geometry shielding factor
%



% There must be at least 3 input arguments
if(nargin < 3)
   error('Insufficient input arguments');
end

% input must be Nx2 array
if(size(input,2) ~= 2)
	error('Input array must be n by 2 array');
end

% Constants (copied from CHLOE). Can be overridden in options.
Ssnow = 1;

m = 2.3; % Angular dependence exponent - commonly used value not actually 
         % given anywhere in the literature

Lf_a = 152; %Apparent fast neutron attenuation coefficient from Sato in 
            %rawattenuationlength.m, at SLHL   
Lf_t = (m+2)/(m+1)*Lf_a; %True fast neutron attenuation coefficient

nvargs = size(varargin,2);
if(nvargs > 0)
   % Must be multiple of 2:
   if(mod(nvargs,2) ~= 0)
      error('Cannot parse optional arguments');
   end
   for k=1:2:nvargs
      switch varargin{k}
         case 'Ssnow'
            Ssnow = varargin{k+1};
         case 'Lf_t'
            Lf_t = varargin{k+1};
         case 'Lf_a'
            Lf_a = varargin{k+1};
         otherwise
            error(['I do not know what optional input argument ' varargin{k} ' means.']);     
      end
   end
end

n_phi = size(input,1);
% n_phi = size(input,1) - 1;

% Modified from Balco et al. (2008) skyline.m function to sort and linearly interpolate
% between inflection points

n_theta = [0:1:360];

if ~(isempty(input(1:end,1)) && isempty(input(1:end,2)));
    
	% step 2. Interpolate linearly between supplied horizon points.

	azR = (input(1:end,1));
	elR = (input(1:end,2));

    % sort in ascending order

	[azR j] = sort(azR);
	elR = elR(j);

	% pad for interpolation;

	azR2(2:(length(azR)+1)) = azR;
	elR2(2:length(elR)+1) = elR;
	azR2(1) = azR(length(azR)) - (360); 
	elR2(1) = elR(length(elR));
	azR2(length(azR)+2) = azR(1) + (360); 
	elR2(length(azR)+2) = elR(1);

	% interpolate;

	horiz = interp1(azR2,elR2,n_theta);

	% flip 

% 	horiz = horiz;
	
else

	horiz = zeros(length(n_theta),1);
	
end;
input2(:,1) = n_theta;
input2(:,2) = horiz;

n_phi = length(n_theta)-1;

values = combined(DD,dip_angle,n_phi,input2);

out.Stopog = values.numer*(m+2)/(2*pi);
out.Lf_calc = Lf_a*(m+2)/(m+1)*values.numer/values.denom;
out.Slambda = Lf_a/out.Lf_calc;
out.ST_calc = out.Slambda*out.Stopog*Ssnow;

% Unlike the version in CHLOE, Shield_num
% and shield_denom functions are combined
% into a single loop for efficiency

function cOut = combined(theta_n,phi_n,n_phi,in)
	if(n_phi == 0) 
	   n_phi=1;
	end

	rad = pi/180;
	phi_nrad = phi_n*rad;
	Cosphi_nrad = cos(phi_nrad);
	Sinphi_nrad = sin(phi_nrad);
	Tanphi_nrad = tan(phi_nrad);

	
	numer_ret = 0;
	denom_ret = 0;
	
	%Set values for this interval
	for k=1:n_phi
      sum_theta_numer = 0;
      sum_theta_denom = 0;
      phihorizon = in(k,2);

      %Calculate the horizon angle in the preferred coord. system
      phi_t = 90-phihorizon;

         %Calculate the surface projection angle, phi_s, for each increment of azimuth
         phi_s = 90+(atan(cos((k-theta_n) * rad) * Tanphi_nrad)/rad);
%          phi_s = 90+round(atan(cos((k-theta_n) * rad) * Tanphi_nrad)/rad);
%           shielding can't be greater than 90 degrees

         %Calculate whether topography or surface projection is the limiting lower angle
         if(phi_t < phi_s) 
            minphi = phi_t;
         else 
            minphi = phi_s;
         end
         %Do the numerical integration over phi, from vertical to limiting lower angle
         sum_phi_numer = 0;
         sum_phi_denom = 0;
         
			for phi=0.5:(minphi - 0.5)
                
                phirad = phi * rad;
                Cosphirad = cos(phirad);
                Sinphirad = sin(phirad);

                cosgamma = Cosphi_nrad * Cosphirad + Sinphi_nrad * Sinphirad * cos((theta_n - k)*rad);
                numer = (Cosphirad^m) * Sinphirad * cosgamma * rad;
                %Accumulate the pieces for integration over phi and multiply by delta phi
                sum_phi_numer = sum_phi_numer + numer;

                denom = (Cosphirad^m) * Sinphirad * rad;

                %Accumulate the pieces for integration over phi and multiply by delta phi
                sum_phi_denom = sum_phi_denom + denom;
			end
			
         %Accumulate the pieces for integration over theta. Multiplication by delta phi
         %will be done in the spreadsheet to minimize computation
         sum_theta_numer = sum_theta_numer + sum_phi_numer*rad;
         sum_theta_denom = sum_theta_denom + sum_phi_denom*rad;
		   
% 		end
		numer_ret = numer_ret + sum_theta_numer;
		denom_ret = denom_ret + sum_theta_denom;
	end

   cOut.numer = numer_ret;
   cOut.denom = denom_ret;

end %combined

end %topo

