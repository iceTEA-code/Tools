function out = get_tdsf_elev(sample,consts,elev_time,elev_vals)

% This function returns scaling factors as a function of time for a
% particular latitude, longitude, and atmospheric pressure.
%
% Syntax: out = get_tdsf(sample,consts);
% 
% Arguments are data structures. The `sample' structure must have the
% following fields:
%   sample.lat -- latitude, decimal degrees, S latitude is negative
%   sample.long -- longitude, decimal degrees. Attempts to deal with both
%   +180/-180 and 0/360 possibilities. 
%   sample.pressure -- atmospheric pressure in mb
%   sample.elevation -- elevation in meters.  
%	sample.scaling -- two letter scaling scheme ('SA','SF','ST','DE','DU','LI','LM') 
%
% Note that both pressure and elevation are now required for
% muons.  
%
% The `consts' structure is intended to be the al_be_consts data structure
% used in the Al-Be calculators, or a subset thereof. The following fields
% are used in this function:
%   consts.t_M 
%   consts.M
%   consts.lon_Rc
%   consts.lat_Rc
%   consts.t_Rc
%   consts.TTRc
%   consts.IHRc
%   consts.lat_pp_KCL
%   consts.lon_pp_KCL
%   consts.MM0_KCL
%   consts.SInf
%   consts.S
%
% See the m-file make_al_be_consts.m for more info on these. 
%
% The output argument has the following fields:
%   out.tv -- vector of time values in years
%
%   out.Rc_De, out.Rc_Du, out.Rc_Li, out.Rc_Lm -- cutoff rigidities at the
%   site (in GV) at the times in out.tv, corresponding to the four scaling
%   schemes
%
%   out.SF_De, out.SF_Du, out.SF_Li, out.SF_Lm -- scaling factors for
%   spallation at the site, at the times in out.tv, corresponding to the
%   four scaling schemes
%
%   out.SF_St -- scalar value -- scaling factor for spallation at the site
%   for the St scaling scheme.
%
%   out.ver -- version number
%
% This version corresponds to version 2.2 of the Al-Be age calculator. 
%
% Written by Greg Balco -- Berkeley Geochronology Center
% balcs@bgc.org
% December, 2008
% Part of the CRONUS-Earth online calculators: 
%      http://hess.ess.washington.edu/math
%
% Copyright 2001-2007, University of Washington
% 2007-2008, Berkeley Geochronology Center
% All rights reserved
% Developed in part with funding from the National Science Foundation.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License, version 2,
% as published by the Free Software Foundation (www.fsf.org).

%This version adapted from Nat Lifton's version of GeomagSX.m in Oct 2010.
%Changes include: deleted thickness correction factor, adjusted the inputs
%to only do one sample at a time. Got rid of the pressure/latitude
%corrections in this part of the code (we do those outside of the scaling
%factors)

% Modified by Shasta Marrero (NMT) in July 2011 to use the new
% nuclide-dependent scaling as implemented by Nat Lifton.  This is the
% Sato/Lifton scaling scheme.  Outputs from this particular scaling scheme
% include one set of scaling factors for each reaction of interest (10
% different scaling factors).  

% Modified by Nico Marrero in March 2014 to not calculate every scaling
% scheme as this is the largest bottleneck in the code - causing trapz to be
% called millions of times unnecessarily. Added sample.scaling to calculate
% specific scaling scheme. Defaults to old behavior if missing or set to 'all'

% Modified by Richard Selwyn Jones in February 2018 to calculate time-
% dependent production for all scaling schemes using an elevation vector.


% Change -180/180 longitude to 0-360
% The interpolation grids for the cutoff rigidities are set up on the basis
% of 0-360 longitude. 

if sample.long < 0
    sample.long = sample.long + 360;
end

% Verify a scaling scheme was requested and convert to lowercase
% otherwise, default to calculating all schemes
if (isfield(sample,'scaling'))
	sample.scaling = lower(sample.scaling);
else
	sample.scaling = 'all';
end

% Create a time vector. The spacing of this vector is variable: in the
% Holocene, it follows the spacing of the paleomagnetic data records, 
% mostly at 500 yrs. Between 12000 and 800,000 yr, the spacing is 1000 yr 
% to agree with the SINT-800 record. After 800,000 yt there is logarithmic
% spacing out to 10 million years. See the hard-copy documentation for 
% get_al_be_age.m for more details. 

% Age Relative to t0=2010
tv = [0:10:50 60:100:50060 51060:1000:2000060 logspace(log10(2001060),7,200)];

LiRc = zeros(1,length(tv));
DeRc = zeros(1,length(tv));
DuRc = zeros(1,length(tv));
LmRc = zeros(1,length(tv));
SaRc = zeros(1,length(tv));

% Need solar modulation for Lifton SF's
this_S = zeros(size(tv)) + consts.SInf;
this_S(1:120) = consts.S;
this_SPhi = zeros(size(tv)) + consts.SPhiInf; % Solar modulation potential for Sato et al. (2008)
this_SPhi(1:120) = consts.SPhi; % Solar modulation potential for Sato et al. (2008)


w = 0.2; % water content for Sato & Niita (2006)
% ethflux0 and thflux0 are the SLHL fluxes


% interpolate an M for tv > 7000...
% first sort the values to be interpolated for
[sorttv,indextv]=sort(tv(77:end));
% use "interpolate" to interpolate for the sorted values
temp_Munsorted = interpolate(consts.t_M,consts.M,sorttv);

% put the sorted values back in the correct order
temp_M = zeros(size(temp_Munsorted));
for v=1:length(sorttv)
    temp_M(indextv(v)) = temp_Munsorted(v);
end


    % Make up the Rc vectors.
    % Start with Lifton et al. 2005
    % First 15 from the data blocks

% Modified to work with new interpolation routines in MATLAB 2012a and later. 09/12
    [longi,lati,tvi] = meshgrid(sample.long,sample.lat,tv(1:76));
    LiRc(1:76) = interp3(consts.lon_Rc,consts.lat_Rc,consts.t_Rc,consts.TTRc,longi,lati,tvi);
    
% Fit to Trajectory-traced GAD dipole field as f(M/M0), as long-term average. This is the one I'm using now...
    dd = [6.89901,-103.241,522.061,-1152.15,1189.18,-448.004;];
    
% I replaced "end" with "length(tv)" because it was incorrectly defining
% the number of elements in the vector.  Not sure why.

    LiRc(77:(length(tv))) = temp_M.*(dd(1)*cos(d2r(sample.lat)) + ...
       dd(2)*(cos(d2r(sample.lat))).^2 + ...
       dd(3)*(cos(d2r(sample.lat))).^3 + ...
       dd(4)*(cos(d2r(sample.lat))).^4 + ...
       dd(5)*(cos(d2r(sample.lat))).^5 + ...
       dd(6)*(cos(d2r(sample.lat))).^6);
    %
    % Check for negative rigidity cutoffs and filter them out.
    %

    if (min(LiRc)<0)
        LiRc(LiRc<0)=0.0;
    end

    % Next, Desilets et al. 2006
    DeRc = LiRc; 
    SaRc = LiRc;
    
    % 1/21/10 Changed Dunai expression to include CALS3K.3
    DuRc(1:76) = interp3(consts.lon_Rc,consts.lat_Rc,consts.t_Rc,consts.IHRc,longi,lati,tvi);
    % The rest from Dunai 2001 equation 1
    DuRc(77:(length(tv))) = 14.9.*temp_M.*((cos(abs(d2r(sample.lat)))).^4);

    % Approximate paleo-pole-positions and field strengths from KC for < 7 ka
    LmLat = abs(90-angdist(sample.lat,sample.long,consts.lat_pp_KCL,consts.lon_pp_KCL));
    LmRc(1:76) = 14.9.*(consts.MM0_KCL).*((cos(abs(d2r(LmLat)))).^4);              
    LmRc(77:(length(tv))) = 14.9.*temp_M.*((cos(abs(d2r(sample.lat)))).^4);

    
    % MODIFICATION
    % Create the elevation vector
    elevation_tv = interp1(elev_time,elev_vals,tv,'linear','extrap'); % Interpolate elevations for tv times
    
    % Convert elevation to pressure
    if sample.lat < -60 % Antarctica-only
        pressure_tv = antatm(elevation_tv);
    else % Elsewhere; Computed from ERA-40 reanalysis data
        pressure_tv = ERA40atm(sample.lat,sample.long,elevation_tv);
    end
    
    % Do the age calculation. Interestingly, because all of the P(t)
    % functions are defined piecewise constant, it's not necessary to have
    % a zero-finding loop. We can just calculate the cumulative integral and
    % reverse-interpolate. 

    % Calculate the unweighted P(t) separately to be sent back in the results.
    % This is the surface production rate taking account of thickness. 
    % P_St is already calculated
	if (strcmpi(sample.scaling,'de') || strcmpi(sample.scaling,'all'))
		SF_De = desilets2006sp(pressure_tv,DeRc);
    end
    if (strcmpi(sample.scaling,'du') || strcmpi(sample.scaling,'all'))
		SF_Du = dunai2001sp(pressure_tv,DuRc);
    end
	if (strcmpi(sample.scaling,'li') || strcmpi(sample.scaling,'all'))
		SF_Li = lifton2006sp(pressure_tv,LiRc,this_S);
	end


% Nuclide-dependent scaling factors
% water content is an input, but it is not important because it gets
% cancelled out in the process of scaling
if (strcmpi(sample.scaling,'sa') || strcmpi(sample.scaling,'sf') || strcmpi(sample.scaling,'LSD') || strcmpi(sample.scaling,'LSDn') || strcmpi(sample.scaling,'all')) 
	scaling=LiftonSatoSX_mod(pressure_tv,SaRc,this_SPhi,0,consts); % MODIFICATION
	
	out.SF_Sf=scaling.sp;
	out.Rc_Sf=SaRc;

	out.SF_Sa10=scaling.Be;
	out.SF_Sa26=scaling.Al;
	out.SF_Sa3=scaling.He;
	out.SF_Sa14=scaling.C;
	out.SF_Sa36Ca=scaling.ClCa;
	out.SF_Sa36K=scaling.ClK;
	out.SF_Sa36Ti=scaling.ClTi;
	out.SF_Sa36Fe=scaling.ClFe;
	out.SF_Saeth=scaling.eth;
	out.SF_Sath=scaling.th;
	out.Rc_Sa=SaRc;
end
if (strcmpi(sample.scaling,'lm') || strcmpi(sample.scaling,'all'))
    SF_Lm = stone2000Rcsp_mod(pressure_tv,LmRc); % MODIFICATION
    out.SF_Lm=SF_Lm;
	out.Rc_Lm=LmRc;
end
if (strcmpi(sample.scaling,'st') || strcmpi(sample.scaling,'all'))
    SF_St = stone2000_mod(sample.lat,pressure_tv,1); % MODIFICATION
	out.SF_St=SF_St;
end
if (strcmpi(sample.scaling,'li') || strcmpi(sample.scaling,'all'))
	out.SF_Li=SF_Li;
	out.Rc_Li=LiRc;
end
if (strcmpi(sample.scaling,'du') || strcmpi(sample.scaling,'all'))
	out.SF_Du=SF_Du;
	out.Rc_Du=DuRc;
end
if (strcmpi(sample.scaling,'de') || strcmpi(sample.scaling,'all'))
	out.SF_De=SF_De;
	out.Rc_De=DeRc;
end

out.tv=tv;
out.elevation = elevation_tv;
out.pressure = pressure_tv;

end

