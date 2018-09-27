function out = get_tdsf(sample,consts);

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

% Version number

out.ver = '2.2-dev';

% Change -180/180 longitude to 0-360
% The interpolation grids for the cutoff rigidities are set up on the basis
% of 0-360 longitude. 

if sample.long < 0;
    sample.long = sample.long + 360;
end;

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

%tv = [0:500:6500 6900 7500:1000:11500 12000:1000:800000 logspace(log10(810000),7,200)];

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

%This needs to get added and checked/tested - Nat email Feb 2016
% > %Per Tatsuhiko Sato, personal communication, 2013, convert annually averaged Usoskin et al. (2011)
% > %solar modulation potential to Sato Force Field Potential due to different
% > %assumed Local Interstellar Spectrum and other factors
% >
% > SPhi = 1.1381076.*SPhi - 1.2738468e-4.*SPhi.^2;
% >
% > consts.SPhi = SPhi;
% > consts.SPhiInf = mean(SPhi);% Changed 12/13/11 to reflect updated SPhi values from Usoskin et al. (2011)

w = 0.2; % water content for Sato & Niita (2006)
% ethflux0 and thflux0 are the SLHL fluxes
%[nflux,ethflux0,thflux0] = NeutronsX(1013.25,0,this_SPhi(1),w);
%[pflux] = ProtonsX(1013.25,0,this_SPhi(1),consts);

%SXref = nflux + pflux;% Sato et al. (2008) Reference hadron flux integral >1 MeV
% SXref = P3n + P3p;% Sato et al. (2008) Reference 3He spallation production rate integral >1 MeV
% SXref = P10n + P10p;% Sato et al. (2008) Reference 10Be spallation production rate integral >1 MeV
% SXref = P14n + P14p;% Sato et al. (2008) Reference 14C spallation production rate integral >1 MeV
% SXref = P26n + P26p;% Sato et al. (2008) Reference 26Al spallation production rate integral >1 MeV


% interpolate an M for tv > 7000...
%first sort the values to be interpolated for
[sorttv,indextv]=sort(tv(77:end));
%use "interpolate" to interpolate for the sorted values
temp_Munsorted = interpolate(consts.t_M,consts.M,sorttv);

%put the sorted values back in the correct order
temp_M = zeros(size(temp_Munsorted));
for v=1:length(sorttv);
    temp_M(indextv(v))=temp_Munsorted(v);
end

%     mt(a) = t_simple .* 1.6; % mt is max time. 

    % Make up the Rc vectors.
    % Start with Lifton et al. 2005
    % First 15 from the data blocks

%   LiRc(1:76) = interp3(consts.lon_Rc,consts.lat_Rc,consts.t_Rc,consts.TTRc,sample.long,sample.lat,tv(1:76));

% Modified to work with new interpolation routines in MATLAB 2012a and later. 09/12
    [longi,lati,tvi] = meshgrid(sample.long,sample.lat,tv(1:76));
    LiRc(1:76) = interp3(consts.lon_Rc,consts.lat_Rc,consts.t_Rc,consts.TTRc,longi,lati,tvi);
    
%   New Equation - unpublished -  Fit to Trajectory-traced GAD dipole field as f(M/M0), as long-term average. This is the one I'm using now...

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
        %        error('LiRc contains negative rigidity cutoffs');
        LiRc(LiRc<0)=0.0;
    end

    % Next, Desilets et al. 2006
    DeRc = LiRc; 
    SaRc = LiRc;
    
%     1/21/10 Changed Dunai expression to include CALS3K.3
%  
    DuRc(1:76) = interp3(consts.lon_Rc,consts.lat_Rc,consts.t_Rc,consts.IHRc,longi,lati,tvi);
%     DuRc(a,1:70) = interp3(consts.lon_Rc,consts.lat_Rc,consts.t_DuRc,consts.IHRc,sample.long,sample.lat,tv(a,1:70));
    % The rest from Dunai 2001 equation 1
    DuRc(77:(length(tv))) = 14.9.*temp_M.*((cos(abs(d2r(sample.lat)))).^4);

    % Finally, paleomagnetically-corrected Lal
    % Same as Old Dunai > 7 ka
    %LmRc(a,:) = DuRc(a,:);
    % Approximate paleo-pole-positions and field strengths from KC for < 7 ka
    LmLat = abs(90-angdist(sample.lat,sample.long,consts.lat_pp_KCL,consts.lon_pp_KCL));
    LmRc(1:76) = 14.9.*(consts.MM0_KCL).*((cos(abs(d2r(LmLat)))).^4);              
    LmRc(77:(length(tv))) = 14.9.*temp_M.*((cos(abs(d2r(sample.lat)))).^4);

    % 4. Do the age calculation. Interestingly, because all of the P(t)
    % functions are defined piecewise constant, it's not necessary to have
    % a zero-finding loop. We can just calculate the cumulative integral and
    % reverse-interpolate. 

    % Calculate the unweighted P(t) separately to be sent back in the results.
    % This is the surface production rate taking account of thickness. 
    % P_St is already calculated
	if (strcmpi(sample.scaling,'de') || strcmpi(sample.scaling,'all'))
		SF_De = desilets2006sp(sample.pressure,DeRc);
    end
    if (strcmpi(sample.scaling,'du') || strcmpi(sample.scaling,'all'))
		SF_Du = dunai2001sp(sample.pressure,DuRc);
    end
	if (strcmpi(sample.scaling,'li') || strcmpi(sample.scaling,'all'))
		SF_Li = lifton2006sp(sample.pressure,LiRc,this_S);
	end
	
%calculate neutron fluxes for total, epithermal, and thermal at the site
%    [nflux,ethflux,thflux] = NeutronsX(sample.pressure,SaRc,this_SPhi,w);
%    [pflux] = ProtonsX(sample.pressure,SaRc,this_SPhi,consts);
%    SF_Sa = ((nflux + pflux))./SXref;
%     SF_Sa(a,:) = ((P3n + P3p).*sample.thickSF.*sample.shielding)./SXref;
%     SF_Sa(a,:) = ((P10n + P10p).*sample.thickSF.*sample.shielding)./SXref;
%     SF_Sa(a,:) = ((P14n + P14p).*sample.thickSF.*sample.shielding)./SXref;
%     SF_Sa(a,:) = ((P26n + P26p).*sample.thickSF.*sample.shielding)./SXref;

%scaling for thermal and epithermal neutrons
%SF_eth=ethflux/ethflux0;
%SF_th=thflux/thflux0;

%Nuclide-dependent scaling factors
% water content is an input, but it is not important because it gets
% cancelled out in the process of scaling
if (strcmpi(sample.scaling,'sa') || strcmpi(sample.scaling,'sf') || strcmpi(sample.scaling,'all')) 
	scaling=LiftonSatoSX(sample.pressure,SaRc,this_SPhi,0,consts);
	
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
	SF_Lm = stone2000Rcsp(sample.pressure,LmRc);
	out.SF_Lm=SF_Lm;
	out.Rc_Lm=LmRc;
end
if (strcmpi(sample.scaling,'st') || strcmpi(sample.scaling,'all'))
	SF_St = stone2000(sample.lat,sample.pressure,1); % scalar
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
	
    
    
    
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% 
% 
% 
% % Pull magnetic field strengths from SINT800 record
% 
% temp_M = interp1(consts.t_M,consts.M,tv(16:end));
% 
% 
% % Cutoff rigidity setup for Lifton 
% 
% LiRc = zeros(size(tv));
% % Geographic interpolation to get Holocene values
% LiRc(1:15) = interp3(consts.lon_Rc,consts.lat_Rc,consts.t_Rc,consts.TTRc,sample.long,sample.lat,tv(1:15));
% % The rest using Equation 6 of Lifton et al. 2005
% LiRc(16:end) = 15.765.*temp_M.*((cos(abs(d2r(sample.lat)))).^3.8);
% 
% % Cutoff rigidity setup for Desilets
% 
% DeRc = LiRc; % Holocene is the same as Lifton;
% % For older times, apply Desilets and Zreda 2003 Equation 19
% % Note that latitude has to be clipped at 55 for this to work. 
% DeLat = min([abs(sample.lat) 55]);
% ee = [-4.3077e-3 2.4352e-2 -4.6757e-3 3.3287e-4 -1.0993e-5 1.7037e-7 -1.0043e-9];
% ff = [1.4792e1 -6.6799e-2 3.5714e-3 2.8005e-5 -2.3902e-5 6.6179e-7 -5.0283e-9];
% DeRc(16:end) = (ee(1)+ff(1).*temp_M) + ...
%                 (ee(2)+ff(2).*temp_M).*(DeLat.^1) + ...
%                 (ee(3)+ff(3).*temp_M).*(DeLat.^2) + ...
%                 (ee(4)+ff(4).*temp_M).*(DeLat.^3) + ...
%                 (ee(5)+ff(5).*temp_M).*(DeLat.^4) + ...
%                 (ee(6)+ff(6).*temp_M).*(DeLat.^5) + ...
%                 (ee(7)+ff(7).*temp_M).*(DeLat.^6);
%   
% % Cutoff rigidity setup for Dunai
% 
% DuRc = zeros(size(tv));
% % Geographic interpolation to get Holocene values
% DuRc(1:15) = interp3(consts.lon_Rc,consts.lat_Rc,consts.t_Rc,consts.IHRc,sample.long,sample.lat,tv(1:15));
% % The rest from Dunai 2001 equation 1
% DuRc(16:end) = 14.9.*temp_M.*((cos(abs(d2r(sample.lat)))).^4);
% 
% % Cutoff rigidity setup for paleomagnetically-corrected Lal
% 
% % Same as Dunai > 7 ka
% LmRc = DuRc;
% % Approximate paleo-pole-positions and field strengths from KC for < 7 ka
% LmLat = abs(90-angdist(sample.lat,sample.long,consts.lat_pp_KCL,consts.lon_pp_KCL));
% LmRc(1:15) = 14.9.*(consts.MM0_KCL').*((cos(abs(d2r(LmLat)))).^4);              
%                 
% % Also need solar modulation for Lifton SF's
% this_S = zeros(size(LiRc)) + consts.SInf;
% this_S(1:19) = consts.S;
% 
% % Make four vectors of scaling factors
% % and one scalar (Stone 2000)
% % Note geographic scaling factor only -- no thickness or erosion
% % integration
% 
% out.SF_De = desilets2006sp(sample.pressure,DeRc);
% out.SF_Du = dunai2001sp(sample.pressure,DuRc);
% out.SF_Li = lifton2006sp(sample.pressure,LiRc,this_S);
% out.SF_Lm = stone2000Rcsp(sample.pressure,LmRc);
% out.SF_St = stone2000(sample.lat,sample.pressure,1); % scalar
% 
% % Also assign cutoff rigidities to out structure
% 
% out.Rc_De = DeRc;
% out.Rc_Du = DuRc;
% out.Rc_Li = LiRc;
% out.Rc_Lm = LmRc;
% out.Rc_Sa = LiRc;
% 
% % And report the time vector
% 
% out.tv = tv;
% 
% % That's it.
% 
% % Muon addition by Shasta Marrero June/July 2010
% 
% % Put the information from Sato's latitude scaling into a table. The table
% % was made using Sato's excel sheet to create the table based on total muon
% % flux (sum of positive and negative flux). The first row of the table is
% % the elevation in meters.  
% SatoTable=[...
% -500    0       500     1000    1500    2000    2500    3000    3500    4000    4500    5000    5500    6000;...
% 1       1       1       1       1       1       1       1       1       1       1       1       1       1;...
% 1       1       1       1       1       1       1       1       1       1       1       1       1       1;...
% 0.993    0.993	0.998	1.000	1.001	1.001	1.001	1.000	1.000	0.999	0.999	0.998	0.998	0.998;...
% 0.988    0.988	0.998	1.003	1.005	1.004	1.003	1.000	0.998	0.996	0.994	0.993	0.992	0.991;...
% 0.984    0.984	1.001	1.007	1.008	1.005	1.001	0.995	0.989	0.984	0.979	0.976	0.973	0.971;...
% 0.980    0.980	1.002	1.009	1.006	0.998	0.986	0.974	0.961	0.950	0.940	0.932	0.925	0.920;...
% 0.985    0.985	1.004	1.005	0.997	0.984	0.968	0.952	0.936	0.921	0.907	0.895	0.884	0.875;...
% 0.946    0.946	0.959	0.965	0.957	0.939	0.917	0.893	0.870	0.847	0.826	0.807	0.790	0.775;...
% 0.920    0.920	0.929	0.933	0.922	0.900	0.874	0.846	0.818	0.793	0.769	0.748	0.730	0.713;...
% 0.898    0.898	0.892	0.896	0.882	0.858	0.829	0.799	0.771	0.744	0.720	0.698	0.679	0.663;...
% 0.882    0.869	0.859	0.860	0.844	0.819	0.789	0.758	0.729	0.702	0.678	0.656	0.638	0.621;...
% 0.862    0.862	0.854	0.854	0.838	0.812	0.782	0.751	0.722	0.695	0.671	0.650	0.631	0.615;...
% 0.856    0.856	0.855	0.856	0.840	0.814	0.784	0.754	0.724	0.697	0.673	0.652	0.633	0.617;...
% 0.854    0.854	0.864	0.866	0.850	0.825	0.795	0.765	0.735	0.708	0.684	0.663	0.644	0.628;...
% 0.854    0.854   0.864   0.866   0.850   0.825   0.795   0.765   0.735   0.708   0.684   0.663   0.644   0.628];
% 
% 
% % Interpolate to the correct elevation.  The elevation should have been 
% % passed to the code in sample.elevation 
% %create a new table to store the new scaling factors for the sample
% %elevation and the rigidity cutoffs
% 
% newtable=zeros(16,2);
% % set up the rigidity cutoffs in the new table (first column)
% newtable(:,1)=[0,1,1.7,2.7,4.2,5.5,7.7,9.3,11.5,12.8,13.7,14.4,14.8,15.1,15.2,25]';
% 
% for i=1:15;
%         newtable(i,2)=interpolate(SatoTable(1,:), SatoTable((i+1),:),sample.elevation);
% end
% % Interpolate to all the rigidity cutoffs for each of the scaling schemes
% 
% %need to input the values for the times in sorted order, so sort them here
% [sortDe,indexDe]=sort(DeRc);
% [sortDu,indexDu]=sort(DuRc);
% [sortLi,indexLi]=sort(LiRc);
% [sortLm,indexLm]=sort(LmRc);
% 
% % Do the interpolation using the 'interpolate' function (written by Brian).
% % In this function, the cutoff rigidities must be sorted.
% 
% unsortedDe=interpolate(newtable(:,1), newtable(:,2), sortDe);
% unsortedDu=interpolate(newtable(:,1), newtable(:,2), sortDu);
% unsortedLi=interpolate(newtable(:,1), newtable(:,2), sortLi);
% unsortedLm=interpolate(newtable(:,1), newtable(:,2), sortLm);
% 
% % Need to resort the results based on the initial index for the vector
% for k=1:length(sortDe);
%     SFmu_De(indexDe(k))=unsortedDe(k);
%     SFmu_Du(indexDu(k))=unsortedDu(k);
%     SFmu_Li(indexLi(k))=unsortedLi(k);
%     SFmu_Lm(indexLm(k))=unsortedLm(k);
% end
% 
% 
% % Assign the out structure so these will be included in the output
% out.SFmu_De=SFmu_De;
% out.SFmu_Du=SFmu_Du;
% out.SFmu_Li=SFmu_Li;
% out.SFmu_Lm=SFmu_Lm;

% Yay! Muon scaling factors are done!

% % Adding code for Sato neutron/proton scaling
% tvSa = [0:100:50000 51000:1000:2000000 logspace(log10(2001000),7,200)];
% temp_MSa = interp1(sato_consts.t_M,sato_consts.M,tvSa(71:end));
% 
% %Calculate the cutoff rigidity to get Sato (Lifton's new values) 
% SaRc = zeros(length(tvSa),1);
% SaRc(1:70) = interp3(sato_consts.lon_Rc,sato_consts.lat_Rc,...
%     sato_consts.t_Rc,sato_consts.TTRc,sample.long,sample.lat,tvSa(1:70));
% ddSa = [6.89901,-103.241,522.061,-1152.15,1189.18,-448.004;];
% 
%     SaRc(71:end) = temp_MSa.*(ddSa(1)*cos(d2r(sample.lat)) + ...
%        ddSa(2)*(cos(d2r(sample.lat))).^2 + ...
%        ddSa(3)*(cos(d2r(sample.lat))).^3 + ...
%        ddSa(4)*(cos(d2r(sample.lat))).^4 + ...
%        ddSa(5)*(cos(d2r(sample.lat))).^5 + ...
%        ddSa(6)*(cos(d2r(sample.lat))).^6);
%    
% % Solar modulation potential for Sato et al. (2008)
% this_SPhi = zeros(size(tvSa)) + sato_consts.SPhiInf; 
% % Solar modulation potential for Sato et al. (2008)
% this_SPhi(1:114) = sato_consts.SPhi; 
% 
% %
% % Truncate the table of rigidity cutoffs and the time vector to go
% % back only 800Kyr to match the muon scaling.
% %
% tvSa=tvSa(1:1300);
% SaRc=SaRc(1:1300);
% this_SPhi=this_SPhi(1:1300);
% %Set up outputs
% 
% out.SF_Sa = sato2008(sample.pressure,SaRc,this_SPhi);
% out.Rc_Sa = SaRc;
% out.tvSa = tvSa;




