function out = LiftonSatoSX_mod(P,Rc,SPhi,w,consts)

% Implements the Lifton Sato et al scaling scheme for spallation.
%
% Syntax: scalingfactor = LiftonSatoSX_mod(P,Rc,SPhi,w,consts);
%
% Where:
%   P = atmospheric pressure (hPa)
%   Rc = cutoff rigidity (GV)
%   SPhi = solar modulation potential (Phi, see source paper)
%   w = fractional water content of ground (nondimensional)
%   
% Vectorized. Send in scalars or vectors of common length. 
% 
% Modified by Richard Jones in February 2018 to calculate production using 
% a time-dependent pressure vector.
%
% Modified by Shasta Marrero (NMT) to include the reactions for Ti & Fe to
% produce chlorine-36. July 2011.
%
% Written by Nat Lifton 2011, Purdue University
% Based on code by Greg Balco -- UW Cosmogenic Nuclide Lab
% balcs@u.washington.edu
% April, 2007
% Part of the CRONUS-Earth online calculators: 
%      http://hess.ess.washington.edu/math
%
% Copyright 2001-2007, University of Washington
% All rights reserved
% Developed in part with funding from the National Science Foundation.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License, version 2,
% as published by the Free Software Foundation (www.fsf.org).
%

% convert pressure to atmospheric depth
%X = P.*1.019716;

% Sref = 509.6; %1950-1850 mean
% Sref = 416.518; %Long-term mean
Sref = SPhi(1); %587.4 MV mean 2001-2010 from Usoskin et al 2011. USE THIS!
% Sref = 400; %Idealized solar minimum value, similar to long-term mean (11.4 ka)
% Changed from SPhi(1) 12/8/11 - more consistent with other calculations
% assuming idealized reference states such as SLHL pressure and Rc = 0. BUT
% does not normalize flux at t0 to 1 at SLHL - more like 0.92 - results in
% really high production predictions. Need to use the value from the 2001-2010.
% Go back to SPhi(1)... 12/15/11

Rcref = 0;
Href = 1013.25;
%Xref = Href.*1.019716;

% % Reduced versions don't do cross sections but includes thermal and
% % epithermal fluxes

% [nflux,ethflux,thflux] = NeutronsX(Href,Rcref,Sref,w);
% [pflux] = ProtonsX(Href,Rcref,Sref);
% SpRef = nflux + pflux;% Sato et al. (2008) Reference hadron flux integral >1 MeV
% EthRef = ethflux;
% ThRef = thflux;
% 
% % Reduced versions don't do cross sections but includes thermal and
% % epithermal fluxes
% 
% [nflux,ethflux,thflux] = NeutronsX(h,Rc,SPhi,w);
% [pflux] = ProtonsX(h,Rc,SPhi);
% Site.sp = ((nflux + pflux))./SpRef;
% Site.eth = ethflux./EthRef;
% Site.th = thflux./ThRef;


% Full version with cross sections and includes thermal and
% epithermal fluxes - reference
[nflux,P3n,P10n,P14n,P26n,P36Can,P36Kn,P36Tin,P36Fen] = NeutronsXS(Href,Rcref,Sref,w,consts);
[ethflux,thflux] = NeutronsLowE(Href,Rcref,Sref,w);
[pflux,P3p,P10p,P14p,P26p,P36Cap,P36Kp,P36Tip,P36Fep] = ProtonsXS(Href,Rcref,Sref,consts);
SpRef = nflux + pflux; % Sato et al. (2008) Reference hadron flux integral >1 MeV
HeRef = P3n + P3p;
BeRef = P10n + P10p;
CRef = P14n + P14p;
AlRef = P26n + P26p;
ClCaRef = P36Can + P36Cap;
ClKRef = P36Kn + P36Kp;
ClTiRef = P36Tin + P36Tip;
ClFeRef = P36Fen + P36Fep;
EthRef = ethflux;
ThRef = thflux;


% Full version with cross sections and includes thermal and
% epithermal fluxes - MODIFICATION
[nflux,P3n,P10n,P14n,P26n,P36Can,P36Kn,P36Tin,P36Fen] = NeutronsXS_mod(P,Rc,SPhi,w,consts);
[ethflux,thflux] = NeutronsLowE_mod(P,Rc,SPhi,w);
[pflux,P3p,P10p,P14p,P26p,P36Cap,P36Kp,P36Tip,P36Fep] = ProtonsXS_mod(P,Rc,SPhi,consts);


Site.sp = ((nflux + pflux))./SpRef;
Site.He = (P3n + P3p)./HeRef;
Site.Be = (P10n + P10p)./BeRef;
Site.C = (P14n + P14p)./CRef;
Site.Al = (P26n + P26p)./AlRef;
Site.ClCa = (P36Can + P36Cap)./ClCaRef;
Site.ClK = (P36Kn + P36Kp)./ClKRef;
Site.ClTi = (P36Tin + P36Tip)./ClTiRef;
Site.ClFe = (P36Fen + P36Fep)./ClFeRef;
Site.eth = ethflux./EthRef;
Site.th = thflux./ThRef;

out = Site;

end

