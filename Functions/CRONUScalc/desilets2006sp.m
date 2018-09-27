function out = desilets2006sp(h,Rc)

% This implements the Desilets et al. 2006 scaling scheme for spallation.
%
% Syntax: scalingfactor = desilets2006(h,Rc);
%
% where
%       h is atmospheric pressure (hPa)
%       Rc is cutoff rigidity (GV)
%
% Vectorized. Arguments should be scalars or vectors of equal length.
%
% Written by Greg Balco -- UW Cosmogenic Nuclide Lab
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
% convert pressure to atmospheric depth

x = h.*1.019716;

% enforce rigidity knee at 2 GV;

lowRc = find(Rc < 2);
Rc(lowRc) = 2 + zeros(size(lowRc));

% reference atmospheric depth at SL

sld = 1033;

% get altitude scaling factor

% get attenuation length

L = (sld - x)./(intOfBeta(sld,Rc) - intOfBeta(x,Rc));

% apply attenuation length

fofx = exp((sld-x)./L);

% get latitude scaling factor

alpha = 10.275;
k = 0.9615;

fofRc = 1 - exp(-alpha.*(Rc.^(-k)));

% total scaling factor 

out = fofRc.*fofx;


% -------- function definition for denominator of DZ2006 Eqns. 4 and 7 -------
% that is, the integral of the expression for Beta

function out = intOfBeta(x,Rc)

% coefficients from 2006 paper

n = 1.0177e-2;
alpha = 1.0207e-1;
k = -3.9527e-1;
a0 = 8.5236e-6;
a1 = -6.3670e-7;
a2 = -7.0814e-9;
a3 = -9.9182e-9;
a4 = 9.9250e-10;
a5 = 2.4925e-11;
a6 = 3.8615e-12;
a7 = -4.8194e-13;
a8 = -1.5371e-14;

% coefficients from 2003 paper

%n = 9.9741e-3;
%alpha = 4.5318e-1;
%k = -8.1613e-2;
%a0 = 6.3813e-6;
%a1 = -6.2639e-7;
%a2 = -5.1187e-9;
%a3 = -7.1914e-9;
%a4 = 1.1291e-9;
%a5 = 1.7400e-11;
%a6 = 2.5816e-12;
%a7 = -5.8588e-13;
%a8 = -1.2168e-14;


out = n.*x./(1 + exp(-alpha.*(Rc.^(-k)))) + (1/2).*(a0+a1.*Rc+a2.*Rc.*Rc).*(x.^2) + ...
    (1/3).*(a3+a4.*Rc+a5.*Rc.*Rc).*(x.^3) + (1/4).*(a6+a7.*Rc+a8.*Rc.*Rc).*(x.^4);

