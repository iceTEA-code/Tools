function out = stone2000Rcsp_mod(P,Rc)

% This implements a cutoff-rigidity rather than latitude based scaling
% scheme based on the Lal spallation polynomials. For use in
% paleomagnetically-corrected exposure-age calculations. 
%
% Syntax: scalingfactor = stone2000Rcsp_mod(P,Rc);
%
% Where
%   P = atmospheric pressure (hPa)
%   Rc = cutoff rigidity (GV)
%
% Vectorized in Rc and P.
%
% See the hard-copy documentation for details of the calculation.
%
% Written by Greg Balco -- UW Cosmogenic Nuclide Lab
% balcs@u.washington.edu
% March, 2007
% Part of the CRONUS-Earth online calculators: 
%      http://hess.ess.washington.edu/math
%
% Copyright 2001-2007, University of Washington
% All rights reserved
% Developed in part with funding from the National Science Foundation.
% 
% Modified by Richard Jones in February 2018 to calculate production using 
% a time-dependent pressure vector.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License, version 2,
% as published by the Free Software Foundation (www.fsf.org).

% 0. Input checks

% if length(P) > 1
%     error('stone2000Rcsp.m -- not vectorized in h');
% end

if max(Rc) > 21
    error('Rc greater than 21 in stone2000Rcsp.m');
end


% 1. Build the scaling factor = f(rigidity) function up to 14.9 GV

ilats_d = [0 10 20 30 40 50 60];
ilats_r = ilats_d.*pi./180; % convert to radians

% Convert latitude to rigidity using Elsasser formula (from Sandstrom)
% Rigidity Rc = Rc(0)cos^4(latitude)
% where Rc(0) = rigidity at equator, that is, 14.9 GV

iRcs = 14.9.*((cos(ilats_r)).^4);
for i = 1:length(P)
    iRcs_P(i,:) = iRcs;
end

% Spallogenic production at index rigidities;

% Constants from Table 1 of Stone(2000)

a = [31.8518 34.3699 40.3153 42.0983 56.7733 69.0720 71.8733];
b = [250.3193 258.4759 308.9894 512.6857 649.1343 832.4566 863.1927];
c = [-0.083393 -0.089807 -0.106248 -0.120551 -0.160859 -0.199252 -0.207069];
d = [7.4260e-5 7.9457e-5 9.4508e-5 1.1752e-4 1.5463e-4 1.9391e-4 2.0127e-4];
e = [-2.2397e-8 -2.3697e-8 -2.8234e-8 -3.8809e-8 -5.0330e-8 -6.3653e-8 -6.6043e-8];

% Apply Eqn. (2) of Stone (2000); MODIFICATION
sf = zeros(numel(P'),numel(a));
sf = a + (b .* exp(P'./(-150))) + (c.*P') + (d.*(P'.^2)) + (e.*(P'.^3));

% Extend to zero rigidity - scaling factor does not change from that at 60
% degrees --

iRcs_P(:,8) = 0; sf(:,8) = sf(:,7); % MODIFICATION

% Extend to 21 GV by fitting a log-log line to the latitude 0-20 values,
% i.e. where rigidity is greater than 10 GV. According to Quemby and Wenk, 
% as summarized  in Sandstrom, log(rigidity) vs. log(nucleon intensity) 
% ought to be linear above 10 GV. Note that this is speculative, but 
% relatively unimportant, as the approximation is pretty much only used 
% for low latides for a short time in the Holocene. 

% MODIFICATION
fits = zeros(length(iRcs_P(:,1)),2);
add_sf = zeros(length(iRcs_P(:,1)),6);
add_iRcs = zeros(length(iRcs_P(:,1)),6);
for ii = 1:length(iRcs_P(:,1))
    fits(ii,:) = polyfit(log(iRcs_P(ii,1:3)),log(sf(ii,1:3)),1 );
    add_sf(ii,:) = exp( log(sf(ii,1)) + fits(ii,1).*( log(21:-1:16) - log(iRcs_P(ii,1)) ) ) ;
    add_iRcs(ii,:) = [21 20 19 18 17 16];
end
iRcs_P = [add_iRcs iRcs_P];
sf = [add_sf sf];
sortsf = zeros(size(sf));
sortiRcs = zeros(size(iRcs_P));
%indexiRcs = zeros(1,length(iRcs_P(1,:)));
for iii = 1:length(sf(:,1))
    [sortiRcs(iii,:),indexiRcs] = sort(iRcs_P(iii,:));
    for m = 1:length(sf(1,:))
        sortsf(iii,indexiRcs(m)) = sf(iii,m);
    end
end


% Interpolate, return
[sortLm,indexLm] = sort(Rc);

unsortedLm_matrix = zeros(length(sortsf),length(Rc));
resortedLm_matrix = unsortedLm_matrix;
for iiii = 1:length(sortsf)
    unsortedLm_matrix(iiii,:) = interpolate(sortiRcs(iiii,:),sortsf(iiii,:),sortLm);
    for k = 1:length(Rc)
        resortedLm_matrix(iiii,indexLm(k)) = unsortedLm_matrix(iiii,k);
    end
end
SF_Lm = diag(resortedLm_matrix);

out=SF_Lm;
%out = interp1(iRcs,sf,Rc);

end
