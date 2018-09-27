function out = dunai2001sp(h,Rc);

% This function implements the Dunai, 2001 scaling scheme for spallation.
%
% Syntax: out = dunai2001sp(h,Rc);
%
% where     h is atmospheric pressure (hPa)
%           Rc is cutoff rigidity (GV)
%
% Vectorized. Use scalar or vector-of-equal-length arguments.
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
% convert pressure to atmospheric depth between sea level and site

delz = (1013.25-h).*1.019716;

% constants

A = 0.5221;
B = -1.7211;
C = 0.3345;
X = 4.2822;
Y = 0.4952;

a = 17.183;
b = 2.060;
c = 5.9164;
x = 2.2964;
y = 130.11;

% sea level scaling factor

N1030 = Y + A./( (1 + exp(-(Rc-X)./B)).^C);

% attenuation length

L = y + a./( (1 + exp(-(Rc-x)./b)).^c);

% total scaling factor

out = N1030.*exp(delz./L);


