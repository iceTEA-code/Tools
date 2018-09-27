function out = lifton2006sp(h,Rc,S);

% Implements the Lifton, 2006 scaling scheme for spallation.
%
% Syntax: scalingfactor = lifton2006sp(h,Rc,S);
%
% Where:
%   h = atmospheric pressure (hPa)
%   Rc = cutoff rigidity (GV)
%   S solar modulation factor (nondimensional, see source paper)
%
% Vectorized. Send in scalars or vectors of common length. 
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

X = h.*1.019716;

% flatten low rigidities. The value of 1.907 comes from the source paper.

lowRc = find(Rc < 1.907);
Rc(lowRc) = 1.907 + zeros(size(lowRc));

% define constants

c = [1.8399 -1.1854e2 -4.9420e-2 8.0139e-1 1.2708e-4 9.4647e-1 -3.2208e-2 1.2688];
%delc = [1.0353e-2 2.6567 1.7512e-3 4.2170e-3 4.3896e-5 3.1630e-2 4.6392e-3 4.0327e-2];

t1 = c(1).*log(X.*S);
t2 = -S.*exp( (c(2).*S)./((Rc + 5.*S).^(2.*S)) );
t3 = c(3).*(X.^c(4));
t4 = c(5).*(( (Rc + 4.*S).*X).^c(6));
t5 = c(7).*((Rc + 4.*S).^c(8));

out = exp(t1 + t2 + t3 + t4 + t5);

