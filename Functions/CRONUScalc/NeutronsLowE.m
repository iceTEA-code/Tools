function [ethflux,thflux] = NeutronsLowE(h,Rc,s,w)
%reduced version that does not do cross-sections but includes thermal and
%epithermal spectra
% Sato et al. (2008) Neutron Spectrum
% Analytical Function Approximation (PARMA)
% Implemented in MATLAB by Nat Lifton, May 2009
% University of Arizona, lifton@email.arizona.edu

% Copyright 2011, Purdue University and University of Arizona
% All rights reserved
% Developed in part with funding from the National Science Foundation.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License, version 2,
% as published by the Free Software Foundation (www.fsf.org).

x = h.*1.019716; % Convert pressure (hPa) to atm depth (g/cm2)

% E = logspace(-8,5,1000);
E = logspace(-8,0,200);
% E = [1.1295 11.295 112.95 1129.5 11295];

% Flatten low rigidities.

lowRc = find(Rc < 1.0);
Rc(lowRc) = 1.0 + zeros(size(lowRc));

%nflux = zeros(length(Rc));

%w = 0.2; % water content, from 0-1
%s = 1700;
%Rc = 12;
%x = 1030;
Et = 2.5e-8; % Thermal Neutron Energy in MeV

% Integrated neutron flux <15 MeV

smin = 400; %units of MV
smax = 1200; %units of MV

a6 = 1.8882e-4;
a7 = 4.4791e-1;
a8 = 1.4361e-3;
a12 = 1.4109e-2;

b11min = 2.5702e1;
b11max = -6.9221;
b12min = -5.0931e-1;
b12max = 1.1336;
b13min= 7.4650;
b13max = 2.6961e1;
b14min = 1.2313e1;
b14max = 1.1746e1;
b15min = 1.0498;
b15max = 2.6171;
b21min = 6.5143e-3;
b21max = 5.3211e-3;
b22min = 3.3511e-5;
b22max = 8.4899e-5;
b23min = 9.4415e-4;
b23max = 2.0704e-3;
b24min = 1.2088e1;
b24max = 1.1714e1;
b25min = 2.7782;
b25max = 3.8051;
b31min = 9.8364e-1;
b31max = 9.7536e-1;
b32min = 1.4964e-4;
b32max = 6.4733e-4;
b33min = -7.3249e-1;
b33max = -2.2750e-1;
b34min = -1.4381;
b34max = 2.0713;
b35min = 2.7448;
b35max = 2.1689;
b41min = 8.8681e-3;
b41max = 9.1435e-3;
b42min = -4.3322e-5;
b42max = -6.4855e-5;
b43min = 1.7293e-2;
b43max = 5.8179e-3;
b44min = -1.0836;
b44max = 1.0168;
b45min = 2.6602;
b45max = 2.4504;

b121 = 9.31e-1;
b122 = 3.70e-2;
b123 = -2.02;
b124 = 2.12;
b125 = 5.34;
b131 = 6.67e-4;
b132 = -1.19e-5;
b133 = 1.00e-4;
b134 = 1.45;
b135 = 4.29;

% Basic Spectrum

b51 = 9.7337e-4;
b52 = -9.6635e-5;
b53 = 1.2121e-2;
b54 = 7.1726;
b55 = 1.4601;
b91 = 5.7199e2;
b92 = 7.1293;
b93 = -1.0703e2;
b94 = 1.8538;
b95 = 1.2142;
b101 = 6.8552e-4;
b102 = 2.7136e-5;
b103 = 5.7823e-4;
b104 = 8.8534;
b105 = 3.6417;
b111 = -5.0800e-1;
b112 = 1.4783e-1;
b113 = 1.0068;
b114 = 9.1556;
b115 = 1.6369;

c1 = 2.3555e-1; % lethargy^-1
c2 = 2.3779; % MeV
c3 = 7.2597e-1;
c5 = 1.2391e2; % MeV
c6 = 2.2318; % MeV
c7 = 1.0791e-3; % lethargy^-1
c8 = 3.6435e-12; % MeV
c9 = 1.6595;
c10 = 8.4782e-8; % MeV
c11 = 1.5054;

% Ground-Level Spectrum

h31 = -2.5184e1;
h32 = 2.7298;
h33 = 7.1526e-2;
h51 = 3.4790e-1;
h52 = 3.3493;
h53 = -1.5744;

g1 = -0.023499;
g2 = -0.012938;
g3 = 10.^(h31 + h32./(w + h33));
g4 = 9.6889e-1;
g5 = h51 + h52.*w + h53.*(w.^2);

fG = 10.^(g1 + g2.*log10(E./g3).*(1-tanh(g4.*log10(E./g5))));

% Thermal Neutron Spectrum

h61 = 1.1800e-1;
h62 = 1.4438e-1;
h63 = 3.8733;
h64 = 6.5298e-1;
h65 = 4.2752e1;

g6 = (h61 + h62.*exp(-h63.*w))./(1 + h64.*exp(-h65.*w));

PhiT = g6.*((E./Et).^2).*exp(-E./Et);

% Total Ground-Level Flux

PhiB = zeros(1,length(E));
PhiG = zeros(1,length(E));
PhiGMev = zeros(1,length(E));
ethflux = zeros(1,length(Rc));
thflux = zeros(1,length(Rc));

for a = 1:length(Rc)
    
    a1min = b11min + b12min.*Rc(a) + b13min./(1 + exp((Rc(a) - b14min)./b15min));
    a1max = b11max + b12max.*Rc(a) + b13max./(1 + exp((Rc(a) - b14max)./b15max));
    a2min = b21min + b22min.*Rc(a) + b23min./(1 + exp((Rc(a) - b24min)./b25min));
    a2max = b21max + b22max.*Rc(a) + b23max./(1 + exp((Rc(a) - b24max)./b25max));
    a3min = b31min + b32min.*Rc(a) + b33min./(1 + exp((Rc(a) - b34min)./b35min));
    a3max = b31max + b32max.*Rc(a) + b33max./(1 + exp((Rc(a) - b34max)./b35max));
    a4min = b41min + b42min.*Rc(a) + b43min./(1 + exp((Rc(a) - b44min)./b45min));
    a4max = b41max + b42max.*Rc(a) + b43max./(1 + exp((Rc(a) - b44max)./b45max));
    
    a5 = b51 + b52.*Rc(a) + b53./(1 + exp((Rc(a) - b54)./b55));
    a9 = b91 + b92.*Rc(a) + b93./(1 + exp((Rc(a) - b94)./b95));
    a10 = b101 + b102.*Rc(a) + b103./(1 + exp((Rc(a) - b104)./b105));
    a11 = b111 + b112.*Rc(a) + b113./(1 + exp((Rc(a) - b114)./b115));
    
    b5 = b121 + b122.*Rc(a) + b123./(1 + exp((Rc(a) - b124)./b125));
    b6 = b131 + b132.*Rc(a) + b133./(1 + exp((Rc(a) - b134)./b135));

    c4 = a5 + a6.*x./(1 + a7.*exp(a8.*x)); % lethargy^-1
    c12 = a9.*(exp(-a10.*x) + a11.*exp(-a12.*x)); % MeV

    PhiLmin = a1min.*(exp(-a2min.*x) - a3min.*exp(-a4min.*x)); %length of Rc
    PhiLmax = a1max.*(exp(-a2max.*x) - a3max.*exp(-a4max.*x)); %length of Rc
    
    f3 = b5 + b6.*x;
    f2 = (PhiLmin - PhiLmax)./(smin.^f3 - smax.^f3);
    f1 = PhiLmin - f2.*smin.^f3;

    PhiL = f1 + f2.*s(a).^f3;

    PhiB = (c1.*(E./c2).^c3).*exp(-E./c2) + c4.*exp((-(log10(E) - log10(c5)).^2)./(2.*(log10(c6)).^2))...
    + c7.*log10(E./c8).*(1 + tanh(c9.*log10(E./c10))).*(1 - tanh(c11.*log10(E./c12)));

    PhiG = PhiL.*(PhiB.*fG + PhiT);
    PhiGMev = PhiG./E;
    
    clipindexethmax = find(E <= 1e-1, 1, 'last' ); %Find epithermal neutron E
    clipindexethmmin = find(E >= 5e-7, 1, 'first' ); %Find epithermal neutron E    
    ethflux(a) = trapz(E(clipindexethmmin:clipindexethmax),PhiGMev(clipindexethmmin:clipindexethmax));
    clipindexth = find(E <= 5e-7, 1, 'last' ); %Find thermal neutron E
    thflux(a) = trapz(E(1:clipindexth),PhiGMev(1:clipindexth));

end    

% Plot it

% figure;clf;
% semilogx(E,PhiG);
% 
% figure;clf;
% loglog(E,PhiGMev); hold on;
% ylim([0 .1]);
