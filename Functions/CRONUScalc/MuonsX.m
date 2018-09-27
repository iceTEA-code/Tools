function mflux = MuonsX(h,Rc,s)

% Sato et al. (2008) Muon Spectra
% Analytical Function Approximation (PARMA)
% Implemented in MATLAB by Nat Lifton, May 2009
% University of Arizona, lifton@email.arizona.edu

% Copyright 2009, University of Arizona
% All rights reserved
% Developed in part with funding from the National Science Foundation.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License, version 2,
% as published by the Free Software Foundation (www.fsf.org).

x = h.*1.019716; % Convert pressure (hPa) to atm depth (g/cm2)

E = logspace(1,5.9030,200); %(in MeV)
% E = logspace(1,5.3010,200); %(in MeV)

Emu = 105.658; % in Mev, Rest energy of muon
alpha3 = 3.7; % Muon spectrum power law exponent
Beta = sqrt(1-(Emu./(Emu + E)).^2); % Particle speed relative to light
c = 3.0e8; % speed of light, in m/s
p = sqrt((E.^2 + 2.*E.*Emu)); % in MeV/c

% Flatten low rigidities

lowRc = find(Rc < 1.0);
Rc(lowRc) = 1.0 + zeros(size(lowRc));

smin = 400; %units of MV
smax = 1200; %units of MV

Phimu = zeros(length(Rc),length(E));


% Negative muon coefficients

u1n = 5.8214e9;
u2n = 3.6228e-3;
u3n = 1.0240;
u4n = 4.5141e-3;
u5n = 3.1992e8;

Phimn = u1n.*(exp(-u2n.*x) - u3n.*exp(-u4n.*x)) + u5n;

w111nmin = 2.0899e3;
w112nmin = 1.2110e2;
w113nmin = -9.2925e2;
w114nmin = 6.8558;
w115nmin = 3.2929;
w111nmax = 2.4185e3;
w112nmax = 1.1240e2;
w113nmax = -8.9497e2;
w114nmax = 7.4497;
w115nmax = 3.5522;
w121nmin = -5.6641;
w122nmin = -6.4998e-1;
w123nmin = 3.5830;
w124nmin = 8.8799e-1;
w125nmin = 3.7337;
w121nmax = -5.6115;
w122nmax = -6.5095e-1;
w123nmax = 3.3115;
w124nmax = 7.7616e-1;
w125nmax = 3.7607;
w131nmin = 1.1807e-2;
w132nmin = 1.5847e-3;
w133nmin = -1.2543e-2;
w134nmin = 3.4411;
w135nmin = 3.6455;
w131nmax = 1.1804e-2;
w132nmax = 1.5798e-3;
w133nmax = -1.2480e-2;
w134nmax = 3.4818;
w135nmax = 3.5926;
w141nmin = -2.5853e-6;
w142nmin = -7.9871e-7;
w143nmin = 2.5370e-5;
w144nmin = 4.9450;
w145nmin = 3.7213;
w141nmax = -2.5196e-6;
w142nmax = -7.9341e-7;
w143nmax = 2.5343e-5;
w144nmax = 4.9219;
w145nmax = 3.7354;
w151nmin = 1.8671e-9;
w152nmin = -1.9787e-10;
w153nmin = -1.7061e-8;
w154nmin = 5.1157;
w155nmin = 4.2354;
w151nmax = 1.8602e-9;
w152nmax = -2.0122e-10;
w153nmax = -1.7016e-8;
w154nmax = 5.1424;
w155nmax = 4.2718;

w211nmin = 8.5946e1;
w212nmin = -5.8637;
w213nmin = 3.6872e2;
w214nmin = 4.8178;
w215nmin = 3.2984;
w211nmax = 8.6974e1;
w212nmax = -5.8773;
w213nmax = 3.7230e2;
w214nmax = 4.6802;
w215nmax = 3.2996;
w221nmin = 3.4175;
w222nmin = 7.9022e-2;
w223nmin = -5.2936e-1;
w224nmin = 6.8789;
w225nmin = 1.0647;
w221nmax = 3.4184;
w222nmax = 7.8730e-2;
w223nmax = -5.3162e-1;
w224nmax = 6.8578;
w225nmax = 1.0891;
w231nmin = -3.3253e-3;
w232nmin = -1.4941e-4;
w233nmin = 1.8630e-3;
w234nmin = 7.0358;
w235nmin = 6.0158e-1;
w231nmax = -3.3203e-3;
w232nmax = -1.4962e-4;
w233nmax = 1.8556e-3;
w234nmax = 7.0391;
w235nmax = 6.0068e-1;
w241nmin = -2.6862e-6;
w242nmin = -8.9985e-8;
w243nmin = -2.7068e-6;
w244nmin = 7.0511;
w245nmin = 4.6369e-1;  
w241nmax = -2.6832e-6;
w242nmax = -8.9349e-8;
w243nmax = -2.7056e-6;
w244nmax = 7.0489;
w245nmax = 4.6511e-1;
w251nmin = 2.3372e-9;
w252nmin = 1.5003e-10;
w253nmin = 1.1941e-9;
w254nmin = 7.0490;
w255nmin = 3.5646e-1;
w251nmax = 2.3300e-9;
w252nmax = 1.4973e-10;
w253nmax = 1.1994e-9;
w254nmax = 7.0449;
w255nmax = 3.6172e-1;

w311nmin = 7.8736e-1;
w312nmin = -1.8004e-2;
w313nmin = -3.0414e-1;
w314nmin = 1.4479e1;
w315nmin = 5.6128;
w311nmax = 8.1367e-1;
w312nmax = -2.4784e-2;
w313nmax = -3.1104e-1;
w314nmax = 1.0553e1;
w315nmax = 3.6057;
w321nmin = 2.1362e-3;
w322nmin = 4.9866e-5;
w323nmin = 1.4331e-3;
w324nmin = 8.1043;
w325nmin = 3.4619;
w321nmax = 6.6470e-4;
w322nmax = 1.3546e-4;
w323nmax = 1.8371e-3;
w324nmax = 9.2913;
w325nmax = 2.3906;
w331nmin = -6.0480e-6;
w332nmin = -1.3554e-7;
w333nmin = -3.9433e-6;
w334nmin = 7.8291;
w335nmin = 4.3398;
w331nmax = -3.7978e-6;
w332nmax = -2.9193e-7;
w333nmax = -2.5834e-6;
w334nmax = 9.6668;
w335nmax = 1.3763;
w341nmin = 6.6770e-9;
w342nmin = 1.0885e-12;
w343nmin = 1.5756e-9;
w344nmin = 2.2697e1;
w345nmin = 1.9922;
w341nmax = 2.7492e-9;
w342nmax = 3.3458e-10;
w343nmax = 2.3109e-9;
w344nmax = 1.0281e1;
w345nmax = 1.3660;
w351nmin = -3.0952e-12;
w352nmin = 3.8044e-14;
w353nmin = 7.4580e-13;
w354nmin = 7.8473;
w355nmin = 2.0013;
w351nmax = -1.8076e-12;
w352nmax = -4.1711e-14;
w353nmax = 4.6284e-13;
w354nmax = 4.5439;
w355nmax = 4.7886e-1;

h51n = 5.6500e-1;
h52n = 1.2100e-2;
h53n = -3.5700e-1;
h54n = 4.7300;
h55n = 1.4600;
h61n = 8.8000e-5;
h62n = -3.8900e-6;
h63n = 4.9100e-4;
h64n = 4.5100;
h65n = 1.7200;

% Positive muon coefficients
u1p = 6.2603e9;
u2p = 3.4320e-3;
u3p = 1.0131;
u4p = 4.1817e-3;
u5p = 3.7543e8;

Phimp = u1p.*(exp(-u2p.*x) - u3p.*exp(-u4p.*x)) + u5p;

w111pmin = 2.0538e3;
w112pmin = 1.2598e2;
w113pmin = -1.0131e3;
w114pmin = 6.1791;
w115pmin = 3.4718;
w111pmax = 2.3945e3;
w112pmax = 1.1790e2;
w113pmax = -9.4920e2;
w114pmax = 7.0369;
w115pmax = 3.8446;
w121pmin = -5.6688;
w122pmin = -6.5475e-1;
w123pmin = 3.5933;
w124pmin = 1.3137;
w125pmin = 3.2223;
w121pmax = -5.6246;
w122pmax = -6.5784e-1;
w123pmax = 3.2754;
w124pmax = 1.0604;
w125pmax = 3.3353;
w131pmin = 1.1700e-2;
w132pmin = 1.5748e-3;
w133pmin = -1.2521e-2;
w134pmin = 3.2601;
w135pmin = 3.6451;
w131pmax = 1.1736e-2;
w132pmax = 1.5714e-3;
w133pmax = -1.2383e-2;
w134pmax = 3.3054;
w135pmax = 3.5833;
w141pmin = -2.3130e-6;
w142pmin = -7.5964e-7;
w143pmin = 2.4832e-5;
w144pmin = 4.9409;
w145pmin = 3.7979;
w141pmax = -2.2412e-6;
w142pmax = -7.5644e-7;
w143pmax = 2.4834e-5;
w144pmax = 4.8875;
w145pmax = 3.8034;
w151pmin = 1.7430e-9;
w152pmin = -2.2205e-10;
w153pmin = -1.6916e-8;
w154pmin = 5.1206;
w155pmin = 4.3875;
w151pmax = 1.7462e-9;
w152pmax = -2.2603e-10;
w153pmax = -1.6852e-8;
w154pmax = 5.1768;
w155pmax = 4.3997;

w211pmin = 8.4834e1;
w212pmin = -5.7723;
w213pmin = 3.7035e2;
w214pmin = 4.8084;
w215pmin = 3.3589;
w211pmax = 8.7301e1;
w212pmax = -5.9021;
w213pmax = 3.7664e2;
w214pmax = 4.5920;
w215pmax = 3.3933;
w221pmin = 3.4086;
w222pmin = 7.8728e-2;
w223pmin = -5.2000e-1;
w224pmin = 6.8730;
w225pmin = 1.0869;
w221pmax = 3.4070;
w222pmax = 7.8501e-2;
w223pmax = -5.2268e-1;
w224pmax = 6.8422;
w225pmax = 1.0916;
w231pmin = -3.3162e-3;
w232pmin = -1.4917e-4;
w233pmin = 1.8524e-3;
w234pmin = 7.0237;
w235pmin = 6.0692e-1;
w231pmax = -3.3141e-3;
w232pmax = -1.4904e-4;
w233pmax = 1.8518e-3;
w234pmax = 7.0237;
w235pmax = 6.1137e-1;
w241pmin = -2.6781e-6;
w242pmin = -8.8820e-8;
w243pmin = -2.7098e-6;
w244pmin = 7.0420;
w245pmin = 4.6845e-1;
w241pmax = -2.6774e-6;
w242pmax = -8.8086e-8;
w243pmax = -2.7055e-6;
w244pmax = 7.0422;
w245pmax = 4.7162e-1;
w251pmin = 2.3267e-9;
w252pmin = 1.4896e-10;
w253pmin = 1.2010e-9;
w254pmin = 7.0431;
w255pmin = 3.6378e-1;
w251pmax = 2.3187e-9;
w252pmax = 1.4872e-10;
w253pmax = 1.2045e-9;
w254pmax = 7.0488;
w255pmax = 3.6659e-1;

w311pmin = 7.6040e-1;
w312pmin = -1.8020e-2;
w313pmin = -2.7253e-1;
w314pmin = 1.1292e1;
w315pmin = 5.3901;
w311pmax = 9.2327e-1;
w312pmax = -2.9590e-2;
w313pmax = -4.2838e-1;
w314pmax = 9.6573;
w315pmax = 4.0023;
w321pmin = 2.0613e-3;
w322pmin = 6.1719e-5;
w323pmin = 1.7751e-3;
w324pmin = 7.5508;
w325pmin = 3.9262;
w321pmax = 8.4438e-4;
w322pmax = 1.3392e-4;
w323pmax = 1.8096e-3;
w324pmax = 9.2554;
w325pmax = 2.4406;
w331pmin = -5.9644e-6;
w332pmin = -1.4795e-7;
w333pmin = -4.1301e-6;
w334pmin = 7.5298;
w335pmin = 4.3879;
w331pmax = -3.9078e-6;
w332pmax = -2.8780e-7;
w333pmax = -2.4920e-6;
w334pmax = 9.7445;
w335pmax = 1.4865;
w341pmin = 6.4640e-9;
w342pmin = -9.2764e-12;
w343pmin = 1.7352e-9;
w344pmin = 2.3633e1;
w345pmin = 1.6729;
w341pmax = 1.9852e-9;
w342pmax = 3.5716e-10;
w343pmax = 2.9465e-9;
w344pmax = 1.0431e1;
w345pmax = 1.9364;
w351pmin = -3.2101e-12;
w352pmin = 5.4637e-14;
w353pmin = 9.2092e-13;
w354pmin = 7.5423;
w355pmin = 2.6570;
w351pmax = -1.7751e-12;
w352pmax = -3.1711e-14;
w353pmax = 4.7927e-13;
w354pmax = 4.2050;
w355pmax = 7.4704e-1;

h51p = 5.0600e-1;
h52p = 1.3000e-2;
h53p = -3.9400e-1;
h54p = 4.1200;
h55p = 1.3300;
h61p = 1.3900e-4;
h62p = 6.9500e-6;
h63p = 7.4700e-4;
h64p = 3.7200;
h65p = 1.9700;

for a = 1:length(Rc)
    %Negative Muons
    v11nmin = w111nmin + w112nmin.*Rc(a) + w113nmin./(1 + exp((Rc(a) - w114nmin)./w115nmin));
    v11nmax = w111nmax + w112nmax.*Rc(a) + w113nmax./(1 + exp((Rc(a) - w114nmax)./w115nmax));
    v12nmin = w121nmin + w122nmin.*Rc(a) + w123nmin./(1 + exp((Rc(a) - w124nmin)./w125nmin));
    v12nmax = w121nmax + w122nmax.*Rc(a) + w123nmax./(1 + exp((Rc(a) - w124nmax)./w125nmax));
    v13nmin = w131nmin + w132nmin.*Rc(a) + w133nmin./(1 + exp((Rc(a) - w134nmin)./w135nmin));
    v13nmax = w131nmax + w132nmax.*Rc(a) + w133nmax./(1 + exp((Rc(a) - w134nmax)./w135nmax));
    v14nmin = w141nmin + w142nmin.*Rc(a) + w143nmin./(1 + exp((Rc(a) - w144nmin)./w145nmin));
    v14nmax = w141nmax + w142nmax.*Rc(a) + w143nmax./(1 + exp((Rc(a) - w144nmax)./w145nmax));
    v15nmin = w151nmin + w152nmin.*Rc(a) + w153nmin./(1 + exp((Rc(a) - w154nmin)./w155nmin));
    v15nmax = w151nmax + w152nmax.*Rc(a) + w153nmax./(1 + exp((Rc(a) - w154nmax)./w155nmax));
    v21nmin = w211nmin + w212nmin.*Rc(a) + w213nmin./(1 + exp((Rc(a) - w214nmin)./w215nmin));
    v21nmax = w211nmax + w212nmax.*Rc(a) + w213nmax./(1 + exp((Rc(a) - w214nmax)./w215nmax));
    v22nmin = w221nmin + w222nmin.*Rc(a) + w223nmin./(1 + exp((Rc(a) - w224nmin)./w225nmin));
    v22nmax = w221nmax + w222nmax.*Rc(a) + w223nmax./(1 + exp((Rc(a) - w224nmax)./w225nmax));
    v23nmin = w231nmin + w232nmin.*Rc(a) + w233nmin./(1 + exp((Rc(a) - w234nmin)./w235nmin));
    v23nmax = w231nmax + w232nmax.*Rc(a) + w233nmax./(1 + exp((Rc(a) - w234nmax)./w235nmax));
    v24nmin = w241nmin + w242nmin.*Rc(a) + w243nmin./(1 + exp((Rc(a) - w244nmin)./w245nmin));
    v24nmax = w241nmax + w242nmax.*Rc(a) + w243nmax./(1 + exp((Rc(a) - w244nmax)./w245nmax));
    v25nmin = w251nmin + w252nmin.*Rc(a) + w253nmin./(1 + exp((Rc(a) - w254nmin)./w255nmin));
    v25nmax = w251nmax + w252nmax.*Rc(a) + w253nmax./(1 + exp((Rc(a) - w254nmax)./w255nmax));
    v31nmin = w311nmin + w312nmin.*Rc(a) + w313nmin./(1 + exp((Rc(a) - w314nmin)./w315nmin));
    v31nmax = w311nmax + w312nmax.*Rc(a) + w313nmax./(1 + exp((Rc(a) - w314nmax)./w315nmax));
    v32nmin = w321nmin + w322nmin.*Rc(a) + w323nmin./(1 + exp((Rc(a) - w324nmin)./w325nmin));
    v32nmax = w321nmax + w322nmax.*Rc(a) + w323nmax./(1 + exp((Rc(a) - w324nmax)./w325nmax));
    v33nmin = w331nmin + w332nmin.*Rc(a) + w333nmin./(1 + exp((Rc(a) - w334nmin)./w335nmin));
    v33nmax = w331nmax + w332nmax.*Rc(a) + w333nmax./(1 + exp((Rc(a) - w334nmax)./w335nmax));
    v34nmin = w341nmin + w342nmin.*Rc(a) + w343nmin./(1 + exp((Rc(a) - w344nmin)./w345nmin));
    v34nmax = w341nmax + w342nmax.*Rc(a) + w343nmax./(1 + exp((Rc(a) - w344nmax)./w345nmax));
    v35nmin = w351nmin + w352nmin.*Rc(a) + w353nmin./(1 + exp((Rc(a) - w354nmin)./w355nmin));
    v35nmax = w351nmax + w352nmax.*Rc(a) + w353nmax./(1 + exp((Rc(a) - w354nmax)./w355nmax));

    t1nmin = v11nmin + v12nmin.*x + v13nmin.*x.^2 + v14nmin.*x.^3 + v15nmin.*x.^4;
    t1nmax = v11nmax + v12nmax.*x + v13nmax.*x.^2 + v14nmax.*x.^3 + v15nmax.*x.^4;
    t2nmin = v21nmin + v22nmin.*x + v23nmin.*x.^2 + v24nmin.*x.^3 + v25nmin.*x.^4;
    t2nmax = v21nmax + v22nmax.*x + v23nmax.*x.^2 + v24nmax.*x.^3 + v25nmax.*x.^4;
    t3nmin = v31nmin + v32nmin.*x + v33nmin.*x.^2 + v34nmin.*x.^3 + v35nmin.*x.^4;
    t3nmax = v31nmax + v32nmax.*x + v33nmax.*x.^2 + v34nmax.*x.^3 + v35nmax.*x.^4;

    phimunmin = Phimn.*(E + (t1nmin + t2nmin.*log10(E))./(Beta.^t3nmin)).^-alpha3;
    phimunmax = Phimn.*(E + (t1nmax + t2nmax.*log10(E))./(Beta.^t3nmax)).^-alpha3;

    g5n = h51n + h52n.*Rc(a) + h53n./(1 + exp((Rc(a) - h54n)./h55n));
    g6n = h61n + h62n.*Rc(a) + h63n./(1 + exp((Rc(a) - h64n)./h65n));

    f3n = g5n + g6n.*x;
    f2n = (phimunmin - phimunmax)./(smin.^f3n - smax.^f3n);
    f1n = phimunmin - f2n.*smin.^f3n;

    phimun = f1n + f2n.*s(a).^f3n;
    
    %Positive Muons
    v11pmin = w111pmin + w112pmin.*Rc(a) + w113pmin./(1 + exp((Rc(a) - w114pmin)./w115pmin));
    v11pmax = w111pmax + w112pmax.*Rc(a) + w113pmax./(1 + exp((Rc(a) - w114pmax)./w115pmax));
    v12pmin = w121pmin + w122pmin.*Rc(a) + w123pmin./(1 + exp((Rc(a) - w124pmin)./w125pmin));
    v12pmax = w121pmax + w122pmax.*Rc(a) + w123pmax./(1 + exp((Rc(a) - w124pmax)./w125pmax));
    v13pmin = w131pmin + w132pmin.*Rc(a) + w133pmin./(1 + exp((Rc(a) - w134pmin)./w135pmin));
    v13pmax = w131pmax + w132pmax.*Rc(a) + w133pmax./(1 + exp((Rc(a) - w134pmax)./w135pmax));
    v14pmin = w141pmin + w142pmin.*Rc(a) + w143pmin./(1 + exp((Rc(a) - w144pmin)./w145pmin));
    v14pmax = w141pmax + w142pmax.*Rc(a) + w143pmax./(1 + exp((Rc(a) - w144pmax)./w145pmax));
    v15pmin = w151pmin + w152pmin.*Rc(a) + w153pmin./(1 + exp((Rc(a) - w154pmin)./w155pmin));
    v15pmax = w151pmax + w152pmax.*Rc(a) + w153pmax./(1 + exp((Rc(a) - w154pmax)./w155pmax));
    v21pmin = w211pmin + w212pmin.*Rc(a) + w213pmin./(1 + exp((Rc(a) - w214pmin)./w215pmin));
    v21pmax = w211pmax + w212pmax.*Rc(a) + w213pmax./(1 + exp((Rc(a) - w214pmax)./w215pmax));
    v22pmin = w221pmin + w222pmin.*Rc(a) + w223pmin./(1 + exp((Rc(a) - w224pmin)./w225pmin));
    v22pmax = w221pmax + w222pmax.*Rc(a) + w223pmax./(1 + exp((Rc(a) - w224pmax)./w225pmax));
    v23pmin = w231pmin + w232pmin.*Rc(a) + w233pmin./(1 + exp((Rc(a) - w234pmin)./w235pmin));
    v23pmax = w231pmax + w232pmax.*Rc(a) + w233pmax./(1 + exp((Rc(a) - w234pmax)./w235pmax));
    v24pmin = w241pmin + w242pmin.*Rc(a) + w243pmin./(1 + exp((Rc(a) - w244pmin)./w245pmin));
    v24pmax = w241pmax + w242pmax.*Rc(a) + w243pmax./(1 + exp((Rc(a) - w244pmax)./w245pmax));
    v25pmin = w251pmin + w252pmin.*Rc(a) + w253pmin./(1 + exp((Rc(a) - w254pmin)./w255pmin));
    v25pmax = w251pmax + w252pmax.*Rc(a) + w253pmax./(1 + exp((Rc(a) - w254pmax)./w255pmax));
    v31pmin = w311pmin + w312pmin.*Rc(a) + w313pmin./(1 + exp((Rc(a) - w314pmin)./w315pmin));
    v31pmax = w311pmax + w312pmax.*Rc(a) + w313pmax./(1 + exp((Rc(a) - w314pmax)./w315pmax));
    v32pmin = w321pmin + w322pmin.*Rc(a) + w323pmin./(1 + exp((Rc(a) - w324pmin)./w325pmin));
    v32pmax = w321pmax + w322pmax.*Rc(a) + w323pmax./(1 + exp((Rc(a) - w324pmax)./w325pmax));
    v33pmin = w331pmin + w332pmin.*Rc(a) + w333pmin./(1 + exp((Rc(a) - w334pmin)./w335pmin));
    v33pmax = w331pmax + w332pmax.*Rc(a) + w333pmax./(1 + exp((Rc(a) - w334pmax)./w335pmax));
    v34pmin = w341pmin + w342pmin.*Rc(a) + w343pmin./(1 + exp((Rc(a) - w344pmin)./w345pmin));
    v34pmax = w341pmax + w342pmax.*Rc(a) + w343pmax./(1 + exp((Rc(a) - w344pmax)./w345pmax));
    v35pmin = w351pmin + w352pmin.*Rc(a) + w353pmin./(1 + exp((Rc(a) - w354pmin)./w355pmin));
    v35pmax = w351pmax + w352pmax.*Rc(a) + w353pmax./(1 + exp((Rc(a) - w354pmax)./w355pmax));

    t1pmin = v11pmin + v12pmin.*x + v13pmin.*x.^2 + v14pmin.*x.^3 + v15pmin.*x.^4;
    t1pmax = v11pmax + v12pmax.*x + v13pmax.*x.^2 + v14pmax.*x.^3 + v15pmax.*x.^4;
    t2pmin = v21pmin + v22pmin.*x + v23pmin.*x.^2 + v24pmin.*x.^3 + v25pmin.*x.^4;
    t2pmax = v21pmax + v22pmax.*x + v23pmax.*x.^2 + v24pmax.*x.^3 + v25pmax.*x.^4;
    t3pmin = v31pmin + v32pmin.*x + v33pmin.*x.^2 + v34pmin.*x.^3 + v35pmin.*x.^4;
    t3pmax = v31pmax + v32pmax.*x + v33pmax.*x.^2 + v34pmax.*x.^3 + v35pmax.*x.^4;

    phimupmin = Phimp.*(E + (t1pmin + t2pmin.*log10(E))./(Beta.^t3pmin)).^-alpha3;
    phimupmax = Phimp.*(E + (t1pmax + t2pmax.*log10(E))./(Beta.^t3pmax)).^-alpha3;

    g5p = h51p + h52p.*Rc(a) + h53p./(1 + exp((Rc(a) - h54p)./h55p));
    g6p = h61p + h62p.*Rc(a) + h63p./(1 + exp((Rc(a) - h64p)./h65p));

    f3p = g5p + g6p.*x;
    f2p = (phimupmin - phimupmax)./(smin.^f3p - smax.^f3p);
    f1p = phimupmin - f2p.*smin.^f3p;

    phimup = f1p + f2p.*s(a).^f3p;

    % Total Ground-Level Flux

    Phimu(a,:) = phimun + phimup;

%     clipindex = find(E <= 800, 1, 'last' );
    % clipindex = max(find(E <= 200));
    
    % Integral flux

%     mflux.total(a) = trapz(E(clipindex:end),Phimu(a,clipindex:end));
    % Differential fluxes
    mflux.neg(a,:) = phimun;
    mflux.pos(a,:) = phimup;
    % Integral fluxes for positive and negative muons
%     mflux.nint(a) = trapz(E(clipindex:end),phimun(clipindex:end));
%     mflux.pint(a) = trapz(E(clipindex:end),phimup(clipindex:end));
    mflux.E = E;
    mflux.p = p;
end

% 
% % Plot it
% 
% figure;clf;
% loglog(E,phimun);
% xlabel('E)'); ylabel('mu-');
% 
% figure;clf;
% loglog(E,phimup);
% xlabel('E'); ylabel('mu+');
% 
% figure;clf;
% loglog(E,Phimu); hold on;
% xlabel('E'); ylabel('muTotal');
% % ylim([0 .1]);
