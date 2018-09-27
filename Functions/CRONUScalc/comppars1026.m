
%
% cp=comppars1026(physpars,samppars,scalefactors,maxdepth)
%
% Generates a complete set of parameters for beryllium-10/Al-26 dating
% from the basic physics parameters, sample specific parameters, 
% and scale factors.
%
% maxdepth is an optional parameter.  The default value is 2500
% g/cm^2, or about 10 meters.  Computing muon production profiles
% down this far is extremely expensive.  For surface samples, a max
% depth of 25 g/cm^2 (or about 10cm)  is more appropriate.  
%
function cp=comppars1026(pp,sp1026,sf1026,maxdepth)
%
% First, set maxdepth to a default value if not supplied.
%
if (nargin < 4)
  maxdepth=2500;
end

%
% Set production depth for the sample. This defaults to 200000 because
% production equations change at this depth. 150000 is usually safe, even
% for high erosion samples. 
proddepth=200000;

% Record the maximum depth so that we know how deep it is safe to go.
%
if maxdepth<proddepth;
    cp.maxdepth=maxdepth;
else 
    cp.maxdepth=proddepth;
end

%
% Setup Lambdafe.  
%
cp.Lambdafe=sp1026.Lambdafe;
%

% Zs
%
cp.Zs=sp1026.ls*sp1026.rb;

%
%---------------------------------------------------------------------
%
% The following are modifications for Sato/Heisinger muons.  
%
%
% Pmunf
%
% old formulation: 
%  cp.Pmunf=0.0000058*pp.Yf*pp.Phimuf0*exp(-cp.Zs/pp.Lambdamu);
%
%  New Formulation uses Greg's Heisinger code to calculate the fluxes at 
% a vector of depths to create a table that can be referred to later

% RcEst = 14.9.*((cos(abs(d2r(sp1026.latitude)))).^4); %Dipolar estimate for these purposes
% The RcEst above was originally used and has been replaced by Nat Lifton's
% new model.  He  fit the function below to trajectory-traced cutoffs for the
% modern dipolar field strength, so it should be more accurate overall.

dd = [6.89901,-103.241,522.061,-1152.15,1189.18,-448.004;];
%   Trajectory-traced dipolar estimate for these purposes
   RcEst = (dd(1)*cos(d2r(sp1026.latitude)) + ...
       dd(2)*(cos(d2r(sp1026.latitude))).^2 + ...
       dd(3)*(cos(d2r(sp1026.latitude))).^3 + ...
       dd(4)*(cos(d2r(sp1026.latitude))).^4 + ...
       dd(5)*(cos(d2r(sp1026.latitude))).^5 + ...
       dd(6)*(cos(d2r(sp1026.latitude))).^6);
   
%     % %Use mean Solar Modulation Parameter (SPhiInf)
nnn=ceil(maxdepth/5);
%
deltadepth=maxdepth/nnn;    % this gives the step for each

cp.depthvector=[0:deltadepth:maxdepth];
%
%store the output fluxes that we need

%Also store the muons production rates from the code
%
% Preallocate space for cp.muon1026.  
%
cp.muon1026=zeros(3,length(cp.depthvector));

%
% Store the depth vector.
%
cp.muon1026(1,:)=cp.depthvector;

%
% Pmun0
%
% Old Formulation (eqn 3.49 in G&P)
%   cp.Pmun0=pp.Ys*pp.Psimu0+0.0000058*pp.Yf*pp.Phimuf0;
%
%  New Formulation uses Greg's Heisinger table from above to find the 
%  surface production of neutrons for that particular elevation.  

%cp.Pmun0=pp.Ys*(cp.negflux(1))+0.0000058*(cp.totalflux(1));
%

% Generate Be-10 production rates if there is Be-10 data.
%
if sp1026.concentration10 > 0

  %set the constants to Be-10
  Natoms = pp.Natoms10;
  %  sigma190 = pp.sigma190_10;
  %    pp.delsigma190 = pp.delsigma190_10; % not used
  sigma0=pp.sigma010;
  k_negpartial = pp.k_negpartial10;
  fstar=pp.fstar10;
  %    pp.delk_neg = pp.delk_neg10; % not used
  bemuon=muonfluxsato(cp.depthvector,sf1026.P,RcEst,pp.SPhiInf,pp,'yes');
  % Now calculate the production rates. 
  z=cp.depthvector;

  %store muon scaling factor
  cp.SFmufast=bemuon.SFmufast;
  cp.SFmuslow=bemuon.SFmuslow;
  
  % Depth-dependent parts of the fast muon reaction cross-section
  % Balco original - from Heisinger fast muon paper Sea Level
  % Beta = 0.846 - 0.015 .* log((z./100)+1) + 0.003139 .* (log((z./100)+1).^2);
  
  % % For Beacon heights
  % aalpha = 0.75;
  % Beta =  0.842344 - 0.0107398 log((z./100)+1) + 0.00240182 log((z./100)+1)^2
  % 
  % aalpha = 0.85;
  % Beta =  0.888695 - 0.00716992 log((z./100)+1) + 0.00169676 log((z./100)+1)^2
  % 
  aalpha = 1.0;
  Beta = 1.0;
  % 
  % alpha = 1.15;
  % Beta =  1.18129 + 0.00903804 log((z./100)+1) - 0.00273586 lnlog((z./100)+1)^2
  % 
  % aalpha = 1.30;
  % Beta =  1.45283+0.04615 lnlog((z./100)+1) - 0.0153481 lnlog((z./100)+1)^2 + 0.000672339 lnlog((z./100)+1)^3
  % % 
  
  Ebar = 7.6 + 321.7.*(1 - exp(-8.059e-6.*z)) + 50.7.*(1-exp(-5.05e-7.*z));
  %sigma0 = consts.sigma190./(190.^aalpha);
  
  % fast muon production
  
  P_fast = bemuon.phi.*Beta.*(Ebar.^aalpha).*sigma0.*Natoms;
  
  % negative muon capture
  P_neg = bemuon.R.*k_negpartial.*fstar;
  
  cp.P_total10=P_fast+P_neg;
  
  cp.muon1026(2,:)=cp.P_total10;
  cp.negflux=bemuon.R;
  cp.totalflux=bemuon.phi;
  cp.hasbe=1;
else
  cp.hasbe=0;
end

%
% Generate Al-26 production rates if there is data.
%
%check to see if there is Al data
if sp1026.concentration26 > 0
  %set the constants to Al-26
  Natoms = pp.Natoms26;
  %  sigma190 = pp.sigma190_26;
  %    pp.delsigma190 = pp.delsigma190_10; % not used
  sigma0=pp.sigma026;
  k_negpartial = pp.k_negpartial26;
  fstar=pp.fstar26;
  %    pp.delk_neg = pp.delk_neg10; % not used
  almuon=muonfluxsato(cp.depthvector,sf1026.P,RcEst,pp.SPhiInf,pp,'yes');
  % Now calculate the production rates. 
  z=cp.depthvector;
  % Depth-dependent parts of the fast muon reaction cross-section
  % Balco original - from Heisinger fast muon paper Sea Level
  % Beta = 0.846 - 0.015 .* log((z./100)+1) + 0.003139 .* (log((z./100)+1).^2);
  
  % % For Beacon heights
  % aalpha = 0.75;
  % Beta =  0.842344 - 0.0107398 log((z./100)+1) + 0.00240182 log((z./100)+1)^2
  % 
  % aalpha = 0.85;
  % Beta =  0.888695 - 0.00716992 log((z./100)+1) + 0.00169676 log((z./100)+1)^2
  % 
  aalpha = 1.0;
  Beta = 1.0;
  % 
  % alpha = 1.15;
  % Beta =  1.18129 + 0.00903804 log((z./100)+1) - 0.00273586 lnlog((z./100)+1)^2
  % 
  % aalpha = 1.30;
  % Beta =  1.45283+0.04615 lnlog((z./100)+1) - 0.0153481 lnlog((z./100)+1)^2 + 0.000672339 lnlog((z./100)+1)^3
  % % 

  Ebar = 7.6 + 321.7.*(1 - exp(-8.059e-6.*z)) + 50.7.*(1-exp(-5.05e-7.*z));
  %sigma0 = consts.sigma190./(190.^aalpha);

  % fast muon production

  P_fast = almuon.phi.*Beta.*(Ebar.^aalpha).*sigma0.*Natoms;

  % negative muon capture
  P_neg = almuon.R.*k_negpartial.*fstar;

  cp.P_total26=P_fast+P_neg;
  
  cp.muon1026(3,:)=cp.P_total26;
  
  cp.hasal=1;
  if cp.hasbe == 0
    cp.negflux=almuon.R;
    cp.totalflux=almuon.phi;
    %store muon scaling factor
    cp.SFmufast=almuon.SFmufast;
    cp.SFmuslow=almuon.SFmuslow;
  end
else
  cp.hasal=0;
end


