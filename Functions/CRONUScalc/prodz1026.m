%
% [PtotalBe,PtotalAl,PsBe,PmuBe,PsAl,PmuAl]=prodz1026(z,pp,sf,cp)
%
%
% The prodz function calculates the individual contribution of Al and
% Be to the production from spallation and muon pathways 
% at a particular given depth. 

%  pp,sf,cp               physics, site, sample, and
%                         computed parameters.
%  z                      Depths in g/cm^2.
% 
% Get current scaling factors.
% Then call this function to compute total production for a vector of
% depths.  
% Note: This function has been updated to work with nuclide-dependent
% scaling factors (only applicable to the Sato/Lifton scaling scheme).
%
function [ProdtotalBe,ProdtotalAl,ProdsBe,ProdmuBe,ProdsAl,ProdmuAl]=prodz1026(z,pp,sf,cp)
%
% Make sure that we don't have any depths that are too big.
%
if (max(z) > cp.maxdepth)
  error('Prodz called for depth greater than maxdepth.');
end
%
% Make sure that sf.SelXX is a number.
%
if (isnan(sf.currentsf.Sel10))
    error('sf.currentsf.Sel10 is a NaN');
end
if (isnan(sf.currentsf.Sel26))
    error('sf.currentsf.Sel26 is a NaN');
end
%
% We'll get some slightly negative depths due to roundoff errors.
% Fix them.
%
for i=1:length(z)
  if (z(i)<-1.0e-4)
    error('z(i) is negative!');
  end
  if (z(i)<0.0)
    z(i)=0.0;
  end
end

%
%find the number of depths given in the input vector z
numberdepths=length(z);
%
%For each z, find the appropriate negative muon flux and total flux terms by
%interpolation
%assume depths are always given in an already-sorted vector

negfluxdepth=interpolate(cp.depthvector,cp.negflux,z);
totalfluxdepth=interpolate(cp.depthvector,cp.totalflux,z);


% Ps(z)=Sel*ST*Ps_0*exp(-z/BigLambda_f)
%
% Introduce a new term to speed up computation
%
expfactor=exp(-z/cp.Lambdafe);
%
% New productions for each spallation pathway separately, then deal with
% muons

ProdsBe=sf.currentsf.Sel10*sf.ST*pp.PsBe*expfactor;
%
% Muons are now being multiplied by the terrain shielding factor as well as
% the latitude-dependent muon scaling factor.  
%
ProdmuBe=sf.ST*interpolate(cp.muon1026(1,:),cp.muon1026(2,:),z);
if ((sum(isnan(ProdmuBe)) > 0) & (cp.hasbe==1))
  ProdmuBe(isnan(ProdmuBe))=0;
end

%
% Now, compute the total production.
%
ProdtotalBe=ProdsBe+ProdmuBe;
%
% Look for NaN's in the results.
%
if ((sum(isnan(ProdtotalBe)) > 0) & (cp.hasbe==1))
  ProdsBe
  ProdmuBe
  error('Prodz1026 produced NaN in Be!');
end

% check to see if there is an aluminum concentration.  If there is, then
% calculate production rates for Al (spallation and muon)
%
if (cp.hasal==1)
    ProdsAl=sf.currentsf.Sel26*sf.ST*pp.PsAl*expfactor;
    ProdmuAl=sf.ST*interpolate(cp.muon1026(1,:),cp.muon1026(3,:),z);
    %compute total production
    ProdtotalAl=ProdsAl+ProdmuAl;
    if (sum(isnan(ProdtotalAl)) > 0) 
      ProdsAl
      ProdmuAl
      error('Prodz1026 produced NaN in Al!');
    end
else
    ProdsAl=NaN;
    ProdmuAl=NaN;
    ProdtotalAl=NaN;
end


