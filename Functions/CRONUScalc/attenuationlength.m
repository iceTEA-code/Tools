function l=attenuationlength(lat,lon,elev,pressure)
%
% Uses the Sato model to compute an approximate attenuation length
% for spallation.  This attenuation length is typically within
% about 10% of attenuation lengths that have been obtained by
% fitting to depth profiles.  It should generally be OK for surface
% samples, where a small error in the attenuation length will not
% be significant.  
%
% The attenuation length may need to be further corrected for
% topography.
%
% The attenuation length is computed based on the long term average
% rigidity cutoff for the site.  This will not necessarily be
% appropriate for relatively young sites.  
%
% Input:
%    lat                  Latitude (degrees)
%    lon                  Longitude (degrees east)
%    elev                 Elevation (m)
%    pressure             Atmospheric pressure (hPa)
%
% Output:
%    l                    Attenuation length (g/cm^2)

%
% Construct a fake 1026 sample so that we can call get pars and
% obtain the historical rigidity cutoffs for this site. 
%
sample=zeros(15,1);
sample(1)=lat;
sample(2)=lon;
sample(3)=elev;
sample(4)=pressure;
sample(5)=3;
sample(6)=2.66;
sample(7)=1.0;
sample(8)=0.0;
sample(9)=0.0;
sample(10)=0.0;
sample(11)=0.0;
sample(12)=0.0;
sample(13)=150;
sample(14)=0;
sample(15)=2010;
%
% Run get pars to get geomag information for this site.
%
[pp,sp,sf,cp]=getpars1026(sample,'sa');
meanRc=mean(sf.tdsf.Rc_Sa);
l=rawattenuationlength(pressure,meanRc);