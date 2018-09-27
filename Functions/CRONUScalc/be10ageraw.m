%
%  age=be10ageraw(sampledata,pp,sf,scaling_model)
%
%  Given the data for a sample computes the age of the sample.
%
%
% This inner routine does not handle the uncertainty calculations, 
% which are done by be10age.m.  Instead, this inner routine simply 
% does the basic computation of the age.
%
% The sampledata1026 vector contains the following information:
%
%1. Latitude (decimal degrees, -90(S) to +90(N))
%2. Longitude (decimal degrees, 0-360 degrees east)
%3. Elevation (meters)
%4. Pressure (hPa)      
%5. sample thickness (cm)
%6. bulk density (g/cm^3)
%7. Shielding factor for terrain, snow, etc. (unitless)
%8. Erosion-rate epsilon (g/(cm^2*kyr))
%9. Sample 10-Be concentration (atoms of 10-Be/g of target)
%10. Sample 26-Al concentration (atoms of 26-Al/g of target)
%11. Inheritance for Be (atoms 10-Be/g of target)
%12. Inheritance for Al (atoms 26-Al/g of target)
%13. Effective attenuation length -Lambdafe (g/cm^2)
%14. Depth to top of sample (g/cm^2)
%15. Year sampled (e.g. 2010)
%
% Also requires physical parameters (pp) and scaling factors (sf).
%
% Returns an output vector:
%
% 1. age (kyr)
% 2. age uncertainty (always 0 in this routine)
% 3. Contemporary Elevation/latitude scaling factor for neutrons for Be (unitless)
% 4. Contemporary Elevation/latitude scaling factor for fast muons (unitless)
% 5. Contemporary Elevation/lat scaling factor for slow muons (unitless)
% 6. Contemporary depth avg prod rate, neutron spallation (atoms/g/yr)
% 7. Contemporary depth avg prod rate, muons (atoms/g/yr)
% 8. Qs (unitless)
% 9. Qmu (unitless)
% 10. Inherited 10-Be (atoms/g of target)
% 11. Measured 10-Be (atoms/g of target)
% 12. Internal (analytical only) uncertainty (always 0 in this routine)
%
function output=be10ageraw(sampledata,pp,sf,scaling_model)
%
% Make sampledata a column vectors if it isn't already.
%
if (size(sampledata,1)==1)
  sampledata=sampledata';
end
%
% First, check that the input data is reasonable.
%
if (length(sampledata) ~= 15)
  error('sampledata has wrong size!');
end
%
% Setup the values of sample parameters.
%
sp=samppars1026(sampledata);
%
% Figure out the maximum possible depth at which we'll ever need a
% production rate.  This is depthtotop + maxage * erosion +
% thickness * density + a safety factor.
%
maxage=8160;        % 6* half-life = 8160 ka 
maxdepth=sp.depthtotop+maxage*sp.epsilon+sp.ls*sp.rb+1000; 
%
% Computed parameters.  
%
cp=comppars1026(pp,sp,sf,maxdepth);
%
% Compute the age. 
%
age=computeage10(pp,sp,sf,cp,scaling_model);
%
% Next, compute production rates for various pathways.
% Here, we call prodz once more at time 0 (present) and then
% compute the depth average from the results.
%
%
% In doing this computation, use the contemporary production rate.
%
sf.currentsf=getcurrentsf(sf,0,scaling_model,'be');
%
%
%
nlayers=100;
thickness=sp.ls*sp.rb;
thick=(thickness/nlayers)*ones(nlayers,1);
totalthickness=sum(thick);
midpt=sp.depthtotop+(thickness/nlayers)*((1:nlayers)'-0.5);
[PtotalBe,PtotalAl,PsBe,PmuBe,PsAl,PmuAl]=prodz1026(midpt,pp,sf,cp);
srate=(thick'*PsBe)/totalthickness;
murate=(thick'*PmuBe)/totalthickness;
Qs=srate/PsBe(1);
Qmu=murate/PmuBe(1);
%
% Setup the output vector.
%
output=[age; 0; sf.currentsf.Sel10; cp.SFmufast; cp.SFmuslow; srate; murate; Qs; Qmu; sp.inheritance10;...
	sp.concentration10; 0];



