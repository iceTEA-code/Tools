%
%  [output,times,plotprodbe,derivs]=be10age(sampledata,sampleuncertainties,scaling_model)
%
%  Given the data for a sample and associated one standard
%  deviation uncertainties, computes the age of the sample and uncertainty.
%
% The sampledata vector contains the following information:
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
% A second input vector, sampleuncertainties, contains 1-sigma
% uncertainties for all 15 inputs.  In general, we assume that
% these 15 inputs are uncorrelated. 
%
% scaling_model is one of 'DE','DU','LI','LM','SA','SF','ST' and 
% informs which scaling model is being used
%
% Returns an output vector:
%
% 1. Age (kyr)
% 2. Total (external) age uncertainty (kyr) 
% 3. Contemporary Elevation/latitude scaling factor for neutrons for Be (unitless)
% 4. Contemporary Elevation/latitude scaling factor for fast muons (unitless)
% 5. Contemporary Elevation/lat scaling factor for slow muons (unitless)
% 6. Contemporary depth avg prod rate, neutron spallation (atoms/g/yr)
% 7. Contemporary depth avg prod rate, muons (atoms/g/yr)
% 8. Qs (unitless)
% 9. Qmu (unitless)
% 10. Inherited 10-Be (atoms/g of target)
% 11. Measured 10-Be (atoms/g of target)
% 12. Analytical only (internal) uncertainty (kyr)
%
% An optional output that can be useful for debugging is the vector
% derivs, which contains the derivatives of the age with respect to
% the 15 input parameters.  Note that the derivatives are only
% computed if the corresponding uncertainty is nonzero.  
%
function [output,times,plotprodbe,derivs]=be10age(sampledata,sampleuncertainties,scaling_model)
  expectedOutputs = 12;
  %
  % Make sampledata and uncertainties column vectors if they aren't already.
  %
  if (size(sampledata,1)==1)
    sampledata=sampledata';
  end
  if (size(sampleuncertainties,1)==1)
    sampleuncertainties=sampleuncertainties';
  end
  
  %
  % First, check that the input data is reasonable.
  %
  if (length(sampledata) ~= 15)
    error('sampledata has wrong size!');
  end
  if (length(sampleuncertainties) ~= 15)
    error('sampleuncertainties has wrong size!');
  end
  %
  % Setup the physical parameters.
  %
  pp=physpars(scaling_model);
  %
  % Extract the sample parameters from the sampledatavector.
  %
  sp=samppars1026(sampledata);
  %
  % Get the scale factors.
  %
  sf=scalefacs1026(sp,scaling_model);
  %
  % We need an absolute maximum age for several purposes, including
  % detecting saturated samples and setting the maximum depth for comppars.
  %
  maxage=8160;        % 6* half-life = 8160 ka                            
  %
  % Figure out the maximum possible depth at which we'll ever need a
  % production rate.  This is depthtotop + maxage * erosion (g/cm^2/kyr)+
  % thickness * density + a safety factor.
  %
  maxdepth=sp.depthtotop+maxage*sp.epsilon+sp.ls*sp.rb+1000;
  
  %check for maximum depth: muon formulation is only good down to 2e5 g/cm2.
  % This will likely never happen, but if it does, Matlab will error out saying we have not
  % supplied times and plotprod and we will receive that error and be able to help the user.
  if maxdepth > 2e5
    fprintf(1,'Maximum sample depth (%f) exceeds muon formulation of 2e5 g/cm2. \n Options: Lower the erosion rate, lower the maxage in file be10age, or change muon formulation in muonfluxsato.m',[maxdepth])
    warning('This sample exceeds muon maximum depth!');
    output=NaN*ones(expectedOutputs,1);
    [ times, plotprodbe ] = getPlotData(maxage, sf, 0, 'be', scaling_model);
    return;
  end
  %
  % Computed parameters.
  %
  cp=comppars1026(pp,sp,sf,maxdepth);
  %
  % Get contemporary surface production rates in atoms/g 
  %
  sf.currentsf=getcurrentsf(sf,0,scaling_model,'be');
  [ProdtotalBe,ProdtotalAl,ProdsBe,ProdmuBe,ProdsAl,ProdmuAl]=prodz1026(0,pp,sf,cp);
  %
  % Before we go any further, check to make sure that this sample
  % isn't saturated.  If the measured concentration is greater than
  % 95% of the saturated concentration then return NaN's.
  %
  timeStepForAging=0.1;
  [satN10,satN26]=predN1026(pp,sp,sf,cp,maxage,scaling_model,timeStepForAging);
  if (sp.concentration10 > 0.95*satN10)
    fprintf(1,'Saturated N10=%f, Measured N10=%f \n',[satN10; sp.concentration10]);
    warning('This sample is saturated!');
    output=NaN*ones(expectedOutputs,1);
    [ times, plotprodbe ] = getPlotData(maxage, sf, ProdsBe, 'be', scaling_model);
    return;
  end
  %
  % Compute the nominal results.
  %
  output=be10ageraw(sampledata,pp,sf,scaling_model);
  age=output(1);
  
  [ times, plotprodbe ] = getPlotData(age, sf, ProdsBe, 'be', scaling_model);
  
  %
  % Start off with a sum of 0.
  %
  uncertainty=0.0;
  %
  % Work through all 15 sample parameters.  For each nonzero
  % uncertainty, add in the uncertainty term.
  %
  derivs=zeros(15,1);
  for i=1:15
    if ((sampleuncertainties(i) ~= 0.0) | (nargout > 3))
      if (sampledata(i) ~= 0.0)
        thisdelta=0.01*abs(sampledata(i));
      else
        thisdelta=0.01;
      end
      deltasampledata=sampledata;
      deltasampledata(i)=deltasampledata(i)+thisdelta;
      deltasp=samppars1026(deltasampledata);
      deltasf=scalefacs1026(deltasp,scaling_model);
      deltaoutput=be10ageraw(deltasampledata,pp,deltasf,scaling_model);
      deltaage=deltaoutput(1);
      dagedi=(deltaage-age)/(thisdelta);
      derivs(i)=dagedi;
      if (~isnan(sampleuncertainties(i)))
        uncertainty=uncertainty+(dagedi^2*sampleuncertainties(i)^2);
      end
    end
  end
  %
  % Add in terms for the uncertainty in production rates.  
  %
  %
  % Uncertainty in PsBe.
  %
  deltapp=pp;
  deltapp.PsBe=pp.PsBe+0.01*abs(pp.PsBe);
  deltaoutput=be10ageraw(sampledata,deltapp,sf,scaling_model);
  deltaage=deltaoutput(1);
  dagedpsBe=(deltaage-age)/(0.01*abs(pp.PsBe));
  uncertaintyext=uncertainty+(dagedpsBe^2*pp.sigmaPsBe^2);
  %
  % Finally, take the square root of uncertainty to get a standard deviation.
  %
  %This is the internal uncertainty (analytical only)
  uncertainty=sqrt(uncertainty);
  %This is external/total uncertainty
  uncertaintyext=sqrt(uncertaintyext);
  %
  % Put the uncertainty in the output.
  %
  output(2)=uncertaintyext;
  output(12)=uncertainty;

  if size(output,1) ~= expectedOutputs
    error('Be outputs not expected size');
  end
end
