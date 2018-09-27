%
% [N10,N26]=predN1026depth(pp,sp,sf,cp,age,depths,scaling_model)
%
% Predicts Al and Be concentrations for a point sample at a particular 
% depth (or vector of depths).
%
%  pp,sp,sf,cp            Chloe style physics, site, sample, and
%                         computed parameters.
%  age                    Exposure age in kyr.
%  depths                 Depths in g/cm^2.
%  scaling_model		  The scaling model being used
%  deltat                 Time step in ka (deltat=0.1 is 100 yrs)
%
%
function [N10,N26]=predN1026depth(pp,sp,sf,cp,age,depths,scaling_model,deltat)
%
% Figure out how many depths we're doing this for.
%
ndepths=length(depths);
%
% Get the erosion rate in gm/(cm^2*yr) (sp.epsilon is in gm/cm^2/kyr)  
%
erosionrate=sp.epsilon/1000;
%
% Note that from now on in this function, the erosion rate is in g/cm^2/yr.
%
%
% Adjust contemporary depth to depth at start time.
%
currentdepths=depths+erosionrate*age*1000;
%
% Set the maximum production depth being considered for the samples. This
% should be set to 200,000 g/cm2 or shallower.  
proddepth=200000;
%
% We use a time step of 100 years for the integration in time.  
%
%deltat=1;%0.1;
%
% We integrate until the point in time where the sample was collected.
%
tfinal=sp.tfinal;
%
% Figure out the depth spacing.
%
deltadepth=deltat*erosionrate*1000;
%
% Now, the main integration loop.
%
N10=zeros(ndepths,1);
N26=zeros(ndepths,1);
t=tfinal-age;
while (t+deltat < tfinal)
% Check to see if depths are unreasonable (negative/above surface)
    if (min(currentdepths-deltadepth/2)<0)
    fprintf(1,'currentdepth-deltadepth/2 is %f\n',currentdepths- ...
	    deltadepth/2);
    warning('Negative depth!');
    end
    %
    %
    if (min(currentdepths-deltadepth/2)<proddepth)
        % Update the elevation/latitude scaling factor to the current
        % time.  Note that our t is in kyr, with negative values in the
        % past, while the TDSF stuff is done in terms of years before present.
        %
        interptime=t+deltat/2;
        sf.currentsf=getcurrentsf(sf,interptime,scaling_model,'albe');
        
        %
        % Compute the production rate.  We use the mid point of the range
        % of depths corresponding to the current time interval.
        [pz10,pz26]=prodz1026(currentdepths-deltadepth/2,pp,sf,cp);
        %
        %
        % There are two terms here.  The first term is the old inventory of cosmogenic
        % nuclide, decayed for a time period of deltat.  The second term represents the
        % newly generated cosmogenic nuclide, including the radioactive decay of
        % some of it that was generated early in the time period.  The
        % radioactive decay factor is:
        %
        % f=int(P*exp(-r(deltat-t)),t=0..deltat)=P*(1-exp(-r*deltat))/r
        %
        % Note that for extremely small values of deltat you could get roundoff
        % errors when taking 1-exp(-r*deltat).  This hasn't been a problem for
        % our range of deltat values.
        %
        % The effect of using this term is that predNXX's results are essentially
        % indepedendent of deltat if the production rates are constant in time and
        % the erosion rate is 0.   If the erosion rate is nonzero, then deltat must
        % be small enough that very little erosion occurs during a time period, or
        % else predNXX's result will depend on deltat.
        %
        % Update N10 and N26.
        %
        f10=(1.0-exp(-pp.lambda10Be*deltat*1000))/pp.lambda10Be;
        N10=N10*exp(-pp.lambda10Be*deltat*1000)+...
            pz10*f10;
        
        f26=(1.0-exp(-pp.lambda26Al*deltat*1000))/pp.lambda26Al;
        N26=N26*exp(-pp.lambda26Al*deltat*1000)+...
            pz26*f26;
        
        %
        % Update t.
        %
        t=t+deltat;
        %
        % Update depth
        %
        currentdepths=currentdepths-deltat*1000*erosionrate;
    else
        %if samples are below production depth, calculate the new time and
        %depths for the samples to bring them into the production range. 
        
        depthskip=ceil((min(currentdepths-deltadepth/2)-proddepth)/deltadepth)*deltadepth;
        
        timeskip=depthskip/(1000*erosionrate);
        % Update t.
        %
        t=t+timeskip;
        %
        % Update depth
        %
        currentdepths=currentdepths-depthskip;
    end

end
%
% One final fractional time step.  deltatf is the length of this
% fractional time step.
%
deltatf=tfinal-t;
deltadepth=deltatf*erosionrate*1000;
%
% Update the elevation/latitude scaling factor to the current
% time.  Note that our t is in kyr, with negative values in the
% past, while the TDSF stuff is done in terms of years before present. 
%
interptime=t+deltatf/2;
sf.currentsf=getcurrentsf(sf,interptime,scaling_model,'albe');
%
% Compute the production rates.
[pz10,pz26]=prodz1026(currentdepths-deltadepth/2,pp,sf,cp);
%
% Update N10 and N26.
%
f10=(1.0-exp(-pp.lambda10Be*deltatf*1000))/pp.lambda10Be;
N10=N10*exp(-pp.lambda10Be*deltatf*1000)+...
    pz10*f10;

f26=(1.0-exp(-pp.lambda26Al*deltatf*1000))/pp.lambda26Al;
N26=N26*exp(-pp.lambda26Al*deltatf*1000)+...
    pz26*f26;
