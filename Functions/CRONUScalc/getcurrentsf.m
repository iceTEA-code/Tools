%
% out=getcurrentsf(sf,t,scaling_model,nuclide);
%
% Get the elevation/latitude time dependent scaling factors for
% neutrons by nuclide (SelXX) and thermal (SFth) and epithermal neutrons (SFeth) 
% by interpolating the tables of time dependent scaling factors.  
% This file needs to be edited to select a particular scaling scheme.  For
% scaling schemes that do not compute nuclide-dependent scaling factors,
% all the scaling factors are set equal to the Sel10 (which is the Sel
% factor for all the non-muon reactions in the codes). 
%
function out=getcurrentsf(sf,t,scaling_model,nuclide);
%
% Handle extremely old or extremely young times. 
%
if (t<-max(sf.tdsf.tv)/1000)
    t=-max(sf.tdsf.tv)/1000;
    warning('extremely old time in getcurrentsf!');
end

%
% If scaling_model is not specified, error
%
if (~exist('scaling_model','var'))
	error('CRONUS:InternalError','No scaling_model supplied to getcurrentsf');
else
	scaling_model = lower(scaling_model);
end

%
% If nuclide is not specified, default to 'all'
%
if (~exist('nuclide','var'))
	nuclide = 'all'; 
else
	nuclide = lower(nuclide);
end

%
% Look for times after 2110 AD.  This is most likely an error.
% 
if (t>0.3)
  t
  warning('extremely young time in getcurrentsf!');
  t=0;
end
%
% For all times after t=0 (2010), use the 2010 scaling factors by
% forcing t=0.  We don't let this get too far out of hand- this is
% limited to 2110 (see above.)
%
if (t > 0)
  t=0.0;
end

% All of the following interpolate calls use a converted value for t
% so we do this now to cut down on math operations (this is called
% many many times)
new_t = -t*1000;

switch scaling_model
case 'de'
%
% This version for the Desilets and Zreda scheme.
% No point in checking nuclide since it's the same 
% calculation regardless
%
	out.Sel10=interpolate(sf.tdsf.tv,sf.tdsf.SF_De,new_t);
	out.SelSF=out.Sel10;
    out.Sel26=out.Sel10;
	out.Sel14=out.Sel10;
	out.Sel3=out.Sel10;
	out.Sel36Ca=out.Sel10;
	out.Sel36K=out.Sel10;
	out.Sel36Ti=out.Sel10;
	out.Sel36Fe=out.Sel10;
	out.Sel21=out.Sel10;
	out.SFth=out.Sel10;
	out.SFeth=out.Sel10;
	
case 'du'
%
% This version for the Dunai scheme.
% No point in checking nuclide since it's the same 
% calculation regardless
%
	out.Sel10=interpolate(sf.tdsf.tv,sf.tdsf.SF_Du,new_t);
	out.SelSF=out.Sel10;
    out.Sel26=out.Sel10;
	out.Sel14=out.Sel10;
	out.Sel3=out.Sel10;
	out.Sel36Ca=out.Sel10;
	out.Sel36K=out.Sel10;
	out.Sel36Ti=out.Sel10;
	out.Sel36Fe=out.Sel10;
	out.Sel21=out.Sel10;
	out.SFth=out.Sel10;
	out.SFeth=out.Sel10;
	
case 'li'
%
% This version for the Lifton scheme.
% No point in checking nuclide since it's the same 
% calculation regardless
%
	out.Sel10=interpolate(sf.tdsf.tv,sf.tdsf.SF_Li,new_t);
	out.SelSF=out.Sel10;
    out.Sel26=out.Sel10;
	out.Sel14=out.Sel10;
	out.Sel3=out.Sel10;
	out.Sel36Ca=out.Sel10;
	out.Sel36K=out.Sel10;
	out.Sel36Ti=out.Sel10;
	out.Sel36Fe=out.Sel10;
	out.Sel21=out.Sel10;
	out.SFth=out.Sel10;
	out.SFeth=out.Sel10;
	
case 'lm'
%
% This version for the time dependent Lal scheme.
% No point in checking nuclide since it's the same 
% calculation regardless
%
	out.Sel10=interpolate(sf.tdsf.tv,sf.tdsf.SF_Lm,new_t);
	out.SelSF=out.Sel10;
    out.Sel26=out.Sel10;
	out.Sel14=out.Sel10;
	out.Sel3=out.Sel10;
	out.Sel36Ca=out.Sel10;
	out.Sel36K=out.Sel10;
	out.Sel36Ti=out.Sel10;
	out.Sel36Fe=out.Sel10;
	out.Sel21=out.Sel10;
	out.SFth=out.Sel10;
	out.SFeth=out.Sel10;
	
case 'st'
%
% This version for the Lal/Stone scheme.  
% No point in checking nuclide since it's the same 
% calculation regardless
%
    if length(sf.tdsf.SF_St) < length(sf.tdsf.tv) % MODIFICATION
        out.Sel10=sf.tdsf.SF_St;
    else
        out.Sel10=interpolate(sf.tdsf.tv,sf.tdsf.SF_St,new_t);
    end
	out.SelSF=out.Sel10;
    out.Sel26=out.Sel10;
	out.Sel14=out.Sel10;
	out.Sel3=out.Sel10;
	out.Sel36Ca=out.Sel10;
	out.Sel36K=out.Sel10;
	out.Sel36Ti=out.Sel10;
	out.Sel36Fe=out.Sel10;
	out.Sel21=out.Sel10;
	out.SFth=out.Sel10;
	out.SFeth=out.Sel10;
	
case 'sa'
%
% This version for the Sato nuclide-dependent scheme.  Implements nuclide-dependent scaling
% by incorporating cross-sections.
%
if (strcmpi(nuclide,'all'))
	out.Sel10=interpolate(sf.tdsf.tv,sf.tdsf.SF_Sa10,new_t);
    out.SelSF=interpolate(sf.tdsf.tv,sf.tdsf.SF_Sf,new_t);
	out.Sel21=out.Sel10;
	out.Sel26=interpolate(sf.tdsf.tv,sf.tdsf.SF_Sa26,new_t);
	out.Sel14=interpolate(sf.tdsf.tv,sf.tdsf.SF_Sa14,new_t);
	out.Sel3=interpolate(sf.tdsf.tv,sf.tdsf.SF_Sa3,new_t);
	out.Sel36Ca=interpolate(sf.tdsf.tv,sf.tdsf.SF_Sa36Ca,new_t);
	out.Sel36K=interpolate(sf.tdsf.tv,sf.tdsf.SF_Sa36K,new_t);
	out.Sel36Ti=interpolate(sf.tdsf.tv,sf.tdsf.SF_Sa36Ti,new_t);
	out.Sel36Fe=interpolate(sf.tdsf.tv,sf.tdsf.SF_Sa36Fe,new_t);
	out.SFth=interpolate(sf.tdsf.tv,sf.tdsf.SF_Sath,new_t);
	out.SFeth=interpolate(sf.tdsf.tv,sf.tdsf.SF_Saeth,new_t);
elseif (strcmpi(nuclide,'be') || strcmpi(nuclide,'al') || strcmpi(nuclide,'albe'))
	out.Sel10=interpolate(sf.tdsf.tv,sf.tdsf.SF_Sa10,new_t);
	out.Sel26=interpolate(sf.tdsf.tv,sf.tdsf.SF_Sa26,new_t);
elseif (strcmpi(nuclide,'c'))
	out.Sel14=interpolate(sf.tdsf.tv,sf.tdsf.SF_Sa14,new_t);
elseif (strcmpi(nuclide,'he'))
	out.Sel3=interpolate(sf.tdsf.tv,sf.tdsf.SF_Sa3,new_t);
elseif (strcmpi(nuclide,'ne'))
	%same as out.Sel10
	out.Sel21=interpolate(sf.tdsf.tv,sf.tdsf.SF_Sa10,new_t);
elseif (strcmpi(nuclide,'cl'))
	out.Sel36Ca=interpolate(sf.tdsf.tv,sf.tdsf.SF_Sa36Ca,new_t);
	out.Sel36K=interpolate(sf.tdsf.tv,sf.tdsf.SF_Sa36K,new_t);
	out.Sel36Ti=interpolate(sf.tdsf.tv,sf.tdsf.SF_Sa36Ti,new_t);
	out.Sel36Fe=interpolate(sf.tdsf.tv,sf.tdsf.SF_Sa36Fe,new_t);
	out.SFth=interpolate(sf.tdsf.tv,sf.tdsf.SF_Sath,new_t);
	out.SFeth=interpolate(sf.tdsf.tv,sf.tdsf.SF_Saeth,new_t);
    out.SelSF=interpolate(sf.tdsf.tv,sf.tdsf.SF_Sf,new_t);
end
%
%uncomment and change below if we go back to time-dependent muons
%
%SFmufast=interpolate(sf.tdsf.tv,sf.tdsf.SFmu_Li,new_t);
%SFmuslow=interpolate(sf.tdsf.tv,sf.tdsf.SFmu_Li,new_t);
	
case 'sf'
%
% This version is for the Sato flux scaling scheme (NOT nuclide dependent) 
% No point in checking nuclide since it's the same 
% calculation regardless
%
	out.Sel10=interpolate(sf.tdsf.tv,sf.tdsf.SF_Sf,new_t);
	out.SelSF=out.Sel10;
    out.Sel26=out.Sel10;
	out.Sel14=out.Sel10;
	out.Sel3=out.Sel10;
	out.Sel36Ca=out.Sel10;
	out.Sel36K=out.Sel10;
	out.Sel36Ti=out.Sel10;
	out.Sel36Fe=out.Sel10;
	out.Sel21=out.Sel10;
	out.SFth=out.Sel10;
	out.SFeth=out.Sel10;
end