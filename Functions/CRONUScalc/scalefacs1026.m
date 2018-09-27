%
% sf=scalefacs1026(sample,scaling_model)
%
function sf=scalefacs1026(sp1026,scaling_model)
%
% Check if scaling_model was specified, otherwise set it to default
%
if (~exist('scaling_model','var')), scaling_model = 'all'; end

%
% Setup the scale factors.
%
sf.ST=sp1026.ST;
sf.SLth=1;
sf.SLeth=1;
sf.P=sp1026.P;
sf.elevation=sp1026.elevation;

%
% Use Greg/Nat's code to compute the time dependent scale factors.
% The reason that we have a separate tdsfsample here is that the
% get_tdsf() function wasn't written by us.
%

load pmag_consts
tdsfsample.lat=sp1026.latitude;
tdsfsample.long=sp1026.longitude;
tdsfsample.pressure=sp1026.P;
tdsfsample.elevation=sp1026.elevation;
tdsfsample.scaling=scaling_model;
sf.tdsf=get_tdsf(tdsfsample,pmag_consts);
