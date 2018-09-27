%
% [N10,N26]=predN1026(pp,sp,sf,cp,age,scaling_model)
%
% Predicts N10 and N26 (in atoms/gm) given the parameters and the age in ka.
% This version has been updated to deal with samples at depth (see
% the sp.depthtotop parameter.)  It make use of predN1026depth(). deltat is
% the time step in ka (deltat of 0.1 is 100 years). 
%
function [N10,N26]=predN1026(pp,sp,sf,cp,age,scaling_model,deltat)
%
% Convert the sample thickness to g/cm^2.
%
thickness=sp.ls*sp.rb;
%
% Set the number of depths to interpolate production across the
% thickness of the sample.
%
ndepths=10;
%
% Work out the specific depths at which to compute cumulative
% production.  
%
depths=zeros(ndepths,1);
deltadepth=thickness/ndepths;
for i=1:ndepths
  depths(i)=sp.depthtotop+deltadepth/2+deltadepth*(i-1);
end
%
% Compute N10 and N26 at each depth.
%
N10depths=zeros(ndepths,1);
N26depths=zeros(ndepths,1);
[N10depths,N26depths]=predN1026depth(pp,sp,sf,cp,age,depths,scaling_model,deltat);
%
% average the results.
%
N10=mean(N10depths);
N26=mean(N26depths);
