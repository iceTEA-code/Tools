%
% age=computeage10(pp,sp,sf,cp,scaling_model)
%
% Uses Be-10 concentration to estimate the age of a sample.  
%
function age=computeage10(pp,sp,sf,cp,scaling_model)
%
% Set the maximum possible age.  This routine can fail if you give
% it an older sample, but keeping this reasonably small helps
% performance.  
%
maxage=8160;                     % half-life * 6
%
% The target N10 is sp.N10.  For 10-Be, there's no radiogenic
% concentration, and we assume that any blank has already been
% taken out.
%
targetN10=sp.concentration10-sp.inheritance10;
%
% The age will always be between left and right.  A binary search
% will be used to reduce the width of the interval.
%
left=0.0;
right=maxage;
%
% Main loop, reduce the width of the interval.  We'll settle for an 
% absolute accurate of 1.0e-6.  This is much more precise than we need
% for aging, but it is helpful in doing the calibration of production
% rates.
%
timeStepForAging=0.1;
while ((right-left)/right >1.0e-8)
%
% Take the midpoint of the interval and compute the production there.
%
  midpoint=(right+left)/2;
  [testN10,testN26]=predN1026(pp,sp,sf,cp,midpoint,scaling_model,timeStepForAging);
  if (testN10 < targetN10)
    %
    % it's older then midpoint.
    %
    left=midpoint;
  else
    %
    % it's younger then midpoint.
    %
    right=midpoint;
  end
end
age=midpoint;
