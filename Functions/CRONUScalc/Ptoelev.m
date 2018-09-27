%
%  elev=Ptoelev(lat,long,P)
%
%  Converts a pressure into a corresponding elevation using the
%  ERA40 atmosphere.
%
function elev=Ptoelev(lat,long,P)
left=-500;
right=7000;
while (right-left > 1)
  mid=(left+right)/2;
  Ptest=ERA40atm(lat,long,mid);
  if (Ptest > P)
    left=mid;
  else
    right=mid;
  end
end
elev=(left+right)/2;
