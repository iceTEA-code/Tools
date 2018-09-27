function out = Rv0(z)

% this subfunction returns the stopping rate of vertically traveling muons
% as a function of depth z at sea level and high latitude.

a = exp(-5.5e-6.*z);
b = z + 21000;
c = (z + 1000).^1.66 + 1.567e5;
dadz = -5.5e-6 .* exp(-5.5e-6.*z);
dbdz = 1;
dcdz = 1.66.*(z + 1000).^0.66;

out = -5.401e7 .* (b.*c.*dadz - a.*(c.*dbdz + b.*dcdz))./(b.^2 .* c.^2);

% full depth calculation appears in comments below
%R_1 = -5.401e7 .* (b.*c.*dadz - a.*(c.*dbdz + b.*dcdz))./(b.^2 .* c.^2);
%f = (121100./z).^2;
%g = exp(-z./121100);
%dfdz = (-2.*(121100.^2))./(z.^3);
%dgdz = -exp(-z./121100)./121100;
%R_2 = -1.82e-6.*(g.*dfdz + f.*dgdz);
%out(find(z<200000)) = R_1(find(z<200000));
%out(find(z>=200000)) = R_2(find(z>=200000));

% -------------------------------------------------------------------------
