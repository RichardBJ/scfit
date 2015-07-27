function cdf = emdistcdf(x, t, a)
%EMDISTCDF Returns the cumulative distribution function for an exponential
%mixture
%   x - point at which to return the cdf
%   t - taus of the exponential mixture
%   a - areas of the exponential mixture

flip=false;
[mx,~]=size(x);
if mx==1
    x=x';
    flip=true;
end
[~,nt]=size(t);
if nt==1
    t=t';
    a=a';
end

cdf = 1 - exp(-x*(1./t)) * a';

if flip
    cdf = cdf';
end

end

