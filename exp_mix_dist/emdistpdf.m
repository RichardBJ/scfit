function pdf = emdistpdf (x, t, a)
%EMDISTPDF Probability density function for an exponential mixture
%distribution
%   x - point at which to return the pdf
%   t - taus of the exponential mixture
%   a - areas of the mixture components

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

pdf = exp(-x*(1./t)) * (a./t)';

if flip
    pdf=pdf';
end

end

