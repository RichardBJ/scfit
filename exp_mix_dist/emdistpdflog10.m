function pdf = emdistpdflog10(x, t, w)
%EMDISTPDFLOG10 Probability density function of an exponential mixtures
%distribution for the logarithm of the variable
%   If u is a random variable whose distribution is a mixture of exponential
% distributions, and x = log10(u), then the pdf of x = u*ln(10)*pdf(u)
%   x - logarithm of the value at which to get the exponential mixture pdf

flip=false;
[mx,~]=size(x);
if mx==1
    x=x';
    flip=true;
end
[~,nt]=size(t);
if nt==1
    t=t';
    w=w';
end

% y = zeros(size(x));
u = 10.^x;
% for ii=1:length(t)
%     y = y + w(ii)*exppdf(u,t(ii)).*u*log(10);
% end

pdf = u.*emdistpdf(u,t,w)*log(10);

if flip
    pdf = pdf';
end

end

