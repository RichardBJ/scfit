function pdf = emdistpdfc (x, t, w)
%EMDISTPDFC Conditional pdf of a mixture of exponential distributions
%    on only observing values between min(x) and max(x)
%   x - column vector of values at which to evalutate the pdf
%   t - row vector of the taus of each exponential component
%   w - row vector of the amplitudes of each exponential component

tmax = max(x);
tmin = min(x);
f = exp(-x*(1./t)) * (w./t)';
p_first_last = emdistcdf(tmax,t,w) - emdistcdf(tmin,t,w);
pdf = f ./ p_first_last;

end