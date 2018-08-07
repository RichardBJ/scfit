function [ pdf ] = expmixdistpdf( x, taus, areas, varargin)
%EXPMIXDISTPDF pdf of an exponential mixture distribution
%   Calculate the probability density function of an exponential mixture
%   and return each component of the mixture
%
%   Arguments:
%       x - points at which to calculate the pdf
%       taus - time components comprising the exponential mixture
%           distribution (a column vector)
%       areas - proportion of the mixture distribution corresponding to the
%           taus (a column vector)
%       log10 - a name-value pair of 'log10', val where val is true or
%       false indicating if the values in x are log10 transformed values
%
%   Return:
%       pdf = probability density function as an m-by-n matrix where m is
%       the the number of values in x and n is the number of components in
%       the distribution

assert(size(x, 2) == 1, 'x must be a column vector');
assert(size(taus, 1) == 1, 'taus must be a row vector');
assert(all(size(areas) == size(taus)), 'taus and areas must have the same size');

p = inputParser;
addParameter(p,'log10', false, @islogical);
parse(p,varargin{:});
LOGFLAG = p.Results.log10;

at_multiplier = diag(areas ./ taus);

if LOGFLAG == true
    u = 10.^x;
    pdf = exp(-u*(1./taus)) * at_multiplier;
    % pdf = u.*pdf*log(10);
    pdf = bsxfun(@times, u, pdf*log(10));
else
    pdf = exp(-x*(1./taus)) * at_multiplier;
end

end

