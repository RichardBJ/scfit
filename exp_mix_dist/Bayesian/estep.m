function [ output_args ] = estep( y, theta_old )
%ESTEP Computes the Expectation of the 'complete-data' log posterior
%density for a hierarchical exponential mixture model with two components
%   This is the E-step of the ECM algorithm
%   See p. 526-527 of Bayesian Data Analysis - Third Edition (2014)
%   by Gelman, Carlin, Stern, Dunson, Vehtari, and Rubin. QA279.5.G45 2014
%   ISBN 978-1-4398-4095-5
%
%   INPUT:
%       y - data vector (in this case open times). Should be a
%           two-dimensional matrix, with open times in the rows and
%           channels in the columns
%       theta_old - structure with fields corresponding to the parameters
%           the distribution: theta1, theta2, mu1, mu2, lambda.

assert(isstruct(theta_old));
checkfield = @(x) assert(isfield(theta_old, x));
checkfield('theta1');
checkfield('theta2');
checkfield('mu1');
checkfield('mu2');
checkfield('sigma1');
checkfield('sigma2');
checkfield('lambda');

assert(ismatrix(y));
[maxnumopens, numchannels] = size(y);

z = -1 * ones(size(y));
% ii = 1;
% jj = 1;
% while (jj <= numchannels)
%     while (ii <= maxnumopens && y(ii, jj) ~= -1)
%         numer = 
%     end
% end
numer = theta_old.lambda * exppdf(y, theta_old.theta2);
denom = (1 - theta_old.lambda) * exppdf(y, theta_old.theta1) + numer;
z = numer ./ denom;
z(y == -1) = -1;



end

