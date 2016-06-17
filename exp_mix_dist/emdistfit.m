function [ taushat, weightshat, loglik, exitflag ] = emdistfit(data, taus, weights, varargin)
%EMDISTFIT Maximum likelihood estimates of the components in a mixture of
%exponential distributions
%   data - column vector of data to fit
%   taus - initial guesses for the exponential components, the length of
%       taus determines how many components will be estimated
%   weights - intial guesses for the weights of each individual exponential
%       component; these should sum to 1
%   Parameters:
%       'min' - fit data above min
%       'max' - fit data below max
%       'taulower' - lower bound on taus [default: 0.5 * min(data)]
%       'tauupper' - upper bound on taus [default: 5 * max(data)]
%       'weightlower' - lower bound on weights [default: 1e-4]
%       'conditional' - true/false to use conditional distribution (i.e.
%           p(x) given that x was in the observed range of value)
%           [default: true]

% size the arrays so that data is a column vector and
% init and iniw are row vectors
[mx,~]=size(data);
if mx==1
    data=data';
end
[~,nt]=size(taus);
if nt==1
    taus=taus';
    weights=weights';
end

% Parse inputs
% -------------------------------------------------------------------------
p = inputParser;
addRequired(p, 'data', @(x) isnumeric(x));
addRequired(p, 'taus', @isnumeric);
addRequired(p, 'weights', @isnumeric);
addParameter(p, 'min', min(data), @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'max', max(data), @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'taulower', 0.5 * min(data), @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'tauupper', 5 * max(data), @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'weightlower', 1e-4, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'conditional', true, @(x) islogical(x) && isscalar(x));

parse(p, data, taus, weights, varargin{:});
tmin = p.Results.min;
tmax = p.Results.max;
taulower = p.Results.taulower * ones(size(taus));
tauupper = p.Results.tauupper * ones(size(taus));
weightlower = p.Results.weightlower * ones(size(weights));
weightupper = ones(size(weights));
conditionalfit = p.Results.conditional;

% Set up fitting
% -------------------------------------------------------------------------
datatofit = data(data>=tmin & data<=tmax);

x0 = [taus, weights];

% ChanneLab set the lower bound for the taus to be half of the minimum
% value that is being fit, but I don't think it has a hard limit (i.e. not
% user-defined) on the areas/amplitudes
lb = [taulower , weightlower];

% ChanneLab also places a hard upper-limit (i.e. written into the code and
% not entered by the user) on the taus of five times the maximum value
% being fit.
ub = [tauupper, weightupper];

%constrain the sum of the weights to equal 1
Aeq = [zeros(size(taus)), ones(size(taus))];
beq = 1;

% Contrain the mean of the exponential components to be in increasing order
% If there are n components: tau1, tau2, ..., tau_n, weight1, weight2, ..., weight_n
% Then this constraint is specified as 
% [1 -1  0  0 ... 0 0 ... 0;      [0;
%  0  1 -1  0 ... 0 0 ... 0;  <=   0;
%  0  0  1 -1 ... 0 0 ... 0]       0]
ntaus = length(taus);
nweights = length(weights);
assert(ntaus == nweights);
A = zeros(ntaus - 1, ntaus + nweights);
b = zeros(ntaus - 1, 1);
A(sub2ind(size(A), 1:(ntaus - 1), 1:(ntaus - 1))) = 1;
A(sub2ind(size(A), 1:(ntaus - 1), 2:ntaus)) = -1;

% the interior-point algorithm allows both linear equalities and
% bound constraints
% options=optimset('Algorithm','interior-point','MaxFunEvals',5e3);
options = optimoptions('fmincon','Algorithm','interior-point',...
    'MaxFunEvals',5e3,'MaxIter',5e3,'Display','final');

% -------------------------------------------------------------------------
% Find the maximum likelihood parameters for the
% exponential mixture distribution
% -------------------------------------------------------------------------
[params, loglik, exitflag] = fmincon (@emlik, x0, A, b, Aeq, beq, lb, ...
                                      ub, [], options);
taushat = params(1:ntaus);
weightshat = params((ntaus+1):end);                                  

% -------------------------------------------------------------------------
% likelihood function for exponential mixture distribution
% -------------------------------------------------------------------------
    function loglik = emlik(x)
        t=x(1:ntaus);
        w=x((ntaus+1):end);
        
        if conditionalfit
            pdf = emdistpdfc(datatofit, t, w);
        else
            pdf = emdistpdf(datatofit, t, w);
        end
%         pdf = emdistpdf(datatofit,t,w)./(1-emdistcdf(tmin,t,w));

        % minimizing -1 * log-likelihood is equivalent to maximizing
        % log-likelihood, so multiply by -1
        loglik = -1*sum(log(pdf));
    end
        
end

