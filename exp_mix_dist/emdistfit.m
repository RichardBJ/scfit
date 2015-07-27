function [ taus, areas, loglik, N, histvals, bins, pdfx, f, g, ...
    loghistvals, logbins, dlog, pdfxlog, flog, glog, exitflag] = emdistfit( data, init, iniw, varargin)
%EMDISTFIT Maximum likelihood estimates of the components in a mixture of
%exponential distributions
%   data - column vector of data to fit
%   init - initial guesses for the exponential components, the length of
%       init determines how many components will be estimated
%   iniw - intial guesses for the weights of each individual exponential
%       component; these should sum to 1
%   min (optional) - fit data above min
%   max (optional) - fit data below max
%   nbins (optional) - number of bins for histogram
%   plotnormal (optional) - plot a histogram of the values (rather than,
%       say their logarithm

%size the arrays so that data is a column vector and init and iniw are row
%vectors
[mx,~]=size(data);
if mx==1
    data=data';
end
[~,nt]=size(init);
if nt==1
    init=init';
    iniw=iniw';
end

tmax = max(data);
tmin = min(data);
nbins = 50;
switch nargin
    case 4
        tmin = varargin{1};
    case 5
        tmin = varargin{1};
        tmax = varargin{2};
    case 6
        tmin = varargin{1};
        tmax = varargin{2};
        nbins = varargin{3};
end

if nargin==7 && ~isempty(varargin{4})
    PlotNormalHist = varargin{4};
else
    PlotNormalHist = false;
end

datatofit = data(data>=tmin & data<=tmax);

x0 = [init, iniw];
lb = [1e-3*ones(size(init)), 1e-4*ones(size(iniw))];
ub = [inf(size(init)), ones(size(iniw))];
%constrain the sum of the weights to equal 1
Aeq = [zeros(size(init)), ones(size(iniw))];
beq = 1;
%the interior-point algorithm allows both linear equalities and bound
%constraints
% options=optimset('Algorithm','interior-point','MaxFunEvals',5e3);
options = optimoptions('fmincon','Algorithm','interior-point',...
    'MaxFunEvals',5e3,'MaxIter',5e3,'Display','final');

[params, loglik, exitflag] = fmincon (@emlik, x0, [], [], Aeq, beq, lb, ub, [], options);

taus = params(1:length(init));
areas = params(length(init)+1:end);

[ histvals, bins, g, pdfx, f, N] = emhist(datatofit,taus,areas,...
    'nbins',nbins,'PlotNormal',true,'display','off');

if PlotNormalHist
    h=figure ('Name', 'Maximum Likelihood Dwell Times', 'NumberTitle','off');
    bar(bins, histvals);
    hold on;
    plot(pdfx,g);
end

[loghistvals, logbins, glog, pdfxlog, flog, ~, dlog] = ...
    emhist(datatofit,taus,areas,'nbins',nbins,'display','off');

h(2) = figure ('Name', 'Log Dwell Time Histogram', 'NumberTitle', 'off');
bar(logbins, sqrt(loghistvals), 'barwidth', 1, 'facecolor', [0.9 0.9 0.9]);
hold on;
plot(pdfxlog,sqrt(glog),'linewidth',2);
xlabel('Log_{10} Dwell times');
ylabel('Count (square root scale)');

    function loglik = emlik(x)
        t=x(1:length(init));
        w=x(length(init)+1:end);
        pdf = emdistpdfc (datatofit, t, w);
%         pdf = emdistpdf (datatofit, t, w);
%         pdf = emdistpdf(datatofit,t,w)./(1-emdistcdf(tmin,t,w));
        %minimizing -1 * log-likelihood is equivalent to maximizing
        %log-likelihood, so multiply by -1
        loglik = -1*sum(log(pdf));
    end
        
end

