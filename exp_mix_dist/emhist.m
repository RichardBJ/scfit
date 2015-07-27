function [varargout] = emhist( data, t, w, varargin )
%EMHIST Plot histogram of data from exponential mixture distribution with
%the pdf overlaid
%   Input
%       data
%       t - taus of the exponential components in the mixture
%       w - weights of each component pdf = sum(w(ii)*f(ii))
%   Optional Name,Value pairs
%       'nbins' - scalar, number of bins in histogram
%       'display' - 'on' (default) or 'off'
%       'PlotNormal' - logical true or false (default) for plotting of the
%           data itself instead of the base 10 logarithm of the data
%       'npts' - number of points to calculate the pdf at, default = 1000
%       'pdfmin' - point at which to start plotting the pdf, default =  min(data)
%       'pdfmax' - point at which to stop plotting the pdf, default = max(data)

p = inputParser;
addParameter(p,'nbins',50);
addParameter(p,'PlotNormal',false);
addParameter(p,'npts',1e3);
addParameter(p,'pdfmin',min(data));
addParameter(p,'pdfmax',max(data));
addParameter(p,'display','on',@ischar);
parse(p,varargin{:});
nbins = p.Results.nbins;
PlotNormal = p.Results.PlotNormal;
npts = p.Results.npts;
pdfmin = p.Results.pdfmin;
pdfmax = p.Results.pdfmax;
Display = strcmpi(p.Results.display,'on');

[mt,nt] = size(t);
[mw,nw] = size(w);
[md,nd] = size(data);

if md==1
    data = data';
end
if mt~=1
    t = t';
end
if mw~=1
    w = w';
end

tmax = max(data);
tmin = min(data);
%N - estimated total number of events
N = numel(data) / (emdistcdf(tmax,t,w) - emdistcdf(tmin,t,w));

if PlotNormal
    [counts, binctr] = hist(data, nbins);
    % binwidth should be equal
    binwidth = mean(diff(binctr));
    %g(t) - scaled pdf
    
%     x = linspace(binctr(1),binctr(end),npts);
    x = linspace(pdfmin,pdfmax,npts);
    f=emdistpdf(x,t,w);
    g = N*binwidth*f;
    
    if Display
        hnd = figure('Name', 'Maximum Likelihood Dwell Times', 'NumberTitle','off');
        bar(binctr, counts, 'barwidth',1,'facecolor',[0.9,0.9,0.9]);
        hold on;
        plot(x,g,'linewidth',1.5);
        xlabel('Dwell Time');
        ylabel('Count');
    else
        hnd = [];
    end
else
    [counts, binctr] = hist(log10(data),nbins);
    % binwidth should be equal
    binwidth = mean(diff(binctr));

%     x = linspace(binctr(1),binctr(end),npts);
    pdfmin = log10(pdfmin);
    pdfmax = log10(pdfmax);
    x = linspace(pdfmin,pdfmax,npts);
    f = emdistpdflog10(x,t,w);
    g = N*binwidth*f;

    if Display
        hnd = figure ('Name', 'Log Dwell Time Histogram', 'NumberTitle', 'off');
        bar(binctr, sqrt(counts), 'barwidth',1,'facecolor',[0.9,0.9,0.9]);
        hold on;
        plot(x,sqrt(g),'linewidth',1.5);
        xlabel('Log_{10} Dwell Time');
        ylabel('Count (square root scale)');
    else
        hnd = [];
    end
end

if nargout > 1
    varargout = {counts,binctr,g,x,f,N,binwidth,hnd};
else
    varargout = {};
end

end

