function [openbins, openhist, xopen, openpdfhjc, openpdfexact, ...
    shutbins, shuthist, xshut, shutpdfhjc, shutpdfexact ] ...
    = plotqhist (q, A, F, opens, shuts, td, varargin)
%PLOTQHIST Plots the HJC distribution of open and shut times for the
%reaction mechanism give by q and the resolution imposed, td
%   INPUT
%       q - Q matrix
%       A - vector identifying the open states
%       F - vector identifying the shut states
%       opens - vector of open times
%       shuts - vector of shut times
%       td - resolution imposed on the data
%       'nbins'(optional) - number of bins for the histogram. Can also be a
%           vector specifying the center of each bin
%       'OpenInitialVal'(optional) open_ini(optional) - log10 of time at which to start plotting the open time
%           distribution
%       'ShutInitialVal'(optional) shut_ini(optional) - log10 of the time at which to start the shut time
%           distribution
%       'NormMode' normalize_mode(optional) - 'pdf' (default) or 'hist'
%           'pdf' will scale the pdf to the counts in the histogram
%           'hist' will normalize the histogram counts so the area is one;
%             note that the exact pdf will also need to be scaled so
%             the area from the dead time to the maximum dwell time is one.
%       'Plot'(optional) - 'on' (default) or 'off'
%       'BarColor' (optional) - color of histogram bars

p = inputParser;
addParameter(p,'nbins',50,@isnumeric);
addParameter(p,'OpenInitialVal',min(log10(opens)),@isnumeric);
addParameter(p,'ShutInitialVal',min(log10(shuts)),@isnumeric);
addParameter(p,'NormMode','pdf',@(x) ischar(x)&&(strcmpi(x,'pdf')||strcmpi(x,'hist')));
addParameter(p,'Plot','on',@(x) ischar(x)&&(strcmpi(x,'on')||strcmpi(x,'off')));
addParameter(p, 'BarColor', [0.230 0.299 0.754]);
parse(p,varargin{:});
nbins = p.Results.nbins;
xopenini = p.Results.OpenInitialVal;
xshutini = p.Results.ShutInitialVal;
mode = p.Results.NormMode;
PLOTFLAG = p.Results.Plot;
barclr = p.Results.BarColor;

[openhist,openbins] = hist(log10(opens),nbins);
dopen = mean(diff(openbins));
nopen = numel(opens);
xopen = linspace(xopenini,openbins(end),1e3);
openpdfhjc = hjcdist(q,A,F,td,xopen);
[opentimes, openareas, shuttimes,shutareas] = qmatopenshutpdf(q,A,F);
openpdfexact = emdistpdflog10(xopen,opentimes,openareas);
% exact_open_pdf_area = emdistcdf(100,opentimes,openareas) - emdistcdf(td,opentimes,openareas);
exact_open_pdf_area = emdistcdf(max(opens),opentimes,openareas) - emdistcdf(td,opentimes,openareas);
Nopen = nopen ./ exact_open_pdf_area;

if strcmpi(mode,'pdf')
    openpdfhjc = nopen*dopen*openpdfhjc;
    openpdfexact = Nopen*dopen*openpdfexact;
else
    openhist = openhist ./ (nopen*dopen);
    openpdfexact = openpdfexact ./ exact_open_pdf_area;
end

if strcmpi(PLOTFLAG,'on')
    figure ('Name','Open Time Histogram', 'NumberTitle', 'off');
    bar(openbins,sqrt(openhist), 1, 'facecolor', barclr);
    hold on;
    plot(xopen,sqrt(openpdfhjc),'k','LineWidth',2);
    plot(xopen,sqrt(openpdfexact),':k','LineWidth',2);
    xlabel ('Log Open Time (ms)');
    ylabel ('Count (square root scale)');
    hold off;
end

nshut=length(shuts);
[shuthist,shutbins]=hist(log10(shuts),nbins);
dshut=mean(diff(shutbins));
xshut=linspace(xshutini,shutbins(end),1001);
% exact_shut_pdf_area = emdistcdf(100000,shuttimes,shutareas)-emdistcdf(td,shuttimes,shutareas);
exact_shut_pdf_area = emdistcdf(max(shuts),shuttimes,shutareas)-emdistcdf(td,shuttimes,shutareas);
Nshut = nshut ./ exact_shut_pdf_area;
shutpdfhjc = hjcdist(q,F,A,td,xshut);
shutpdfexact = emdistpdflog10(xshut,shuttimes,shutareas);

if strcmpi(mode,'pdf')
    shutpdfhjc = nshut*dshut*shutpdfhjc;
    shutpdfexact = Nshut*dshut*shutpdfexact;
else
    shuthist = shuthist ./ (nshut*dshut);
    shutpdfexact = shutpdfexact ./ exact_shut_pdf_area;
end

if strcmpi(PLOTFLAG,'on')
    figure ('Name', 'Shut Time Histogram', 'NumberTitle', 'off');
    bar(shutbins,sqrt(shuthist), 1, 'facecolor', barclr);
    hold on;
    plot(xshut,sqrt(shutpdfhjc),'k','LineWidth',2);
    plot(xshut,sqrt(shutpdfexact),':k','LineWidth',2);
    xlabel('Log Shut Time (ms)');
    ylabel('Count (square root scale)');
    hold off;
end

end
