function [ openvals, shutvals, xopen, xshut ] = plothist( dwells, states, varargin )
%PLOTHIST Summary of this function goes here
%   INPUT
%       dwells - a numerical vector (1-dimension) of the dwell times in
%           milliseconds
%       states - a numerical vector of same length as dwells giving the
%           state of each dwell (should only have zeroes for closed and
%           non-zeroes for open)
%       obw (OPTIONAL) - open bin width in millisecs [default: 0.05]
%       sbw (OPTIONAL) - shut bin width in millisecs [default: 0.1]
%       hFig (OPTIONAL) - handle for plotting [default: [] ]
%       axes (OPTIONAL) - logical indicating whether hFig contains axis
%       handles or figure handles

p = inputParser;
addRequired(p,'dwells',@isnumeric);
addRequired(p,'states',@isnumeric);
addOptional(p, 'obw', 0.05, @isnumeric);
addOptional(p, 'sbw', 0.1, @isnumeric);
validationFxn = @(x) numel(x)==2 && all(ishghandle(x));
addOptional(p,'hFig',[],validationFxn);
addOptional(p, 'axes', false, @islogical);
parse(p,dwells,states,varargin{:});

obw = p.Results.obw;
sbw = p.Results.sbw;

opens = log10(dwells(states~=0));
shuts = log10(dwells(states==0));

% xopen = min(opens):0.05:max(opens);
first = floor(10*min(opens))./10;
last = ceil(10*max(opens))./10;
% xopen = first:0.1:last;
xopen = first:obw:last;
% xshut = min(shuts):0.1:max(shuts);
xshut = min(shuts):sbw:max(shuts);

openvals = hist(opens,xopen);
shutvals = hist(shuts,xshut);

if isempty(p.Results.hFig)
    figure;
elseif ~p.Results.axes
    figure(p.Results.hFig(1));
end
if p.Results.axes && ~isempty(p.Results.hFig)
    bar(p.Results.hFig(1), xopen, sqrt(openvals), 1);
    xlabel(p.Results.hFig(1), 'Log [open times (ms)]');
    ylabel(p.Results.hFig(1), 'Count (square root scale)');
else
    bar(xopen,sqrt(openvals),1);
    xlabel('Log [open times (ms)]');
    ylabel('Count (square root scale)');
end


if isempty(p.Results.hFig)
    figure;
elseif ~p.Results.axes
    figure(p.Results.hFig(2));
end
if ~p.Results.axes
    bar(xshut,sqrt(shutvals),1);
    xlabel('Log [shut times (ms)]');
    ylabel('Count (square root scale)');
elseif ~isempty(p.Results.hFig)
    bar(p.Results.hFig(2), xshut, sqrt(shutvals), 1);
    xlabel(p.Results.hFig(2), 'Log [shut time (ms)]');
    ylabel(p.Results.hFig(2), 'Count (square root scale)');
end

end

