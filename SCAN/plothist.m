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

p = inputParser;
addRequired(p,'dwells',@isnumeric);
addRequired(p,'states',@isnumeric);
addOptional(p, 'obw', 0.05, @isnumeric);
addOptional(p, 'sbw', 0.1, @isnumeric);
validationFxn = @(x) numel(x)==2 && all(ishghandle(x));
addOptional(p,'hFig',[],validationFxn);
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
else
    figure(p.Results.hFig(1));
end
bar(xopen,sqrt(openvals),1);
xlabel('Log [open times (ms)]');
ylabel('Count (square root scale)');

if isempty(p.Results.hFig)
    figure;
else
    figure(p.Results.hFig(2));
end
bar(xshut,sqrt(shutvals),1);
xlabel('Log [shut times (ms)]');
ylabel('Count (square root scale)');

end

