function [ openvals, shutvals, xopen, xshut ] = plothist( dwells, states, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

p = inputParser;
addRequired(p,'dwells',@isnumeric);
addRequired(p,'states',@isnumeric);
validationFxn = @(x) numel(x)==2 && all(ishghandle(x));
addOptional(p,'hFig',[],validationFxn);
parse(p,dwells,states,varargin{:});

opens = log10(dwells(states~=0));
shuts = log10(dwells(states==0));

% xopen = min(opens):0.05:max(opens);
first = floor(10*min(opens))./10;
last = ceil(10*max(opens))./10;
xopen = first:0.1:last,0.1;
xshut = min(shuts):0.1:max(shuts);

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

