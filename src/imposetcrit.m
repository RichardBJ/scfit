function [ bursts, burst_states, stats, nbursts ] = imposetcrit( dwells, states, tcrit, varargin )
%IMPOSERES Imposes resolutions for open and shut durations in an
%   idealized channel recording
%   INPUT:
%       dwells - the duration of sojourns in the states given by the input
%           array states
%       states - identity of the state of each dwell time
%       res - resolution (in milliseconds) of durations in each state
%
%   OUTPUT:
%       dwells - resolved durations
%       states - states corresponding to resolved durations
%       stats - indices of unresolved durations from original list of dwell
%           times

tol = 1e-12;

if ~isnumeric(dwells)
    error('Input dwells is in the wrong format.  It should be a numeric matrix of duration times');
    return
end
if ~isnumeric(states)
    error('Input states is in the wrong format.  It should be a numeric matrix of state identities');
    return
end

if nargin < 4
    zeroAmp=0;
else
    zeroAmp = varargin{1};
end

% to resolve durations less than or equal to the resolution
% idx = find (abs(dwells-res)<=eps | dwells+eps < res);
% idx = find (abs(dwells-res)<=tol | dwells+tol < res);

% to resolve durations strictly less than the resolution try
% idx = find (dwells+eps < res);
idx = find (dwells-tol > tcrit & abs(states)<=zeroAmp);

if isempty(idx)
    bursts = dwells;
    burst_states = states;
    stats = [];
    nbursts = 1;
    return
end

nbursts = numel(idx)+1;
bursts = nan(length(dwells),nbursts);
burst_states = nan(length(dwells),nbursts);

first=1;
for ii=1:nbursts-1
    indices = first:idx(ii)-1;
    ndwell = length(indices);
    bursts(1:ndwell,ii) = dwells(indices);
    burst_states(1:ndwell,ii) = states(indices);
    first=idx(ii)+1;
end
indices = first:length(dwells);
ndwell = length(indices);
bursts(1:ndwell,ii+1) = dwells(indices);
burst_states(1:ndwell,ii+1) = states(indices);

bursts(all(isnan(bursts),2),:) = [];
burst_states(all(isnan(burst_states),2),:) = [];
stats = idx;

end

