function [ dwells, states, stats ] = imposeres( dwells, states, openres, shutres )
%IMPOSERES Imposes resolutions for open and shut durations in an
%   idealized channel recording
%   INPUT:
%       dwells  - the duration of sojourns in the states given by the input
%                 array states
%       states  - identity of the state of each dwell time
%       openres - resolution (in milliseconds) of durations in states other
%                 than 0
%       shutres - shut time resolution (in ms) of durations in state 0
%
%   OUTPUT:
%       dwells - resolved durations
%       states - states corresponding to resolved durations
%       stats  - indices of unresolved durations from original list of
%                dwell times

tol = 1e-12;

if ~isnumeric(dwells)
    error('Input dwells is in the wrong format.  It should be a numeric matrix of duration times');
    return
end
if ~isnumeric(states)
    error('Input states is in the wrong format.  It should be a numeric matrix of state identities');
    return
end
if ~isfloat(states)
    error('Input states is in the wrong format. It should be a floating poing matrix');
    return
end

% Make sure we start with a resolvable dwell
ii=1;
while (dwells(ii)+tol<openres & states(ii)~=0) || (dwells(ii)+tol<shutres & states(ii)==0)
    dwells(ii)=[];
    states(ii)=[];
%     ii=ii+1;
%     if ii>=length(dwells)
    if isempty(dwells)
        break;
    end
end

% to resolve durations less than or equal to the resolution
% idx = find (abs(dwells-res)<=eps | dwells+eps < res);
% idx = find (abs(dwells-res)<=tol | dwells+tol < res);

% to resolve durations strictly less than the resolution try
% idx = find (dwells+eps < res);
% idx = find (dwells+tol < res);
idxopen = find(dwells+tol < openres & states~=0);
idxshut = find(dwells+tol < shutres & states==0);
[m,n] = size(idxopen);
[m2,n2] = size(idxshut);
% if m>n 
if n==n2
    idx = sort([idxopen; idxshut]);
elseif m==m2
    idx = sort([idxopen, idxshut]);
end

ii=1;
tempdwell=0;

while ii<=length(idx)
   tempdwell = dwells(idx(ii));
   dwells(idx(ii))=inf;
   states(idx(ii))=inf;
   jj=ii+1;
   while (jj<=length(idx)) && (idx(jj)-idx(jj-1) == 1)
        tempdwell = tempdwell + dwells(idx(jj));
        dwells(idx(jj))=inf;
        states(idx(jj))=inf;
        jj=jj+1;
   end
   dwells(idx(ii)-1) = dwells(idx(ii)-1) + tempdwell;
   ii=jj;
end

dwells(isinf(dwells)) = [];
states(isinf(states)) = [];

stats = idx;

idx = find(diff(states)==0);
for ii=length(idx):-1:1
    ind = idx(ii)+1;
    dwells(idx(ii)) = dwells(idx(ii)) + dwells(ind);
    dwells(ind)=inf;
    states(ind)=inf;
end
dwells(isinf(dwells)) = [];
states(isinf(states)) = [];

% stats = cell(10,2);
% while ~isempty(idx)
%     stats{ii,1}=idx;
%     ii=ii+1;
%     idx = find(abs(filter(ones(ii,1)/ii,1,dwells)-res)<eps);
% end
% 
% for jj=ii-1:-1:1
%     idx = find(abs(filter(ones(jj,1)/jj,1,dwells)-res)<eps);
%     dwells(idx-jj)=dwells(idx-jj)+dwells(idx);
%     dwells(idx)=[];
%     states(idx)=[];
%     stats{jj,2}=idx;
% %     fprintf('Number of %d consecutive irresolvable durations: %d\n', jj, length(idx));
% end

end

