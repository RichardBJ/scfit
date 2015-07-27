function [ cd, cs ] = concatdwells( dwells, states, varargin )
%CONCATDWELLS Concatenate contiguous open and closed durations
%   INPUT
%       dwells
%       states
%       tol (optional) - pA for real difference.  If ommitted, CONCATDWELLS
%           will combine any two dwells with different amplitudes
%       zeroAmp (optional) - if abs(amplitude) < zeroAmp, the dwell is
%           considered to be a shutting. Default is 0.

if nargin < 3 || isempty(varargin{1})
%     tol = 1e-12;
    tol = inf;
else
    tol = varargin{1};
end

if nargin < 4
    zeroAmp = 0;
else
    zeroAmp = varargin{2};
end

ii=1;
nstates = length(states);
while ii<nstates
    % skip the shut states and handle those at the end
    if abs(states(ii)) <= zeroAmp
        states(ii)=0;
        ii=ii+1;
        continue;
    end
    temp = dwells(ii);
    jj = ii+1;    
    while jj<=nstates && (abs(states(jj)-states(ii))<tol)
        if abs(states(jj)) <= zeroAmp
            states(jj) = 0;
            jj=jj+1;
            break;
        end
        temp = temp + dwells(jj);
        dwells(jj) = inf;
        states(jj) = inf;
        jj=jj+1;
    end
    dwells(ii) = temp;
    ii=jj;
end
dwells(isinf(dwells)) = [];
states(isinf(states)) = [];

idx = find(diff(states)==0);
for ii=length(idx):-1:1
    ind = idx(ii)+1;
    dwells(idx(ii)) = dwells(idx(ii)) + dwells(ind);
    dwells(ind)=inf;
    states(ind)=inf;
end
dwells(isinf(dwells)) = [];
states(isinf(states)) = [];

cd = dwells;
cs = states;

end

