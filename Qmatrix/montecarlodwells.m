function [dwells, states] = montecarlodwells (q,A,F,n,iniState)
%MONTECARLODWELLS Monte Carlo simulation of a continuous time Markov
%process given by the transition rate matrix q
%
%INPUT
%   q - Q matrix -- cannot vary with time
%   A - vector of indices of the open (or up) state
%   F - vector of indices of the closed (or down) state
%   n - number of transitions
%   iniState (optional) - initial state of the system
%
%OUTPUT
%   dwells - dwell times
%   states - states corresponding to the dwell times
%       0 = shut (or down)
%       1 = open (or up)

dwells = zeros(n,1);
states = nan(n,1);

if nargin < 5
    p0 = eqoccupy(q);
    state = find(rand()<=cumsum(p0),1);
else
    state = iniState;
end

for ii=1:n
    if any(state==A)
        while any(state==A)
            time = exprnd(-1./q(state,state),1);
            dwells(ii) = dwells(ii) + time;
            idx = find(state~=1:length(q));
            pt = q(state,idx)./sum(q(state,idx));
            state = idx(find(rand()<=cumsum(pt),1));
        end
        states(ii) = 1;
    else
        while any(state==F)         
            time = exprnd(-1./q(state,state),1);
            dwells(ii) = dwells(ii) + time;
            idx = find(state~=1:length(q));
            pt = q(state,idx)./sum(q(state,idx));
            state = idx(find(rand()<=cumsum(pt),1));
        end
        states(ii) = 0;
    end
end

end
