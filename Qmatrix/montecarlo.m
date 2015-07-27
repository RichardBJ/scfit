function [ states, iniState ] = montecarlo( q, c, time, wfm, nReceptors, iniState )
%MONTECARLO Run a Monte-Carlo simulation using an ion channel activation
%   mechanism described by the Q matrix, q
%   MONTECARLO creates a probability transition matrix = expm(q*adinterval)
%   and rolls the die at every point in the time vector to determine what
%   state to go to
% INPUTS
%     q - Q matrix; Q(i,j) is the rate constant to go from state i to j
%         while the diagonal elements are such that the sum of each row
%         is zero
%     c - matrix representing concentration dependence of each
%         transition; a value of 1 means the transition is
%         concentration-dependent while a value of 0 means it is not
%     time - an equally-spaced vector of time points to simulate
%     wfm (optional) - the concentration of ligand at each point in the vector
%         time. If wfm is omitted, the concentration is assumed to be 0.
%         Thus, c should probably be a zero matrix.
%     nReceptors (optional) - number of receptors to simulate
%     iniState (optional) - starting state for the Monte-Carlo simulation
% OUTPUTS
%    states
%    iniState - useful in the iniState was not given

constWfm=false;
if nargin<5 || isempty(nReceptors)
    nReceptors=1;
end
if nargin<4 || isempty(wfm)
    wfm=0;
    constWfm = true;
end
numStates = size(q,1);
conc=ones(numStates);
interval = time(2)-time(1);
states=zeros(length(time),nReceptors);
koOnes = ones(numStates+1,1);

[C, ia, ic] = unique(wfm,'R2012a');
prob = zeros(numStates+1, numStates, length(C));

for ii=1:length(ia)
    conc(c~=0) = C(ii);
    Q = conc.*q;
    Q = Q+diag(-sum(Q,2));
    A = expm(Q*interval);
    prob(2:end,:,ii) = cumsum(A,2)';
end

if nargin<6
    conc(c~=0) = wfm(1);
    Q = conc.*q;
    Q = Q-diag(sum(Q,2));
    peq = eqoccupy(Q);
    iniState = find(rand()<cumsum(peq),1);
end

p = prob(:,iniState,ic(1))*ones(1,nReceptors);
r = koOnes*rand(1,nReceptors);
states(1,:) = sum(p<r);

if constWfm
    for ii=2:length(time)
        p = prob(:,states(ii-1,:),ic(1));
        r = koOnes*rand(1,nReceptors);
        states(ii,:) = sum(p<r);
    end
else
    for ii=2:length(time)
        p = prob(:,states(ii-1,:),ic(ii));
        r = koOnes*rand(1,nReceptors);
        states(ii,:) = sum(p<r);
    end
end

end

