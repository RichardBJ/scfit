%% HJCFIT_demo
% Simulate single channel dwell times for mechanism from Q-matrix cookbook
% chapter of single channel book

% Pack Q matrix
nstates=5;
q=zeros(nstates);
% idx = sub2ind(size(q),[2 4 1 3 2 4 1 3 5 4], [1 1 2 2 3 3 4 4 4 5]);
% q(idx) = [2/3, 15, 500, 15e3, 500, 50 3e3 4e3, 10 2e3];
idxvary = sub2ind(size(q),[4 1 3 2 4 1 3 5 4], [1 2 2 3 3 4 4 4 5]);
q(idxvary) = [15, 50, 15e3, 500, 50 3e3 4e3, 10 2e3];
q(2,1) = q(1,2)*q(2,3)*q(3,4)*q(4,1) ./ (q(1,4)*q(4,3)*q(3,2));
q=q-diag(sum(q,2));
q=q*1e-3;

%open states
A=1:2;
%shut states
F=3:5;

%% Simulate data
%number of dwells to simulate
n=50e3;
%simulate data
[dwells, states] = montecarlodwells(q,A,F,n);

%resolution (aka dead time)
td=0.1; % in ms
[rd, rs] = imposeres(dwells, states, td, td);
if rs(1)~=1
    rd(1)=[];
    rs(1)=[];
end
if rs(end)~=0
    rd(end)=[];
    rs(end)=[];
end

%% Fitting

% Set constraints and provide an initial guess
idxconstrain = sub2ind(size(q),2,1);
idxall = union(idxvary,idxconstrain);
qguess = zeros(nstates);
rts = [10000 1000 5e7 0 100 30000 1000 5e7 1000 1e7]; %rates in s^-1
indices = sub2ind(size(qguess), [1 4 1 2 2 3 3 4 4 5], [4 1 2 1 3 2 4 3 5 4]);
qguess(indices)=rts;
qguess(1,2) = qguess(1,2)*0.1e-6;
qguess(4,3)=qguess(4,3)*0.1e-6;
qguess(5,4)=qguess(5,4)*0.1e-6;
qguess = qguess*1e-3; %convert rates to ms

src = sub2ind(size(qguess),[1 2 3 4], [2 3 4 1]);
tgt = sub2ind(size(qguess),[1 4 3 2], [4 3 2 1]);
[gamma,xi] = constrainrate(qguess,idxall,'loop',src,tgt);

[rates, ll, qnew] = hjcfit(rd,qguess,A,F,td,idxall,idxvary,gamma,xi);

%% Plot results
qnew = qnew - diag(sum(qnew,2));
plotqhist (qnew,A,F,rd(rs==1),rd(rs==0),td);
