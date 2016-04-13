%% HJCFIT_demo
% Simulate single channel dwell times for mechanism from Q-matrix cookbook
% chapter of single channel book

% Pack Q matrix
%GluN2A Rates in uM^-1*s^-1
ksp = 230;
ksm = 178;
kfp = 3140;
kfm = 174;
kd1p = 85.1;
kd1m = 29.7;
kd2p = 230;
kd2m = 1.01;
kon = 31.6;
koff = 1010;

q=zeros(8);
q(1,2)=2*kon;
q(2,1)=koff;
q(2,3)=kon;
q(3,2)=2*koff;
q(3,4)=kfp;
q(4,3)=kfm;
q(3,5)=ksp;
q(5,3)=ksm;
q(4,6)=ksp;
q(6,4)=ksm;
q(5,6)=kfp;
q(6,5)=kfm;
q(3,7)=kd1p;
q(7,3)=kd1m;
q(3,8)=kd2p;
q(8,3)=kd2m;

%openstates
A = 6;
%shutstates
F = [1:5, 7:8];

c=zeros(8);
c(1,2)=1;
c(2,3)=1;

%% Set the constraints for the model
idxall = sub2ind(size(Q), [1 2 2 3 3 3 3 3 4 4 5 5 6 6 7 8], [2 1 3 2 4 5 7 8 3 6 3 6 4 5 3 3])';
idxconstrain = sub2ind(size(Q), [2 3 4 6 5 6], [3 2 6 4 6 5]);
idxvary = setdiff(idxall,idxconstrain);
src = sub2ind(size(Q), [1 2], [2 1]);
tgt = sub2ind(size(Q), [2 3], [3 2]);
[gamma,xi] = constrainrate (Q, idxall, 'constrain', src, tgt, [0.5 2]);
src = sub2ind(size(Q), [3 3 4 5], [4 5 3 3]);
tgt = sub2ind(size(Q), [5 4 6 6], [6 6 5 4]);
[tmpgamma, tmpxi] = constrainrate(Q,idxall,'constrain',src,tgt,[1 1 1 1]);
gamma=[gamma;tmpgamma];
xi=[xi;tmpxi];

% Following (Qin et al. 1996 Biophys J) or (Golub and Van Loan 1989)
% We can formulate the linear constraints on mu, which are represented in the
% matrix gamma, into linear combinations of unconstrained variables using
% the QR factorization of gamma
%                           gamma * theta = xi
%   where theta = (...mu(i,j)...)' and xi is a constant vector
%   assume the constraints are independent of one another, which means
%   gamma has full rank
[U,R] = qr(gamma);
R1 = R(:,idx3);
R2 = R(:,idx2);

%% Simulate data
q(c==1) = q(c==1)*1e3;
q = q-diag(sum(q,2));
q = q*1e-3;
%number of dwells to simulate
n=50e3;
%simulate data
[dwells, states] = montecarlodwells(q,A,F,n);

%resolution (aka dead time)
td=0.1; % in ms
[rd, rs] = imposeres(dwells,states,td);
if rs(1)~=1
    rd(1)=[];
    rs(1)=[];
end
if rs(end)~=0
    rd(end)=[];
    rs(end)=[];
end

%% Fitting
qguess = q.*2.^randn(size(q));
% qguess = qnew;

[rates, ll, qnew] = hjcfit(rd,qguess,A,F,td,idxall,idxvary,gamma,xi);

%% Plot results
qnew = qnew - diag(sum(qnew,2));
plotqhist (qnew,A,F,rd(rs==1),rd(rs==0),td);
