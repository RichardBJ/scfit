%% HJCFIT_demo
% Simulate single channel dwell times for mechanism from Q-matrix cookbook
% chapter of single channel book

% Pack Q matrix
%rates in s^-1 from Vance et al (2012) J Physiol for GluN1-1a/GluN2D
bp = 0.25; % in uM^-1*s^-1
bm = 0.21;
k1p = 130;
k1m = 990;
k2p = 880;
k2m = 11000;
k3p = 6800;
k3m = 6200;
k4p = 7400;
k4m = 3600;
k5p = 4.6;
k5m = 3.4;
k6p = 140;
k6m = 480;
k7p = 410;
k7m = 9100;
k8p = 820;
k8m = 12000;
k9p = 5600;
k9m = 4500;

q=zeros(14);
q(1,2)=2*bp;
q(1,8)=k5p;
q(2,1)=bm;
q(2,3)=bp;
q(2,9)=k5p;
q(3,2)=2*bm;
q(3,4)=k1p;
q(3,10)=k5p;
q(4,3)=k1m;
q(4,5)=k2p;
q(4,11)=k5p;
q(5,4)=k2m;
q(5,6)=k3p;
q(5,12)=k5p;
q(6,5)=k3m;
q(6,7)=k4p;
q(7,6)=k4m;
q(8,1)=k5m;
q(8,9)=2*bp;
q(9,2)=k5m;
q(9,8)=bm;
q(9,10)=bp;
q(10,3)=k5m;
q(10,9)=2*bm;
q(10,11)=k6p;
q(11,4)=k5m;
q(11,10)=k6m;
q(11,12)=k7p;
q(12,5)=k5m;
q(12,11)=k7m;
q(12,13)=k8p;
q(13,12)=k8m;
q(13,14)=k9p;
q(14,13)=k9m;

c = zeros(size(q));
c(1,2)=1;c(2,3)=1;c(8,9)=1;c(9,10)=1;

A = [6 7 13 14];
F = setdiff(1:14,A);

%% Set the constraints for the model
% We only need to vary 20 of the 34 rate constants in Q. If we use the
% logarithm of the rate constants then all of the constraints become linear
% and can easily be manipulated
idxall=find(q);
idxvary = sub2ind(size(q), [1 1 2 3 4 4 5 5 8 10 11 11 12 12 6 6 7 13 13 14], [2 8 1 4 3 5 4 6 1 11 10 12 11 13 5 7 6 12 14 13])';

src = sub2ind(size(q), [1 1 1 2 2 2], [2 2 2 1 1 1]);
tgt = sub2ind(size(q), [2 8 9 3 9 10], [3 9 10 2 8 9]);
[gamma,xi] = constrainrate(q,idxall,'constrain',src,tgt,[0.5 1 0.5 2 1 2]);
src = sub2ind(size(q), [1 1 1 1 8 8 8 8], [8 8 8 8 1 1 1 1]);
tgt = sub2ind(size(q), [2 3 4 5 9 10 11 12], [9 10 11 12 2 3 4 5]);
[tmpg,tmpxi] = constrainrate(q,idxall,'constrain',src,tgt,[1 1 1 1 1 1 1 1]);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];

nRates = length(idxall);
idx2 = interp1(idxall, 1:nRates, idxvary);
[~, idx3] = setdiff (idxall, idxvary); 

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
