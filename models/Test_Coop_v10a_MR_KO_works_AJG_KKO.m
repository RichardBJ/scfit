%% NMDA Receptor Model - Test_Coop_v10a_MR_KO_works_AJG
%
% Build a model of NMDA receptor activation
%
alpha3=	2000	;
beta3=	1000	;
alpha1=	500	;
beta1=	10000	;
alpha2=	5000	;
beta2=	1000	;
kfp=	500	;
kfm=	1500	;
ksp=	5	;
ksm=	30	;
kftp=	1000	;
kftm=	2000	;
kstp=	50	;
kstm=	200	;
%
q=zeros(13);
%
q(1,4) = alpha1;
q(4,1) = beta1;
q(2,5) = alpha2;
q(5,2) = beta2;
q(13,6) = alpha3;
q(6,13) = beta3;
q(13,1) = kftp;
q(1,13) = 2*kftm;
q(2,1) = kstp;
q(1,2) = kstm;
q(4,5) = 2*kstm;
q(5,4) = kstp;
q(4,6) = 2*kftm;
q(6,4) = kftp;
q(5,7) = kstm;
q(7,5) = 2*kstp;
q(5,8) = kfm;
q(8,5) = kfp;
q(6,8) = ksm;
q(8,6) = ksp;
q(6,9) = kftm;
q(9,6) = 2*kftp;
q(7,10) = 2*kfm;
q(10,7) = kfp;
q(8,10) = kstm;
q(10,8) = kstp;
q(8,11) = kftm;
q(11,8) = kftp;
q(9,11) = 2*ksm;
q(11,9) = ksp;
q(10,12) = kfm;
q(12,10) = 2*kfp;
q(11,12) = ksm;
q(12,11) = 2*ksp;
q(11,3) = kfp;
q(3,11) = kfm;
q(10,3) = ksp;
q(3,10) = ksm;
q(6,3) = kstp;
q(3,6) = kstm;
q(5,3) = kftp;
q(3,5) = kftm;
%
% openstates
A = [1,2,13];

% shutstates
F = 3:12;

% scale rates to ms^-1
q = q*1e-3;

% q = q-diag(sum(q,2));

%% Constraints

idxall = find(q);

% First define the rates to be constrained and then set the variable rates
% at the end

% declare the idxconstrain variable. It's empty now, but we'll add rates to it as we go.
idxconstrain = []; 
gamma = [];
xi = [];

% Physical constraints

% Constrain rates k+s
src = sub2ind(size(q), 11*ones(1,3), 9*ones(1,3));
tgt = sub2ind(size(q), [10 8 12], [3 6 11]);
[tmpg,tmpxi] = constrainrate (q, idxall, 'constrain', src, tgt, [1 1 2]);
gamma = [gamma; tmpg];
xi = [xi; tmpxi];
% need to transpose the tgt vector to add it to idxconstrain
idxconstrain = [idxconstrain; tgt'];

% Constrain rates k-s
src = sub2ind(size(q), 11*ones(1,1), 12*ones(1,1));
tgt = sub2ind(size(q), [9], [11]);
[tmpg,tmpxi] = constrainrate (q, idxall, 'constrain', src, tgt, [2]);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];
% need to transpose the tgt vector to add it to idxconstrain
idxconstrain = [idxconstrain; tgt'];

% Constrain rates k+s'
src = sub2ind(size(q), 10*ones(1,3), 8*ones(1,3));
tgt = sub2ind(size(q), [2 3 7], [1 6 5]);
[tmpg,tmpxi] = constrainrate (q, idxall, 'constrain', src, tgt, [1 1  2]);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];
% need to transpose the tgt vector to add it to idxconstrain
idxconstrain = [idxconstrain; tgt'];

% % Constrain rates k-s'
src = sub2ind(size(q), 5*ones(1,1), 7*ones(1,1));
tgt = sub2ind(size(q), [4], [5]);
[tmpg,tmpxi] = constrainrate (q, idxall, 'constrain', src, tgt, [2]);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];
% need to transpose the tgt vector to add it to idxconstrain
idxconstrain = [idxconstrain; tgt'];

% Constrain rates k+f
src = sub2ind(size(q), 10*ones(1,3), 7*ones(1,3));
tgt = sub2ind(size(q), [8 11 12], [5 3 10]);
[tmpg,tmpxi] = constrainrate (q, idxall, 'constrain', src, tgt, [1 1 2]);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];
% need to transpose the tgt vector to add it to idxconstrain
idxconstrain = [idxconstrain; tgt'];

% Constrain rates k-f
src = sub2ind(size(q), 10*ones(1,3), 12*ones(1,3));
tgt = sub2ind(size(q), [3 5 7], [11 8 10]);
[tmpg,tmpxi] = constrainrate (q, idxall, 'constrain', src, tgt, [1 1 2]);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];
% need to transpose the tgt vector to add it to idxconstrain
idxconstrain = [idxconstrain; tgt'];

% Constrain rates k+f'
src = sub2ind(size(q), 6*ones(1,4), 4*ones(1,4));
tgt = sub2ind(size(q), [3 11 9 13], [5 8 6 1]);
[tmpg,tmpxi] = constrainrate (q, idxall, 'constrain', src, tgt, [1 1 2 1]);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];
% need to transpose the tgt vector to add it to idxconstrain
idxconstrain = [idxconstrain; tgt'];

% Constrain rates k-f'
src = sub2ind(size(q), 6*ones(1,3), 9*ones(1,3));
tgt = sub2ind(size(q), [8 4 1], [11 6 13]);
[tmpg,tmpxi] = constrainrate (q, idxall, 'constrain', src, tgt, [1 2 2]);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];
% need to transpose the tgt vector to add it to idxconstrain
idxconstrain = [idxconstrain; tgt'];

% Microscopic Reversibility
[gamma,xi,idxconstrain,idxMR] = mrconstraints(q,idxall,idxconstrain,[],gamma,xi);

% create idxvary from the rates that are not constrained
idxvary = setdiff(idxall,idxconstrain);

% reorder idxall and the columns of gamma so that the constrained rates are
% first
% [~,iAllToTheta] = ismember([idxconstrain;idxvary],idxall);
% gamma = gamma(:,iAllToTheta);
% idxall = [idxconstrain;idxvary];
% clear iAllToTheta;

% Use reduced row echelon form to select a linearly independent set of
% columns from gamma -> these will be a basis for the range of gamma
% Then rearrange gamma so that the independent columns are the first
% columns
[~,jb] = rref(gamma);
jc = setdiff(1:length(gamma),jb);
idxconstrain = idxall(jb);
idxvary = idxall(jc);
idxall = idxall([jb,jc]);
gamma = gamma(:,[jb,jc]);

% % Solve the underdetermined system for the minimum norm solution (i.e. find
% % [one of the infinitely many sets of] rates satisfying the constraints)
% [Q,R] = qr(gamma');
% [~,R1] = qr(gamma',0);
% x = Q*[R1'\xi;zeros(length(Q)-length(R1),1)];
% 
% q(idxall) = 10.^x;