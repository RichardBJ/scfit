%% NMDA Receptor Model - Test_Coop_v10a_MR_KO_works_AJG
%
% Build a model of NMDA receptor activation

kA = 500;
kB = 5;
kAm = 1500;
kBm = 30;
kA_withB = 1000;
kB_withA = 50;
kAm_withB = 2000;
kBm_withA = 200;
alpha1 = 500;
alpha2 = 5000;
beta1 = 10000;
beta2 = 1000;

q=zeros(18);

q(1,2) = kB;
q(1,3) = kB;
q(1,5) = kA;
q(1,9) = kA;
q(2,1) = kBm;
q(2,4) = kB;
q(2,6) = kA_withB;
q(2,10) = kA;
q(3,1) = kBm;
q(3,4) = kB;
q(3,7) = kA;
q(3,11) = kA_withB;
q(4,2) = kBm;
q(4,3) = kBm;
q(4,8) = kA_withB;
q(4,12) = kA_withB;
q(5,1) = kAm;
q(5,6) = kB_withA;
q(5,7) = kB;
q(5,13) = kA;
q(6,2) = kAm_withB;
q(6,5) = kBm_withA;
q(6,8) = kB;
q(6,14) = kA;
q(7,3) = kAm;
q(7,5) = kBm;
q(7,8) = kB_withA;
q(7,15) = kA_withB;
q(8,4) = kAm_withB;
q(8,6) = kBm;
q(8,7) = kBm_withA;
q(8,16) = kA_withB;
q(9,1) = kAm;
q(9,10) = kB;
q(9,11) = kB_withA;
q(9,13) = kA;
q(10,2) = kAm;
q(10,9) = kBm;
q(10,12) = kB_withA;
q(10,14) = kA_withB;
q(11,3) = kAm_withB;
q(11,9) = kBm_withA;
q(11,12) = kB;
q(11,15) = kA;
q(12,4) = kAm_withB;
q(12,10) = kBm_withA;
q(12,11) = kBm;
q(12,16) = kA_withB;
q(13,5) = kAm;
q(13,9) = kAm;
q(13,14) = kB_withA;
q(13,15) = kB_withA;
q(14,6) = kAm;
q(14,10) = kAm_withB;
q(14,13) = kBm_withA;
q(14,16) = kB_withA;
q(15,7) = kAm_withB;
q(15,11) = kAm;
q(15,13) = kBm_withA;
q(15,16) = kB_withA;
q(16,8) = kAm_withB;
q(16,12) = kAm_withB;
q(16,14) = kBm_withA;
q(16,15) = kBm_withA;
q(16,17) = beta1;
q(17,16) = alpha1;
q(17,18) = beta2;
q(18,17) = alpha2;

% openstates
A = [17, 18];

% shutstates
F = 1:16;

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

% % Physical constraints
% 
% % Constrain rates k+s
% src = sub2ind(size(q), 1*ones(1,7), 9*ones(1,7));
% tgt = [1, 5;
%     2, 10;
% 	3, 7;
% 	5, 13;
% 	6, 14;
% 	9, 13;
% 	11, 15];
% tgt = sub2ind(size(q), tgt(:,1), tgt(:,2));
% [tmpg,tmpxi] = constrainrate (q, idxall, 'constrain', src, tgt, ones(1,7));
% gamma = [gamma; tmpg];
% xi = [xi; tmpxi];
% idxconstrain = [idxconstrain; tgt];
% 
% % Constrain rates k-s
% src = sub2ind(size(q), 1*ones(1,7), 3*ones(1,7));
% tgt = [1, 2;
%     2, 4;
%     3, 4;
%     5, 7;
%     6, 8;
%     9, 10;
%     11, 12];
% tgt = sub2ind(size(q), tgt(:,1), tgt(:,2));
% [tmpg,tmpxi] = constrainrate (q, idxall, 'constrain', src, tgt, ones(1,7));
% gamma=[gamma;tmpg];
% xi=[xi;tmpxi];
% idxconstrain = [idxconstrain; tgt];
% 
% % Constrain rates k+s'
% src = sub2ind(size(q), 3*ones(1,7), 11*ones(1,7));
% tgt = [2, 6;
% 	4, 8;
% 	4, 12;
% 	7, 15;
% 	8, 16;
% 	10, 14;
% 	12, 16];
% tgt = sub2ind(size(q), tgt(:,1), tgt(:,2));
% [tmpg,tmpxi] = constrainrate (q, idxall, 'constrain', src, tgt, ones(1,7));
% gamma=[gamma;tmpg];
% xi=[xi;tmpxi];
% idxconstrain = [idxconstrain; tgt];
% 
% % % Constrain rates k-s'
% src = sub2ind(size(q), 9*ones(1,7), 11*ones(1,7));
% tgt = [	5, 6;
% 	7, 8;
% 	10, 12;
% 	13, 14;
% 	13, 15;
%     14, 16;
% 	15, 16];
% tgt = sub2ind(size(q), tgt(:,1), tgt(:,2));
% [tmpg,tmpxi] = constrainrate (q, idxall, 'constrain', src, tgt, ones(1,7));
% gamma=[gamma;tmpg];
% xi=[xi;tmpxi];
% idxconstrain = [idxconstrain; tgt];
% 
% % Constrain rates k+f
% src = sub2ind(size(q), 9*ones(1,7), 1*ones(1,7));
% tgt = [	5, 1;
% 	7, 3;
% 	10, 2;
% 	13, 5;
% 	13, 9;
% 	14, 6;
% 	15, 11];
% tgt = sub2ind(size(q), tgt(:,1), tgt(:,2));
% [tmpg,tmpxi] = constrainrate (q, idxall, 'constrain', src, tgt, ones(1,7));
% gamma=[gamma;tmpg];
% xi=[xi;tmpxi];
% idxconstrain = [idxconstrain; tgt];
% 
% % Constrain rates k-f
% src = sub2ind(size(q), 11*ones(1,7), 3*ones(1,7));
% tgt = [	6, 2;
% 	8, 4;
% 	12, 4;
% 	14, 10;
% 	15, 7;
% 	16, 8;
% 	16, 12];
% tgt = sub2ind(size(q), tgt(:,1), tgt(:,2));
% [tmpg,tmpxi] = constrainrate (q, idxall, 'constrain', src, tgt, ones(1,7));
% gamma=[gamma;tmpg];
% xi=[xi;tmpxi];
% idxconstrain = [idxconstrain; tgt];
% 
% % Constrain rates k+f'
% src = sub2ind(size(q), 3*ones(1,7), 1*ones(1,7));
% tgt = [2, 1;
%     4, 2;
%     4, 3;
%     7, 5;
%     8, 6;
%     10, 9;
%     12, 11];
% tgt = sub2ind(size(q), tgt(:,1), tgt(:,2));
% [tmpg,tmpxi] = constrainrate (q, idxall, 'constrain', src, tgt, ones(1,7));
% gamma=[gamma;tmpg];
% xi=[xi;tmpxi];
% % need to transpose the tgt vector to add it to idxconstrain
% idxconstrain = [idxconstrain; tgt];
% 
% % Constrain rates k-f'
% src = sub2ind(size(q), 11*ones(1,7), 9*ones(1,7));
% tgt = [6, 5;
%     8, 7;
%     12, 10;
%     14, 13;
%     15, 13;
%     16, 14;
%     16, 15];
% tgt = sub2ind(size(q), tgt(:,1), tgt(:,2));
% [tmpg,tmpxi] = constrainrate (q, idxall, 'constrain', src, tgt, ones(1,7));
% gamma=[gamma;tmpg];
% xi=[xi;tmpxi];
% idxconstrain = [idxconstrain; tgt];

%% Microscopic Reversibility
[gamma,xi,idxconstrain,idxMR] = mrconstraints(q,idxall,idxconstrain,[],gamma,xi);

% %% Choose the rates to be constrained
% % create idxvary from the rates that are not constrained
% idxvary = setdiff(idxall,idxconstrain);
% 
% % reorder idxall and the columns of gamma so that the constrained rates are
% % first
% % [~,iAllToTheta] = ismember([idxconstrain;idxvary],idxall);
% % gamma = gamma(:,iAllToTheta);
% % idxall = [idxconstrain;idxvary];
% % clear iAllToTheta;
% 
% % Find rates to constrain or vary
% % This code works, but it might be better to use reduced row echelon form
% % [U,R,e] = qr(gamma,'vector');
% % idxvary = idxall(e(end-12:end));
% 
% [~,idx] = rref(gamma');
% gamma = gamma(idx,:);
% xi = xi(idx);
% 
% % Use reduced row echelon form to select a linearly independent set of
% % columns from gamma -> these will be a basis for the range of gamma
% % Then rearrange gamma so that the independent columns are the first
% % columns
% [~,jb] = rref(gamma);
% % gamma may have more rows than columns, meaning the user tried to set too
% % many physical constraints so there are now some non-independent
% % constraints. However, the rank of gamma should still be less than the
% % number of rates in the gating mechanism, otherwise there will be nothing
% % to fit because all the rates will be constrained.
% % jc = setdiff(1:length(gamma),jb); % This won't work in the above case
% jc = setdiff(1:numel(idxall),jb);
% idxconstrain = idxall(jb);
% idxvary = idxall(jc);
% idxall = idxall([jb,jc]);
% gamma = gamma(:,[jb,jc]);
% 
% % % Solve the underdetermined system for the minimum norm solution (i.e. find
% % % [one of the infinitely many sets of] rates satisfying the constraints)
% % [Q,R] = qr(gamma');
% % [~,R1] = qr(gamma',0);
% % x = Q*[R1'\xi;zeros(length(Q)-length(R1),1)];
% % 
% % q(idxall) = 10.^x;
% 
% [Q,R] = qr(gamma);
% R1 = R(:,1:rank(gamma));
% R2 = R(:,(rank(gamma)+1):numel(idxall));
% x2 = log10(q(idxvary));
% x1 = R1\(Q'*xi - R2*x2);
% q(idxall) = 10.^[x1;x2];