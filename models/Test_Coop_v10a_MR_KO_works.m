%% NMDA Receptor Model - Test_Coop_v10a_MR
%
% Build a model of NMDA receptor activation
%
alpha3=	2000	;
beta3=	1000	;
alpha1=	800	;
beta1=	6000	;
alpha2=	4000	;
beta2=	1000	;
kfp=	500	;
kfm=	1500	;
ksp=	10	;
ksm=	20	;
kftp=	444.444	;
kftm=	2000	;
kstp=	200	;
kstm=	600	;
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
% Description of first code block
idxall = find(q);
% First define the rates to be constrained and then set the variable rates
% at the end

% declare the idxconstrain variable. It's empty now, but we'll add rates to it as we go.
idxconstrain = []; 

% Constrain the loop containing states 1, 4, 5, and 2
src = sub2ind(size(q),[1 2 5 4],[2 5 4 1]); %Anti-clockwise
tgt = sub2ind(size(q),[1 4 5 2],[4 5 2 1]); %Clockwise
[gamma,xi] = constrainrate(q,idxall,'loop',src,tgt);

% Constrain rate from states 4 to 1
idxconstrain = [idxconstrain; sub2ind(size(q),4,1)];

% Constrain the loop containing states 11, 3, 6, and 9
src = sub2ind(size(q),[11 3 6 9],[3 6 9 11]); %Anti-clockwise
tgt = sub2ind(size(q),[11 9 6 3],[9 6 3 11]); %Clockwise
[tmpg,tmpxi] = constrainrate(q,idxall,'loop',src,tgt);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];
% Constrain rate from state 11 to 3
idxconstrain = [idxconstrain; sub2ind(size(q),11,3)];

% Constrain the loop containing states 6, 4, 5, and 3
src = sub2ind(size(q),[6 4 5 3],[4 5 3 6]); %Anti-clockwise
tgt = sub2ind(size(q),[6 3 5 4],[3 5 4 6]); %Clockwise
[tmpg,tmpxi] = constrainrate(q,idxall,'loop',src,tgt);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];
% Constrain rate from state 6 to 4
idxconstrain = [idxconstrain; sub2ind(size(q),6,4)];

% Constrain the loop containing states 8, 5, 4, and 6
src = sub2ind(size(q),[8 5 4 6],[5 4 6 8]); %Anti-clockwise
tgt = sub2ind(size(q),[8 6 4 5],[6 4 5 8]); %Clockwise
[tmpg,tmpxi] = constrainrate(q,idxall,'loop',src,tgt);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];
% Constrain rate from state 8 to 5
idxconstrain = [idxconstrain; sub2ind(size(q),8,5)];

% Constrain the loop containing states 9, 6, 8 and 11
src = sub2ind(size(q),[9 6 8 11],[6 8 11 9]); %Anti-clockwise
tgt = sub2ind(size(q),[9 11 8 6],[11 8 6 9]); %Clockwise
[tmpg,tmpxi] = constrainrate(q,idxall,'loop',src,tgt);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];
% Constrain rate from state 9 to 6
idxconstrain = [idxconstrain; sub2ind(size(q),9,6)];

% Constrain the loop containing states 13, 1, 4, and 6
src = sub2ind(size(q),[13 1 4 6],[1 4 6 13]); %Anti-clockwise
tgt = sub2ind(size(q),[13 6 4 1],[6 4 1 13]); %Clockwise
[tmpg,tmpxi] = constrainrate(q,idxall,'loop',src,tgt);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];
% Constrain rate from state 13 to 6
idxconstrain = [idxconstrain; sub2ind(size(q),13,6)];

% Constrain the loop containing states 10, 7, 5, and 3
src = sub2ind(size(q),[10 7 5 3],[7 5 3 10]); %Anti-clockwise
tgt = sub2ind(size(q),[10 3 5 7],[3 5 7 10]); %Clockwise
[tmpg,tmpxi] = constrainrate(q,idxall,'loop',src,tgt);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];
% Constrain rate from state 10 to 7
idxconstrain = [idxconstrain; sub2ind(size(q),10,7)];

% Constrain the loop containing states 10, 8, 6, and 3
src = sub2ind(size(q),[10 8 6 3],[8 6 3 10]); %Anti-clockwise
tgt = sub2ind(size(q),[10 3 6 8],[3 6 8 10]); %Clockwise
[tmpg,tmpxi] = constrainrate(q,idxall,'loop',src,tgt);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];
% Constrain rate from state 10 to 8
idxconstrain = [idxconstrain; sub2ind(size(q),10,8)];

% Constrain the loop containing states 12, 10, 3, 11
src = sub2ind(size(q),[12 10 3 11],[10 3 11 12]); %Anti-clockwise
tgt = sub2ind(size(q),[12 11 3 10],[11 3 10 12]); %Clockwise
[tmpg,tmpxi] = constrainrate(q,idxall,'loop',src,tgt);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];
% Constrain rate from state 12 to 10
idxconstrain = [idxconstrain; sub2ind(size(q),12,10)];

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
src = sub2ind(size(q), 5*ones(1,3), 4*ones(1,3));
tgt = sub2ind(size(q), [2 3 7], [1 6 5]);
[tmpg,tmpxi] = constrainrate (q, idxall, 'constrain', src, tgt, [1 1 1 2]);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];
% need to transpose the tgt vector to add it to idxconstrain
idxconstrain = [idxconstrain; tgt'];

% % Constrain rates k-s'
src = sub2ind(size(q), 5*ones(1,4), 7*ones(1,4));
tgt = sub2ind(size(q), [6 8 4 1], [3 10 5 2]);
[tmpg,tmpxi] = constrainrate (q, idxall, 'constrain', src, tgt, [1 1 2 2]);
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
src = sub2ind(size(q), 6*ones(1,3), 4*ones(1,3));
tgt = sub2ind(size(q), [3 11 13], [5 8 1]);
[tmpg,tmpxi] = constrainrate (q, idxall, 'constrain', src, tgt, [1 1 1 2]);
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

% create idxvary from the rates that are not constrained
idxvary = setdiff(idxall,idxconstrain);

% I found the following rates by using the QR decomposition of gamma and
% found a set of linearly independent columns using
% [U,R,e] = qr(gamma,'vector');
% [ii,jj] = ind2sub(size(q),idxall(e)); tmp=[ii,jj];

idxvary = sub2ind(size(q),...
    [1	5	5	1	2	9	10	6	10	12	3	12	6],...
    [2	2	3	4	5	6	7	8	8	10	11	11	13]);

idxconstrain = sub2ind(size(q),...
    [2	4	13	6	10	11	5	6	3	4	7	8	3	4	8	13	5	5	11	6	11	3	7	8	8	9	10	11	1],...
    [1	1	1	3	3	3	4	4	5	5	5	5	6	6	6	6	7	8	8	9	9	10	10	10	11	11	12	12	13] );