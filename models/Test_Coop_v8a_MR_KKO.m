%% NMDA Receptor Model - Test_Coop_v8a_MR
%
% Build a model of NMDA receptor activation
%
k12=	2.5	;%MR
k21=	40	;
alpha1=	3000	;
beta1=	10000	;
alpha2=	12000	;
beta2=	5000	;
kfp=	500	;
kfm=	1500	;
ksp=	8.3333	;
ksm=	20	;
kftp=	800	;
kftm=	2000	;
kstp=	5	;
kstm=	10	;
%
q=zeros(12);
%
q(1,4) = alpha1;
q(4,1) = beta1;
q(2,5) = alpha2;
q(5,2) = beta2;
q(2,1) = k21;
q(1,2) = k12;
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
A = 1:2;

% shutstates
F = 3:12;

% scale rates to ms^-1
q = q*1e-3;

% q = q-diag(sum(q,2));

%% Constraints
% Description of first code block
idxall = find(q);
% It may help to write out the rates that vary in a different way.
% In this syntax, the alignment is just for visualization purposes. The
% important things are that there is a blank space or a comma between the
% rates in a pair and that between pairs there is a semicolon.
ratesToVary = [ 1 4;
                4 1;
                1 2;
%                 2 1;  This rate is constrained in lines 201-205
                2 5;
                5 2;
                5 4;
                6 9;
                6 4;
                10 7;
                10 12;
                11 9;
                11 12 ];

% It may be easier to write this variable as on the next line
% ratesToVary = [1 4; 4 1; 1 2; 2 1; 2 5; 5 2; 5 4; 6 9; 6 4; 10 7; 10 12; 11 9; 11 12];

idxvary = sub2ind(size(q), ratesToVary(:,1), ratesToVary(:,2));

% They both are equivalent to what we were using before
% idxvary = sub2ind(size(q), [1 4 1 2 2 5 5 6 6 10 10 11 11], [4 1 2 1 5 2 4 9 4 7 12 9 12]);


% Constrain rates k+s
src = sub2ind(size(q), 11*ones(1,3), 9*ones(1,3));
tgt = sub2ind(size(q), [10 8 12], [3 6 11]);
[gamma,xi] = constrainrate (q, idxall, 'constrain', src, tgt, [1 1 2]);

% Constrain rates with microscopic reversibility
% KKO 11 Aug 2014
% Constrain the loop containing states 8, 11, 12, and 10
% src = sub2ind(size(q),[8 10 12 11],[10 12 11 8]); %Anti-clockwise
% tgt = sub2ind(size(q),[8 11 12 10],[11 12 10 8]); %Clockwise
% [tmpg,tmpxi] = constrainrate(q,idxall,'loop',src,tgt);
% gamma=[gamma;tmpg];
% xi=[xi;tmpxi];

% Constrain the loop containing states 1, 4, 5, and 2
src = sub2ind(size(q),[1 2 5 4],[2 5 4 1]); %Anti-clockwise
tgt = sub2ind(size(q),[1 4 5 2],[4 5 2 1]); %Clockwise
[tmpg,tmpxi] = constrainrate(q,idxall,'loop',src,tgt);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];

% Constrain the loop containing states 5, 7, 8, and 10
src = sub2ind(size(q),[5 7 10 8],[7 10 8 5]); %Anti-clockwise
tgt = sub2ind(size(q),[5 8 10 7],[8 10 7 5]); %Clockwise
[tmpg,tmpxi] = constrainrate(q,idxall,'loop',src,tgt);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];

% Constrain the loop containing states 4, 6, 8, and 5
% src = sub2ind(size(q),[4 6 8 5],[6 8 5 4]); %Anti-clockwise
% tgt = sub2ind(size(q),[4 5 8 6],[5 8 6 4]); %Clockwise
% [tmpg,tmpxi] = constrainrate(q,idxall,'loop',src,tgt);
% gamma=[gamma;tmpg];
% xi=[xi;tmpxi];

% Constrain the loop containing states 6, 9, 11, and 8
src = sub2ind(size(q),[6 8 11 9],[8 11 9 6]); %Anti-clockwise
tgt = sub2ind(size(q),[6 9 11 8],[9 11 8 6]); %Clockwise
[tmpg,tmpxi] = constrainrate(q,idxall,'loop',src,tgt);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];

% Constrain the loop containing states 6, 3, 11, and 9
% src = sub2ind(size(q),[6 9 11 3],[9 11 3 6]); %Anti-clockwise
% tgt = sub2ind(size(q),[6 3 11 9],[3 11 9 6]); %Clockwise
% [tmpg,tmpxi] = constrainrate(q,idxall,'loop',src,tgt);
% gamma=[gamma;tmpg];
% xi=[xi;tmpxi];

% Constrain the loop containing states 6, 3, 5, and 4
% src = sub2ind(size(q),[6 3 5 4],[3 5 4 6]); %Anti-clockwise
% tgt = sub2ind(size(q),[6 4 5 3],[4 5 3 6]); %Clockwise
% [tmpg,tmpxi] = constrainrate(q,idxall,'loop',src,tgt);
% gamma=[gamma;tmpg];
% xi=[xi;tmpxi];

% Constrain the loop containing states 3, 5, 7, and 10
% src = sub2ind(size(q),[3 10 7 5],[10 7 5 3]); %Anti-clockwise
% tgt = sub2ind(size(q),[3 5 7 10],[5 7 10 3]); %Clockwise
% [tmpg,tmpxi] = constrainrate(q,idxall,'loop',src,tgt);
% gamma=[gamma;tmpg];
% xi=[xi;tmpxi];

% Constrain the loop containing states 3, 10, 12, and 11
src = sub2ind(size(q),[3 11 12 10],[11 12 10 3]); %Anti-clockwise
tgt = sub2ind(size(q),[3 10 12 11],[10 12 11 3]); %Clockwise
[tmpg,tmpxi] = constrainrate(q,idxall,'loop',src,tgt);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];

% Constrain the interior loop containing states 3, 6, 8, and 10
src = sub2ind(size(q),[3 6 8 10],[6 8 10 3]); %Anti-clockwise
tgt = sub2ind(size(q),[3 10 8 6],[10 8 6 3]); %Clockwise
[tmpg,tmpxi] = constrainrate(q,idxall,'loop',src,tgt);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];

% Constrain the interior loop containing states 3, 11, 8, and 5
src = sub2ind(size(q),[3 5 8 11],[5 8 11 3]); %Anti-clockwise
tgt = sub2ind(size(q),[3 11 8 5],[11 8 5 3]); %Clockwise
[tmpg,tmpxi] = constrainrate(q,idxall,'loop',src,tgt);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];

% Constrain the outside loop containing states 4,5,7,10,12,11,9 and 6
src = sub2ind(size(q),[4 5 7 10 12 11 9 6],[5 7 10 12 11 9 6 4]); %Anti-clockwise
tgt = sub2ind(size(q),[4 6 9 11 12 10 7 5],[6 9 11 12 10 7 5 4]); %Clockwise
[tmpg,tmpxi] = constrainrate(q,idxall,'loop',src,tgt);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];

% Constrain rates k-s
src = sub2ind(size(q), 11*ones(1,1), 12*ones(1,1));
tgt = sub2ind(size(q), [9], [11]);
[tmpg,tmpxi] = constrainrate (q, idxall, 'constrain', src, tgt, [2]);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];

% Constrain rates k+s'
src = sub2ind(size(q), 5*ones(1,4), 4*ones(1,4));
tgt = sub2ind(size(q), [2 3 10 7], [1 6 8 5]);
[tmpg,tmpxi] = constrainrate (q, idxall, 'constrain', src, tgt, [1 1 1 2]);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];

% % Constrain rates k-s'
% src = sub2ind(size(q), 5*ones(1,3), 7*ones(1,3));
% tgt = sub2ind(size(q), [6 8 4], [3 10 5]);
% [tmpg,tmpxi] = constrainrate (q, idxall, 'constrain', src, tgt, [1 1 2]);
% gamma=[gamma;tmpg];
% xi=[xi;tmpxi];

% Constrain rates k+f
src = sub2ind(size(q), 10*ones(1,3), 7*ones(1,3));
tgt = sub2ind(size(q), [8 11 12], [5 3 10]);
[tmpg,tmpxi] = constrainrate (q, idxall, 'constrain', src, tgt, [1 1 2]);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];

% Constrain rates k-f
src = sub2ind(size(q), 10*ones(1,3), 12*ones(1,3));
tgt = sub2ind(size(q), [3 5 7], [11 8 10]);
[tmpg,tmpxi] = constrainrate (q, idxall, 'constrain', src, tgt, [1 1 2]);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];

% Constrain rates k+f'
src = sub2ind(size(q), 6*ones(1,3), 4*ones(1,3));
tgt = sub2ind(size(q), [3 11 9], [5 8 6]);
[tmpg,tmpxi] = constrainrate (q, idxall, 'constrain', src, tgt, [1 1 2]);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];

% Constrain rates k-f'
src = sub2ind(size(q), 6*ones(1,2), 9*ones(1,2));
tgt = sub2ind(size(q), [8 4], [11 6]);
[tmpg,tmpxi] = constrainrate (q, idxall, 'constrain', src, tgt, [1 2]);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];