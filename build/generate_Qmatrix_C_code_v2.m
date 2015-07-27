% Qmatrixloglik (params, Q, idxtheta, M, b, A, F, td, dwells)
% Qloglik_bursts (params, Q, idxtheta, M, b, A, F, ...
%     td, tcrit, dwells)
ms=100;
ms2 = ms*ms;
params = coder.typeof(1,[ms2 1],1);
Q = coder.typeof(1,[ms ms],[1 1]);
% idxAll = coder.typeof(1,[ms2 1],1);
% idxvarytoall = coder.typeof(1,[ms2 1],1);
% idxcontoall = coder.typeof(1,[ms2 1],[1 1]);
% R1 = coder.typeof(1,[ms ms],[1 1]);
% R2 = coder.typeof(1,[ms ms],[1 1]);
% U = coder.typeof(1,[ms ms],[1 1]);
idxtheta = coder.typeof(1,[ms2 1],[1 0]);
M = coder.typeof(1,[ms2 ms2],[1 1]);
b = coder.typeof(1,[ms2 1],[1 0]);
xi = coder.typeof(1,[ms2 1],[1 1]);
A = coder.typeof(1,[1 ms],[0 1]);
F = coder.typeof(1,[1 ms],[0 1]);
td = 1;

% for Qmatrixloglik
fprintf('Generating Qmatrixloglik ...\n');
dwells = coder.typeof(1,[inf 1],[1 0]);
variables = {params, Q, idxtheta, M, b, A, F, td, dwells};
codegen Qmatrixloglik -args variables
fprintf('Qmatrixloglik Done!\n');

% for Qloglik_bursts
fprintf('Generating Qloglik_bursts ...\n');
dwells = coder.typeof(1,[inf inf],[1 1]);
variables = {params, Q, idxtheta, M, b, A, F, td, td, dwells};
codegen Qloglik_bursts -args variables
fprintf('Qloglik_bursts Complete!\n');

% for hjcdist
% function [ pdf ] = hjcdist( q, A, F, td, t, varargin )
fprintf('Generating hjcdist ...\n');
t = coder.typeof(1,[1e5 1e5],[1 1]);
variables = {Q,A,F,td,t,true};
codegen hjcdist -args variables 