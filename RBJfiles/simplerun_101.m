%%Demo by Kevin

%% Load the model
[q, A, F, idxall, gamma, xi, idxvary, idxconstrain, modelfile, ...
    statenames] = loadmodel;

%% Read QuB idealizations in dwt format
[idl_file, path] = uigetfile('*.dwt', 'MultiSelect', 'on');
[durations, amp] = dwtread([path, idl_file]);

%% Impose the resolution
td = 0.05;
[rd, rs] = imposeres(durations, amp, td, td);

%% Concatenate dwell times into contiguous open periods
% pA for real difference is needed for SCN data but not DWT
pA_for_real_diff = 2 * max(abs(rs));
zero_amp = 0;
[rd, rs] = concatdwells(rd, rs, pA_for_real_diff, zero_amp);

%% Fit rates
tic;
[rates, ll, qnew] = hjcfit(rd, q, A, F, td, idxall, idxvary, gamma, xi);
runtime = toc;

%% Calculate channel properties
eq_Po = eqoccupy(qnew);
[opentaus, openareas, shuttaus, shutareas] = qmatopenshutpdf(qnew, A, F);

%% Save results to CSV file
write_qmatrix_results(modelfile, [path, idl_file], qnew, A, F, ...
    idxall, idxvary, td, Inf, nDwells, runtime, ll, ...
    eq_Po, opentaus, openareas, shuttaus, shutareas);
