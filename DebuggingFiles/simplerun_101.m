%%Demo by Kevin Ogden
% Need to add the subfolders to path first
% use prepare.m in project root 

%% Load the model
[q, A, F, idxall, gamma, xi, idxvary, idxconstrain, modelfile, ...
    statenames] = loadmodel;

%% Read QuB idealizations in dwt format
idl_file = uigetfile('*.dwt; *.txt; *.csv; *.scn', 'MultiSelect', 'on');

[path,file,ext] = fileparts(idl_file);
file=[fullfile(path,file),ext];
if ext==".dwt"
    [durations, amp] = dwtread(file);
elseif ext==".scn"
        disp("simplerun not working for scn.... I don't use these to difficult to debug");
        disp("sorry, stopping here!");
        return
        [ durations, amplitudes, properties, calibration, varargout] = scanread( file );
else 
    %need to feed it sample interval/dt or si
    si=0.01; %ms
    [durations, amp] = scfitcsvread(file,si);  
end

%% Impose the resolution
td = 0.05;
[rd, rs] = imposeres(durations, amp, td, td);

%% Concatenate dwell times into contiguous open periods
% pA for real difference is needed for SCN data but not DWT
pA_for_real_diff = 2 * max(abs(rs));
zero_amp = 0;
[rd, rs] = concatdwells(rd, rs, pA_for_real_diff, zero_amp);

% I was getting an error with my csv's if last state was 0, so crop that
% off if it happens.  Not sure if this would happen with a dwt/scn.
if rs(end)==0
    rd(end)=[];
    rs(end)=[];
end

%% Fit rates
tic;
[rates, ll, qnew] = hjcfit(rd, q, A, F, td, idxall, idxvary, gamma, xi);
runtime = toc;

%% Calculate channel properties
eq_Po = eqoccupy(qnew);
[opentaus, openareas, shuttaus, shutareas] = qmatopenshutpdf(qnew, A, F);

%% Save results to CSV file
nDwells=length(rs);
write_qmatrix_results(modelfile, file, qnew, A, F, ...
    idxall, idxvary, td, Inf, nDwells, runtime, ll, ...
    eq_Po, opentaus, openareas, shuttaus, shutareas);
