%% run_Gibb_model_Grid_2o_v4
Gibb_model_Grid_2o_v4

%% Get Idealized Files
[file, path, filterIndex] = uigetfile({'*.dwt','QuB Idealized Data';'*.scn','SCAN idealized data'},'MultiSelect','on');
if filterIndex==1
    DWTFILE=true;
elseif filterIndex==2
    DWTFILE=false;
end
if ~iscellstr(file)
    if ischar(file)
        file = cellstr(file);
    else
        error('The file name should be a cell string. Try entering it manually');
    end
end

%% Get idealized data
% Read in idealized data from all files and concatenate into a single list
% of dwell times
% 
% Assume that the transitions from each file are continuous in time

durations = [];
amp = [];
for ii=1:length(file)
    if DWTFILE
        [tmp_duration, tmp_amp] = dwtread([path, file{ii}]);
        if isinteger(tmp_amp)
            tmp_amp = double(tmp_amp);
        end
    else
        [tmp_duration, tmp_amp, ~, cal] = scanread([path, file{ii}]);
        tmp_amp = tmp_amp*cal;
    end
    
    durations = [durations; tmp_duration];
    amp = [amp; tmp_amp];
end

%% Impose the resoltion
openResolution = 0.045; %ms
shutResolution = 0.045;
[resolvedDwells, resolvedStates] = imposeres(durations, amp, ...
    openResolution, shutResolution);

%% Concatenate dwell times into contiguous open periods
% Needed for SCN data, but not for DWT
if ~DWTFILE
    pA_for_real_diff = 2*max(abs(resolvedStates));
    zeroAmp = 0;
    [resolvedDwells, resolvedStates] = concatdwells(resolvedDwells, ...
        resolvedStates, pA_for_real_diff, zeroAmp);
end

%% Make sure idealization starts with an opening and ends with a closing
% Eventually, I should re-write the hjcfit function so accomodate data that
% starts in any state (not necessarily open) and ends with a resolved dwell
% in any state (not necessarily closed)

if (resolvedStates(1) == 0)
    resolvedDwells(1) = [];
    resolvedStates(1) = [];
end
if (resolvedStates(end) ~= 0)
    resolvedDwells(end) = [];
    resolvedStates(end) = [];
end
nDwells = numel(resolvedDwells);
assert(mod(nDwells,2) == 0);

%% Fit data
tic;

% Fit without tcrit
[rates, ll, qnew, hessian, covarmat, corrmat] = hjcfit(resolvedDwells, ...
    q, A, F, openResolution, idxall, idxvary, gamma, xi);

runtime = toc;

%% Save fitting results
eq_Po = eqoccupy(qnew);
[opentaus, openareas, shuttaus, shutareas] = qmatopenshutpdf(qnew, A, F);
write_qmatrix_results('Gibb_model_Grid_2o_v4',[path file{1}],qnew,A,F,idxall,idxvary,...
    openResolution,Inf(),nDwells,runtime,ll,eq_Po,opentaus, openareas, shuttaus, shutareas);
