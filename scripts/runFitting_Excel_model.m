%% Load the model
[q, A, F, idxall, idxvary, gamma, xi, idxvary, idxconstrain, modelfile, ...
    statenames] = loadmodel;
% The rate constants need to correspond to the same units as the intervals
% (loaded later). Usually the intervals are in milliseconds, which means the
% rates should be in ms^-1 (per millisecond). If the rates in the Excel model
% need to be scaled to match the intervals, do so now, e.g.
% q = q*1e-3;

%% Get file names
% QuB dwt files
% [file, path] = uigetfile('*.dwt','MultiSelect','on');

% SCAN scn files
[idl_file, path] = uigetfile({'*.scn','SCAN idealized data'},'MultiSelect','on');

if ~iscellstr(idl_file)
    if ischar(idl_file)
        idl_file = cellstr(idl_file);
    else
        error('The idl_file name should be a cell string. Try entering it manually');
    end
end
rates = cell(length(idl_file),1);
ll = zeros(length(idl_file),1);
qnew = cell(length(idl_file),1);
runtime = zeros(length(idl_file),1);
nDwells = zeros(length(idl_file),1);
tcrit = zeros(length(idl_file),1);
eq_Po = cell(length(idl_idl_file), 1);
opentaus = cell(length(idl_file), 1);
openareas = cell(length(idl_file), 1);
shuttaus = cell(length(idl_file), 1);
shutareas = cell(length(idl_file), 1);

%% Run loop

for ii=1:length(idl_file);
%% Read idealized data
% Read QuB dwt idl_files    
% [durations, amp] = dwtread([path, idl_file{ii}]);
% if isinteger(amp)
%     amp = double(amp);
% end

% Read SCAN.scn idl_files
[durations, amp, properties, cal] = scanread([path, idl_file{ii}]);
amp=amp*cal;

%% Impose the resolution
td=0.05;
% Use the same resolution for open and shut durations
[rd,rs] = imposeres (durations,amp,td,td);

%% Concatenate dwell times into contiguous open periods
% Needed for SCN data, but not for DWT
pA_for_real_diff = 2*max(abs(rs));
zeroAmp = 0;
[rd, rs] = concatdwells(rd, rs, pA_for_real_diff, zeroAmp);


%% Calculate tcrit
% shuts = rd(rs==0);
% [t,w,ll] = emdistfit(shuts,[0.1 1 100],(1/3)*[1 1 1]);
% tcrit(ii) = findtcrit(t,w);

%% Split the activity into bursts
% tcrit=100;
% [bursts, bstates]=imposetcrit(rd,rs,tcrit);
% bursts(:,all(isnan(bursts)))=[];
% bstates(:,all(isnan(bstates)))=[];
% nDwells(ii) = sum(sum(~isnan(bursts)));

% Make sure idealization starts with opening and ends with a closing
if(rs(1)~=1)
    rd(1)=[];rs(1)=[];
end
if(rs(end)~=0)
    rd(end)=[];rs(end)=[];
end
assert(mod(numel(rd),2)==0);

%% Fit data
tic;

% Fit without tcrit
[rates{ii}, ll(ii), qnew{ii}] = hjcfit(rd,q,A,F,td,idxall,idxvary,gamma,xi);
    
% Fit bursts separated by at least tcrit
% [rates{ii}, ll(ii), qnew{ii}, hessian, covarmat, corrmat] = ...
%     hjcfit(bursts,q,A,F,td,idxall,idxvary,gamma,xi,[],[],tcrit(ii));

runtime(ii) = toc;

%% Save results
eq_Po{ii} = eqoccupy(qnew{ii});
[opentaus{ii}, openareas{ii}, shuttaus{ii}, shutareas{ii}] = ...
    qmatopenshutpdf(qnew{ii}, A, F);
write_qmatrix_results(modelfile, [path, idl_file{ii}], qnew{ii}, A, F, ...
    idxall, idxvary, td, Inf, nDwells(ii), runtime(ii), ll(ii), ...
    eq_Po{ii}, opentaus{ii}, openareas{ii}, shuttaus{ii}, shutareas{ii});
end
