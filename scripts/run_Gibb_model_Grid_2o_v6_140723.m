%% Load model
Gibb_model_Grid_2o_v6a;
q(13,12)=1;
q(14,13)=2;

%% Get dwt files
% [file, path] = uigetfile('*.dwt','MultiSelect','on');
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

%% Run loop

rates = cell(length(file),1);
ll = zeros(length(file),1);
qnew = cell(length(file),1);
runtime = zeros(length(file),1);
nDwells = zeros(length(file),1);
eq_Po = cell(length(file),1);
opentaus = cell(length(file),1);
openareas = cell(length(file),1);
shuttaus = cell(length(file),1);
shutareas = cell(length(file),1);

for ii=1:length(file)
    %% Read idealized data
if DWTFILE
    [durations, amp] = dwtread([path file{ii}]);
    if isinteger(amp)
        amp = double(amp);
    end
else 
    [durations, amp, ~, cal] = scanread([path file{ii}]);
    amp=amp*cal;
end

%% Impose the resolution
td=0.04;
[rd,rs] = imposeres (durations,amp,td,td);

%% Concatenate dwell times into contiguous open periods
% Needed for SCN data, but not for DWT
if ~DWTFILE
    pA_for_real_diff = 2*max(abs(rs));
    zeroAmp = 0;
    [rd, rs] = concatdwells( rd, rs, pA_for_real_diff, zeroAmp);
end

%% Make sure idealization starts with opening and ends with a closing
if(rs(1)==0)
    rd(1)=[];rs(1)=[];
end
if(rs(end)~=0)
    rd(end)=[];rs(end)=[];
end
assert(mod(numel(rd),2)==0);

nDwells = numel(rd);

%% Fit data
tic;

% Fit without tcrit
[rates{ii}, ll(ii), qnew{ii}, hessian, covarmat, corrmat] = hjcfit(rd,q,A,F,td,idxall,idxvary,gamma,xi);
    
runtime(ii) = toc;

%% Save results
eq_Po{ii} = eqoccupy(qnew{ii});
[opentaus{ii}, openareas{ii}, shuttaus{ii}, shutareas{ii}] = ...
    qmatopenshutpdf(qnew{ii},A,F);
write_qmatrix_results('Gibb_model_Grid_2o_v6a.m',[path file{ii}],qnew{ii},A,F,idxall,idxvary,...
    td,Inf,nDwells(ii),runtime(ii),ll(ii),eq_Po{ii},opentaus{ii},...
    openareas{ii}, shuttaus{ii}, shutareas{ii});

end