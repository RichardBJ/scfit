%% Load model
Gibb_model_Trihet_AB;

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

tcrit = zeros(length(file),1);
rates_tcrit = cell(length(file),1);
ll_tcrit = zeros(length(file),1);
qnew_tcrit = cell(length(file),1);
runtime_tcrit = zeros(length(file),1);
nDwells_tcrit = zeros(length(file),1);
eq_Po_tcrit = cell(length(file),1);
opentaus_tcrit = cell(length(file),1);
openareas_tcrit = cell(length(file),1);
shuttaus_tcrit = cell(length(file),1);
shutareas_tcrit = cell(length(file),1);

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

%% Calculate tcrit
shuts = rd(rs==0);
[t,w,ll] = emdistfit(shuts,[0.1 1 100],(1/3)*[1 1 1]);
tcrit(ii) = findtcrit(t,w);

%% Split the activity into bursts
% tcrit(ii)=35;
[bursts, bstates]=imposetcrit(rd,rs,tcrit(ii));
bursts(:,all(isnan(bursts)))=[];
bstates(:,all(isnan(bstates)))=[];

%% Fit data
nDwells(ii) = numel(rd);

tic;
% Fit without tcrit
[rates{ii}, ll(ii), qnew{ii}, hessian, covarmat, corrmat] = hjcfit(rd,q,A,F,td,idxall,idxvary,gamma,xi);

runtime(ii) = toc;

nDwells_tcrit(ii) = sum(sum(~isnan(bursts)));
tic;
% Fit bursts separated by at least tcrit
[rates_tcrit{ii}, ll_tcrit(ii), qnew_tcrit{ii}, hessian, covarmat, corrmat] = ...
    hjcfit(bursts,q,A,F,td,idxall,idxvary,gamma,xi,[],[],tcrit(ii));
runtime_tcrit(ii) = toc;

%% Save results
eq_Po{ii} = eqoccupy(qnew{ii});
[opentaus{ii}, openareas{ii}, shuttaus{ii}, shutareas{ii}] = ...
    qmatopenshutpdf(qnew{ii},A,F);
write_qmatrix_results('Gibb_model_Trihet_AB.m',[path file{ii}],qnew{ii},A,F,idxall,idxvary,...
    td,Inf,nDwells(ii),runtime(ii),ll(ii),eq_Po{ii},opentaus{ii},...
    openareas{ii}, shuttaus{ii}, shutareas{ii});
eq_Po_tcrit{ii} = eqoccupy(qnew_tcrit{ii});
[opentaus_tcrit{ii}, openareas_tcrit{ii}, shuttaus_tcrit{ii}, shutareas_tcrit{ii}] = ...
    qmatopenshutpdf(qnew_tcrit{ii},A,F);
pause(2);
write_qmatrix_results('Gibb_model_Trihet_AB.m',[path file{ii}],qnew_tcrit{ii},A,F,idxall,idxvary,...
    td,tcrit(ii),nDwells_tcrit(ii),runtime_tcrit(ii),ll_tcrit(ii),eq_Po_tcrit{ii},...
    opentaus_tcrit{ii}, openareas_tcrit{ii}, shuttaus_tcrit{ii}, shutareas_tcrit{ii});

end