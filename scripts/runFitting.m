%% Load the model
% Gibb_model;
[q, A, F, idxall, idxvary, gamma, xi, fname] = qmatxlsread;
q = q*1e-3;

%% Get file names
% QuB dwt files
% [file, path] = uigetfile('*.dwt','MultiSelect','on');

% SCAN scn files
[file, path] = uigetfile({'*.scn','SCAN idealized data'},'MultiSelect','on');

if ~iscellstr(file)
    if ischar(file)
        file = cellstr(file);
    else
        error('The file name should be a cell string. Try entering it manually');
    end
end
rates = cell(length(file),1);
ll = zeros(length(file),1);
qnew = cell(length(file),1);
runtime = zeros(length(file),1);
nDwells = zeros(length(file),1);
tcrit = zeros(length(file),1);

%% Run loop

for ii=1:length(file);
%% Read idealized data
% Read QuB dwt files    
% [durations, amp] = dwtread([path file{ii}]);
% if isinteger(amp)
%     amp = double(amp);
% end

% Read SCAN.scn files
[durations, amp, properties, cal] = scanread([path file{ii}]);
amp=amp*cal;

%% Impose the resolution
td=0.05;
[rd,rs] = imposeres (durations,amp,td,td);

%% Concatenate dwell times into contiguous open periods
% Needed for SCN data, but not for DWT
pA_for_real_diff = 2*max(abs(rs));
zeroAmp = 0;
[rd, rs] = concatdwells( rd, rs, pA_for_real_diff, zeroAmp);


%% Calculate tcrit
shuts = rd(rs==0);
[t,w,ll] = emdistfit(shuts,[0.1 1 100],(1/3)*[1 1 1]);
tcrit(ii) = findtcrit(t,w);

%% Split the activity into bursts
% tcrit=100;
[bursts, bstates]=imposetcrit(rd,rs,tcrit);
bursts(:,all(isnan(bursts)))=[];
bstates(:,all(isnan(bstates)))=[];
nDwells(ii) = sum(sum(~isnan(bursts)));

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
%[rates{ii}, ll(ii), qnew{ii}] = hjcfit(rd,q,A,F,td,idxall,idxvary,gamma,xi);
    
% Fit bursts separated by at least tcrit
[rates{ii}, ll(ii), qnew{ii}, hessian, covarmat, corrmat] = ...
    hjcfit(bursts,q,A,F,td,idxall,idxvary,gamma,xi,[],[],tcrit(ii));

runtime(ii) = toc;
end