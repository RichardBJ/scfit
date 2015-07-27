%% Load Qub Model File (qmf)
[q, A, F, idxall, idxvary, gamma, xi, numstates, ~, ~, ~, ~, filename] = qmfread;
% q = q*1e-3;
% q(idxvary) = exprnd(5,size(idxvary));

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

    %% Read idealized data
if DWTFILE
    [durations, amp] = dwtread([path file{1}]);
    if isinteger(amp)
        amp = double(amp);
    end
else 
    [durations, amp, ~, cal] = scanread([path file{1}]);
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
% shuts = rd(rs==0);
% [t,w,ll] = emdistfit(shuts,[0.1 1 100],(1/3)*[1 1 1]);
% tcrit(ii) = findtcrit(t,w);

%% Split the activity into bursts
tcrit=35;
[bursts, bstates]=imposetcrit(rd,rs,tcrit);
bursts(:,all(isnan(bursts)))=[];
bstates(:,all(isnan(bstates)))=[];
nDwells = sum(sum(~isnan(bursts)));

%% Fit data
numrestarts=50;
seed = zeros(numrestarts,1);
guess = zeros(numel(idxvary),numrestarts);
rates = cell(numrestarts,1);
ll = zeros(numrestarts,1);
qnew = cell(numrestarts,1);
hessian = cell(numrestarts,1);
covarmat = cell(numrestarts,1);
corrmat = cell(numrestarts,1);
runtime = zeros(numrestarts,1);
eq_Po = cell(numrestarts,1);
opentaus = cell(numrestarts,1);
openareas = cell(numrestarts,1);
shuttaus = cell(numrestarts,1);
shutareas = cell(numrestarts,1);

options=optimset('LargeScale','off', 'MaxFunEvals',5e3, 'Display','final');

for ii=1:numrestarts
    seed(ii) = normrnd(5,1);
    guess(:,ii) = exprnd(seed(ii),numel(idxvary),1);
    q(idxvary) = guess(:,ii);
    tic;

% Fit without tcrit
% [rates{ii}, ll(ii), qnew{ii}, hessian, covarmat, corrmat] = hjcfit(rd,q,A,F,td,idxall,idxvary,gamma,xi);
    
% Fit bursts separated by at least tcrit
try
    [rates{ii}, ll(ii), qnew{ii}, hessian{ii}, covarmat{ii}, corrmat{ii}] = ...
        hjcfit(bursts,q,A,F,td,idxall,idxvary,gamma,xi,[],options,tcrit);
catch err
%     rethrow(err);
    continue
end
% [~, lltest] = ...
%     hjcfit(bursts,qtest,A,F,td,idxall,idxvary,gamma,xi,false,[],tcrit(ii));
runtime(ii) = toc;

%% Save results
eq_Po{ii} = eqoccupy(qnew{ii});
[opentaus{ii}, openareas{ii}, shuttaus{ii}, shutareas{ii}] = ...
    qmatopenshutpdf(qnew{ii},A,F);
% write_qmatrix_results(filename,[path file{ii}],qnew{ii},A,F,idxall,idxvary,...
%     td,tcrit(ii),nDwells(ii),ll(ii),eq_Po{ii},opentaus{ii}, openareas{ii}, shuttaus{ii}, shutareas{ii});

end