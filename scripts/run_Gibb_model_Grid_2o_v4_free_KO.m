% Gibb_model_Grid_2o_v4
% Test_Coop_v10a_MR_KO_works_AJG_KKO
% q(idxvary) = exprnd(2,size(idxvary));
% q = dlmread('starting rates 2.csv');
KKO_node_model_140922
% q(idxvary) = exprnd(0.75,size(idxvary));

% Find rates to constrain or vary
% This code works, but it might be better to use reduced row echelon form
% [U,R,e] = qr(gamma,'vector');
% idxvary = idxall(e(end-12:end));
% 
% allrates = getconrates(q,gamma,xi,idxall,idxvary);
% q(idxall) = allrates;
% 
% clear U R e idxvary2 allrates;

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

rates = cell(length(file),1);
ll = zeros(length(file),1);
qnew = cell(length(file),1);
runtime = zeros(length(file),1);
nDwells = zeros(length(file),1);
tcrit = zeros(length(file),1);
eq_Po = cell(length(file),1);
opentaus = cell(length(file),1);
openareas = cell(length(file),1);
shuttaus = cell(length(file),1);
shutareas = cell(length(file),1);

%% Run loop

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
% shuts = rd(rs==0);
% [t,w,ll] = emdistfit(shuts,[0.1 1 100],(1/3)*[1 1 1]);
% tcrit(ii) = findtcrit(t,w);

%% Split the activity into bursts
% tcrit(ii)=35;
% [bursts, bstates]=imposetcrit(rd,rs,tcrit(ii));
% bursts(:,all(isnan(bursts)))=[];
% bstates(:,all(isnan(bstates)))=[];
% nDwells(ii) = sum(sum(~isnan(bursts)));

%% Fit data
tic;

% Fit without tcrit
[rates{ii}, ll(ii), qnew{ii}, hessian, covarmat, corrmat, history] = ...
    hjcfit(rd,q,A,F,td,idxall,idxvary,gamma,xi,[],[],[]);
    
% Fit bursts separated by at least tcrit
% [rates{ii}, ll(ii), qnew{ii}, hessian, covarmat, corrmat] = ...
%     hjcfit(bursts,q,A,F,td,idxall,idxvary,gamma,xi,[],[],tcrit(ii));

runtime(ii) = toc;

%% Save results
eq_Po{ii} = eqoccupy(qnew{ii});
[opentaus{ii}, openareas{ii}, shuttaus{ii}, shutareas{ii}] = ...
    qmatopenshutpdf(qnew{ii},A,F);
write_qmatrix_results('KKO_node_model_140920',[path file{ii}],qnew{ii},A,F,idxall,idxvary,...
    td,tcrit(ii),nDwells(ii),runtime(ii),ll(ii),eq_Po{ii},opentaus{ii}, openareas{ii}, shuttaus{ii}, shutareas{ii});

end