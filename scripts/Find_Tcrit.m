%% Get dwt files
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

tcrit = zeros(8,length(file));

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
    
    shuts = rd(rs==0);
    for jj=3:6
        %% Calculate tcrit
        ini_tau = logspace(-1,4,jj);
        ini_weights = ones(1,jj)./jj;
        [t,w,ll] = emdistfit(shuts,ini_tau,ini_weights);
        tcrit(jj-2,ii) = findtcrit(t,w);
    end
end