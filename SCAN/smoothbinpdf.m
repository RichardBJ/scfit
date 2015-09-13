function [ redist_dwells ] = smoothbinpdf( dwells, varargin )
%SMOOTHBINPDF Smooth single channel histograms from QuB dwt
%   REDIST_DWELLS = SMOOTHBINPDF(DWELLS) redistributes the dwells that are
%   less than or equal to 3 times the sampling interval, which is taken as
%   the minimum value in the DWELLS array.
%
%   REDIST_DWELLS = SMOOTHBINPDF(DWELLS, DT) specificies the sampling
%   interval, DT.
%
%   REDIST_DWELLS = SMOOTHBINPDF(DWELLS, DT, THLD) also redistributes
%   dwells less than or equal to THLD times DT. You can pass [] for DT to
%   use the default value.
%
%   From the QuB manual (www.qub.buffalo.edu/wiki/index.php/Modeling:MIL
%
%   This basic binning algorithm gives misleading output for the shortest
%   events, since events appear as exact multiples of the sampling time.
%   Some of the smallest bins contain no exact multiples of ?t, so events
%   with true duration in the bin are counted in neighboring bins. We
%   simulated data and recorded each event's true duration as well as
%   sampled duration, sampling at ?t. We observed that events with sampled
%   duration k?t have true duration distributed between (k ? 1)?t and
%   (k + 1)?t with a peak at k?t and linear fall-off. The "smooth binning"
%   option performs this re-distribution.
%
%   KKO: I am not sure how they simulated the data. Did they increase the 
%   adinterval for simulation and then reduce the number of points (by,
%   for example, only keeping every fifth point)?  Also, I'm not sure what
%   they mean by a linear fall off. Do they mean that it goes to zero at
%   (k-1)*dt and (k+1)*dt? Or perhaps, it jumps down to zero outside of the
%   range (k-1)*dt to (k+1)*dt. Moreover, what about correlations in their
%   simulation? For example, if a dwell was sampled at (k+1)*dt did that
%   affect the sample duration of the following dwell?

% calculate A-D interval
if length(varargin) < 1 || isempty(varargin{1})
    dt = min(dwells);
else
    dt = varargin{1};
end

% set threshold
if length(varargin) < 2
    thld = 3;
else
    thld = varargin{2};
end

% set a threshold below which to re-distribute
% this should depend on the binning method used
% for example, with log-scale binning, each bin spans a different time
% interval, so redistributing a dwell of k*dt to between (k-1)*dt and
% (k+1)*dt won't change the bin a dwell falls in above a certain threshold

redist_dwells = dwells;
dwells = round(dwells / dt);
idx = find(dwells <= thld);

for ii = 1:length(idx)
    k = dwells(idx(ii));
    % can k ever be one, yes, so make sure dead time is 0.05, say, and dt
    % is 0.025 (so you don't get values down to zero.
    f = @(z) dwelltimepdf(z, dt, k);
    g = @(z) normpdf(z, k*dt, dt);
    grnd = @() normrnd(k*dt, dt, 1);
    redist_dwells(idx(ii)) = accrejrnd(f, g, grnd, 2.7, 1, 1);

    % The following lines approximate the re-distribution pdf using a
    % normal distribution, instead of a linear fall off.
%     tmp_dwell = normrnd(k*dt, dt);
%     if abs(tmp_dwell - k*dt) <= dt
%         redist_dwells(idx(ii)) = tmp_dwell;
%     end
end

% TODO -- I think I did this already KKO 150913
% * use rejection sampling to create random samples from a distribution
% with a peak at k*dt and with linear fall off to (k-1)*dt and (k+1)*dt.
% For example, generate random normal number and a random uniform number
% (between 0 and 1). If the random normal number times a constant is less
% than the target distribution (peak at k*dt, etc.) then keep the sample,
% otherwise discard it.

% See MATLAB help on Acceptance-Rejection method for generating
% pseudorandom numbers.

end

