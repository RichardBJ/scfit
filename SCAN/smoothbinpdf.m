function [ output_args ] = smoothbinpdf( input_args )
%SMOOTHBINPDF Smooth single channel histograms from QuB dwt
%   Detailed explanation goes here

% calculate A-D interval
adint = min(diff(dwells));

% set a threshold below which to re-distribute
% this should depend on the binning method used
% for example, with log-scale binning, each bin spans a different time
% interval, so redistributing a dwell of k*dt to between (k-1)*dt and
% (k+1)*dt won't change the bin a dwell falls in above a certain threshold

end

