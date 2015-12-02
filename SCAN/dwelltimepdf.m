function [output] = dwelltimepdf(x, dt, k)
%DWELLTIMEPDF Peak pdf at k*dt with linear fall off to (k-1)dt and (k+1)dt
%   From QuB Manual: We observed that events with sampled
%   duration kdt have true duration distributed between (k - 1)dt and
%   (k + 1)dt with a peak at kdt and linear fall-off. The "smooth binning"
%   option performs this re-distribution.

output = zeros(size(x));

idx = (abs(x - k*dt) <= dt);
output(idx) = 1 - (1/dt)*abs(x(idx) - k*dt);

% Since this is a pdf, the area under the curve must be one
% Because we made the peak == 1, and we assume dt (the sampling interval)
% to be 0.025 (i.e. sampled at 40 kHz) our normalization factor is 40
output = output * 40;

end