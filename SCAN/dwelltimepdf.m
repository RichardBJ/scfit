function [output] = dwelltimepdf(x, dt, k)
%DWELLTIMEPDF details
%   More details

output = zeros(size(x));

idx = (abs(x - k*dt) <= dt);
output(idx) = 1 - (1/dt)*abs(x(idx) - k*dt);
output = output * 40;

end