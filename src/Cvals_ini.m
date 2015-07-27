function [ out ] = Cvals_ini( A, D, ii, m )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% if m==0
%     out = A(:,:,ii);
% else
%     out = D(:,:,ii)*Cvals_ini(A, D, ii, m-1) ./ m;
% end

out = A(:,:,ii);

for kk=1:m
    out = D(:,:,ii)*out ./ kk;
end

end

