function [ y ] = rootsdetW( q,A,F,tau,x0 ) %#codegen
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% f = @(x) detW_mex(x,q,A,F,tau);
% try
    y = fzeroKO(@detW,x0,q,A,F,tau);
% catch error
    if isnan(y)
        y=inf;
    end
%     fprintf ('%s', error.message);
% end    
end

