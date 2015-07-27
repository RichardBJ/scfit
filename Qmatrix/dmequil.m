function [ pie ] = dmequil( P )
%DMEQUIL Calculates the equilibrium distribution for a discrete-time Markov
%chain
% Input
%   P - transition probability matrix of discrete-time Markov chain
% Output
%   pie - row vector of equilibrium distribution
% 
% For discrete-time Markov chains
%       pie = pie*P subject to sum(pie) = 1;

pie = [zeros(1,length(P)), 1] / [P-eye(size(P)), ones(length(P),1)];
end

