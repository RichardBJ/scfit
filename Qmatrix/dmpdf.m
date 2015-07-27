function [ f ] = dmpdf( P, A, F, t )
%DMPDF probability density function for dwell times of a discrete Markov
%chain
%   P - transition probability matrix
%   A - vector identifying the open states
%   F - vector identifying the shut states
%   t - values at which to calculate the pdf

peq = dmequil(P);
uA = ones(length(A),1);
uF = ones(length(F),1);
P_AA = P(A,A);
P_AF = P(A,F);
n = max(t);
% [tsrt, idx] = sort(t,'descend');
% pie is the probability of being in an open state given you observe one
% closed interval
pie = peq(F)*P(F,A) ./ (peq(F)*P(F,A)*uA);

f = zeros(size(t));
% tmp = zeros(tsrt(1),1);
c = zeros(1,n-1);
% alpha = zeros(length(A),1,tsrt(1));
alpha = zeros(n-1,length(A));

if n>1
    alpha(1,:) = pie*P_AA;
    c(1) = alpha(1,:)*uA;
    alpha(1,:) = alpha(1,:)./c(1);
    for jj=2:n-1
        alpha(jj,:) = alpha(jj-1,:)*P_AA;
        c(jj) = alpha(jj,:)*uA;
        alpha(jj,:) = alpha(jj,:)./c(jj);
    end
end

for ii=1:length(t)
    if t(ii)>1
        f(ii) = prod([c(1:t(ii)-1), alpha(t(ii)-1,:)*P_AF*uF]);
    else
        f(ii) = pie*P_AF*uF;
    end
end

% f(ii) = pie*(P(A,A)^(t(ii)-1))*P(A,F)*ones(length(F),1);

% if size(t) ~= size(f)
%     f=f';
% end

end

