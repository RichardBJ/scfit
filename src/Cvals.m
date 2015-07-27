function [ C, lambda ] = Cvals( q, A, F, td, mMax )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% [~, lambda, sA] = qmatvals(q);
[D, ~, lambda, sA] = Dvals (q,A,F,td);
nstates=length(q);
% C = cell(nstates,mMax,mMax);
C = zeros(length(A) ,length(A) , nstates, mMax+1, mMax+1);
sAaa = sA(A,A,:);
% sAaa = mat2cell(sAaa(:,:),length(A),length(A)*ones(1,nstates));

for m = 0:mMax
    for n = 0:m
        for ii=1:nstates
            tmp=zeros(length(A),length(A));
            if n == m
%                 C{ii,m+1,n+1} = Cvals_ini(sAaa, D, ii, m);
                C(:,:,ii,m+1,n+1) = Cvals_ini(sAaa, D, ii, m);
            elseif n == 0
                for jj = 1:nstates
                    if ii==jj
                        continue;
                    end
                    for r = 0:m-1
                        tmp = tmp + real(D(:,:,ii)*C(:,:,jj,m,r+1)*factorial(r) ./ (lambda(jj)-lambda(ii)).^(r+1)) ...
                              - real(D(:,:,jj)*C(:,:,ii,m,r+1)*factorial(r) ./ (lambda(ii)-lambda(jj)).^(r+1));
                    end
                end
                C(:,:,ii,m+1,n+1) = tmp;
            else
                for jj = 1:nstates
                    if ii==jj
                        continue;
                    end
                    for r = n:m-1
                        tmp = tmp + real(D(:,:,jj)*C(:,:,ii,m,r+1)*factorial(r) ./ (n .* (lambda(ii)-lambda(jj)).^(r-n+1)));
                    end
                end
                C(:,:,ii,m+1,n+1) = D(:,:,ii)*C(:,:,ii,m,n)./n - tmp;
            end
        end
    end
end
        
end

