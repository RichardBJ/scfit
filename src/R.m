function [ f ] = R( C, lambda, tau, s, areaR, mMax, t )
%R Calculates the value of the matrix function R(t)
%   R(t) is a kind of reliability or survivor function, in which R(i,j)
%   gives the probability that a resolved open time, starting in state i,
%   has not yet finished and is currently in state j.  Another way to state
%   this is that R(i,j)[t] = the probability that 1)you're in state j and 
%   2)there has been no resolvable shut time during the interval 0 to t
%   given that you were in state i at time 0.
%
%   For details see Hawkes et al (1990, 1992)
%
%   C - cell array returned from the function Cvals. C is used to calculate
%   the exact value of R for times less than or equal to twice the imposed
%   resolution
%
%   lambda - eigenvalues of the Q matrix. Returned from the functions Cvals
%   or qmatvals
%
%   tau - the imposed resolution
%
%   s - generalized eigenvalues returned from the function asymptoticRvals
%
%   areaR - matrix returned from asymptoticRvals that has the area of each
%   exponential component given in s
%
%   t - the time at which to return R(t)

tol = 1e-12;

nr = size(areaR,1);
nc = size(areaR,2);
f=zeros(nr,nc);

if t<0
    return;
end

if t<=mMax*tau
    kA = length(lambda);
    m = ceil(t./tau)-1;
    if abs(t) < tol
        m = 0;
    end
%     m = floor(t./tau);

    for ii = 1:(m+1)
        for jj = 1:kA
            for kk = 1:ii
                % Beware of the indexing!!!! KKO 140923
                f = f + real( (-1).^(ii-1) * C(:,:,jj,ii,kk)*(t-(ii-1)*tau).^(kk-1) * exp(-lambda(jj)*(t-(ii-1)*tau)) );
            end
        end
    end
else
    kA = length(s);
    for ii=1:kA
        f = f + real(exp(s(ii)*t)*areaR(:,:,ii));
    end
end

end

