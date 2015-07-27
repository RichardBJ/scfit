function [taus, lambda, A, areas] = qmatvals( q, p0 )
%QMATVALS Calculates the time constants and spectral matrices (which are
%used to calculate amplitudes) for a given Q matrix 
%   q: Q matrix, a k-by-k matrix where k is the number of states and 
%       q(i,j) is the transition rate from state i to state j and
%       q(i,i) is a number such that the sum of each row is zero
%   OUTPUT
%       taus - time constant for macroscopic fluctuations
%       lambda - eigenvalues of q; lambda(i) corresponds to A(:,:,i)
%       A - spectral matrices of q
%       areas - if a p0 is supplied, the areas of each tau is returned

tol = 1e-12;

[v, lambda] = eig(-q);
y = inv(v);
A = zeros(length(v),length(v),size(q,1));
for k=1:size(q,1)
    A(:,:,k) = real(v(:,k)*y(k,:));
end
lambda = diag(lambda);
lambda(abs(lambda)<tol)=0;
[lambda, idx] = sort(lambda);
taus = 1./lambda;
A = A(:,:,idx);

%if an eigenvalue of Q is zero, then there is not a corresponding tau
taus(lambda==0)=[];

if nargin==2
    if (size(p0,2) == 1)
        p0=p0';
    end
    areas = zeros(size(q));
    for aa = 1:size(q,1)
        for bb = 2:size(q,1)
            areas(aa,bb) = p0*A(:,aa,bb);
        end
    end
    areas(:,1)=[];
end

end

