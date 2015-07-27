function [ allrates, conrates, idxconstrain] = getconrates( q, gamma, xi, idxAll, idxVary )
%GETCONRATES Get constrained rates from 
%   Detailed explanation goes here

    % Following (Qin et al. 1996 Biophys J) or (Golub and Van Loan 1989)
    % If the linear constraints on all the parameters mu, where mu(i,j) = log10(q(i,j)),
    % are represented in the matrix gamma and
    %                           gamma * theta = xi
    % where theta = (...mu(i,j)...)' and xi is a constant vector
    % Then theta can be formulated into linear combinations of unconstrained variables
    % using the QR factorization of gamma  
    %   assume the constraints are independent of one another, which means
    %   gamma has full rank
    
    mn = size(idxAll);
    if mn(1) == 1
        idxAll = idxAll.';
    end
    mn = size(idxVary);
    if mn(1) == 1
        idxVary = idxVary.';
    end
    
    [U,R1,R2,idxtheta,Gamma,R,idxconstrain,nConstrain] = constrainqr(gamma,idxAll,idxVary);
    M = [ -R1\R2; eye(numel(idxVary)) ];
    b = [ R1\(U\xi) ; zeros(numel(idxVary),1)];
    % now q(idxtheta) == 10.^(M*x + b); where x are the variables to be optimized
    theta2 = log10(q(idxVary));
    theta = M*theta2 + b;
    allrates = 10.^theta;
    conrates = allrates(1:nConstrain);
    
    [~,iThetaToAll] = ismember(idxAll,idxtheta);
    allrates = allrates(iThetaToAll);
    
%     theta1 = -R1\R2*theta2 + R1\(U\xi);
%     conrates = 10.^theta1;


end

