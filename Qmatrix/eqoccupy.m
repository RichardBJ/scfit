function peq = eqoccupy (q)
%EQOCCUPY Calculates the equilibrium state occupancies for the Q matrix
%   q: a k-by-k matrix where k is the number of states and 
%       q(i,j) is the transition rate from state i to state j and
%       q(i,i) is a number such that the sum of each row is zero
%       -1/q(i,i) also happens to be the mean lifetime of a sojourn 
%       in the state i

n=length(q);

% At equilibrium, the transition probabilities do not change with time
% Hence, dp/dt = 0, but dp/dt = p(t)*Q so
% p(inf)*Q = 0, subject to the sum of the probabilities must be 1, i.e.
% p(inf)*u = 1, where u is a column vector of ones
%
% Hence p(inf)*(Q|u) = (0|1)
% or p(inf) = (0|1) / (Q|u)

peq = [zeros(1,n), 1] / [q, ones(n,1)];

% The following code is from Colquhoun and Hawkes 1995 - works fine
% u = ones(1,size(q,1));
% S = [q, ones(size(q,1),1)];
% % Test for singularity
% M = S*S';
% if rank(M) < min(size(M))
%     %singluar
%     %so calculate the states at some long time at which equilibrium
%     %is likely to have been reached (assume time is in seconds)
% %     peq = u*expm(q*100);
% %     fprintf ('Singular');
% %     try
% %         peq = u/M;
% %     catch
%         warning ('Matrix is Singular');
%         q;
%         %well, output something, I guess
%         peq = zeros(1,size(q,1));
%         peq(1) = 1;
% %     end
% else
% %     peq = u/(S*S');
%     peq = u/M;
% end

end