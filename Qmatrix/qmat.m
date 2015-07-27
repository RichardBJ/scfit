function pt = qmat(q, t, p0)
%QMAT Calculates the occupancies of each state at time t given by the
%Q matrix, q
%   q: Q matrix, a k-by-k matrix where k is the number of states and 
%       q(i,j) is the transition rate from state i to state j and
%       q(i,i) is a number such that the sum of each row is zero
%   t: a vector of times
%
%   pt: a length(t)-by-k vector of occupancies of state k at time t
%   p0 (optional): starting occupancies

pt = zeros(length(t),size(q,1));
if nargin < 3
    p0 = eqoccupy(q);
end

for k=1:length(t)
    pt(k,:) = p0*expm(q*t(k));
end

end

