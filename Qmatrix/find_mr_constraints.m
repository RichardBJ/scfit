function [ idx, ii, jj, transitions ] = find_mr_constraints( q )
%FIND_MR_CONSTRAINTS Finds rates in Q matrix to fix for microscopic
%revesibility
%   q - connectivity matrix representing single channel activation scheme.
%   The diagonal elements of q should be zero.
%
%   Output
%   idx - indices of rates in Q matrix to fix for microscopic reversibility
%       use ind2sub() to convert these to subscripts, e.g. Q(i,j)
%   ii, jj - subscripts of constrained rates to use in Q(ii,jj).
%       [ii,jj] = ind2sub(size(q),idx);
%   transitions - String desribing the transitions to be fixed for
%       microscopic reversibility, e.g. 1 -> 2

Q = sparse(q);
if ~graphisspantree(Q)
   Tree = graphminspantree(Q);
   idx = setdiff(find(tril(q)), find(Tree));
else
    idx = [];
end

[ii, jj] = ind2sub(size(q),idx);
transitions = num2str([ii,jj],'%d -> %d');

end

