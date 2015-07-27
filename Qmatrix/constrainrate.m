function [ A, B ] = constrainrate (q, idxAll, type, sourceIdx, targetIdx , c)
%CONSTRAINRATE This function returns a matrix, A, of coefficients giving the
%linear constraints on specified rate constants in the Q matrix and the
%equivalence column vector, B, such that A*theta = B where theta is a
%column vector of the rate constants
%   q - the Q matrix
%   idxAll - indices of all the rate constants that make up Q
%   type - the type of constraint, either 'fix', 'constrain', or 'loop'
%   sourceIdx - the indices of the source rate constants.  These are the
%        rate constants that are free to vary.  When type is 'fix', sourceIdx is
%        the indices of the rates that should be fixed to a certain value
%   targetIdx - indices of rate constants that will be constrained
%   c (optional) - vector of the constants for the constraints.  CONSTRAINRATE will
%       return the log10(c) in B.  If c is not given, it is assumed to be 1

if nargin < 6
    c=ones(1,length(sourceIdx));
end
if numel(c) == 1
    c = repmat(c,numel(sourceIdx),1);
end

switch type
    case 'fix'
        A = zeros(length(sourceIdx),length(idxAll));
        B = zeros(length(sourceIdx),1);
%         idx2 = interp1(idxAll, 1:length(idxAll), sourceIdx);
        [~,idx2] = ismember(sourceIdx, idxAll);
        for ii=1:length(sourceIdx)
            A(ii,idx2(ii)) = 1;
            B(ii) = log10(c(ii));
        end
    case 'constrain'
        A = zeros(length(sourceIdx),length(idxAll));
        B = zeros(length(sourceIdx),1);
%         idx2 = interp1(idxAll, 1:length(idxAll), sourceIdx);
%         idx3 = interp1(idxAll, 1:length(idxAll), targetIdx);
        [~,idx2] = ismember(sourceIdx, idxAll);
        [~,idx3] = ismember(targetIdx, idxAll);
        for ii=1:length(sourceIdx)
            A(ii, idx2(ii)) = -1;
            A(ii, idx3(ii)) = 1;
            B(ii) = log10(c(ii));
        end
    case 'loop'
        A = zeros(1,length(idxAll));
%         idx2 = interp1(idxAll,1:length(idxAll),[sourceIdx, targetIdx]);
        [~,idx2] = ismember([sourceIdx,targetIdx], idxAll);
        A(idx2) = [ones(1,length(sourceIdx)), -1*ones(1,length(targetIdx))];
        B = log10(c(1));
end

end

