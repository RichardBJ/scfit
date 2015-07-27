function [Gamma, Xi, idxConstrain, idxMR, MST] = mrconstraints(q, idxall, idxconstrain, idxmr, gamma, xi)
%MRCONSTRAINTS Find rates to fix for microscopic reversibility in a gating
%mechanism represented by a Q matrix
%
% 1) Create an undirected graph of Q (with the physical constrains from
%idxconstrain -- this may not be needed)
%   a) Need to convert a directed graph to an undirected graph
%       i) Matlab in the example in graphisomorphism uses G2 = G1 + G1';
%       ii) One issue is that the rates from idxconstrain will be in either
%       the upper or lower triangle of the Q matrix, but for an undirected
%       graph, the Matlab algorithms ignore the upper triangle, so we'll
%       need to be sure to set the weights for both transitions in the edge
%       containing the rates in idxconstrain
%
% 2) Assign weights to the edges in the graph G (which is really the
%constrained Q matrix [undirected])
%
%   a) Do so in such a way that
%       weight = 1 for edges (i.e. connections) to be included in the
%       minimum spanning tree. These connections will be set by physical
%       constraints (and not microscopic reversibility, if possible)
%
%       weights = s where s is the # of states for edges that the user 
%       wants to be set by microscopic reversibility (and therefore these
%       are to be excluded from the minimum spanning tree)
%
%       weights = 2 for any other edge
%
% 3) Find the minimum spanning tree and from this select rates to constrain
% by microscopic reversibility (the weights from step 2a) should make it so
% that rates the user wants set by MR are excluded from the MST if possible
%
% 4) Find the unique shortest path (I assume this is guaranteed to exist in
% a minimum spanning tree) between each pair of states connected by an edge
% that is not in the MST
%   a) These paths will determine which cycle to use for the constrain

if ~isempty(idxconstrain) && ~isempty(gamma)
    [nrows, ncols] = size(gamma);
    gammaRank = rank(gamma);
    if nrows > gammaRank
        warning('The constraints on the rates in Q are not all independent');
    end
    idxConstrain = idxconstrain;
    Gamma = gamma;
    Xi = xi;
else
    idxConstrain = [];
    Gamma = [];
    Xi = [];
end

assert(isempty(intersect(idxconstrain,idxmr)), ...
    'Some constraints are set by physical constraints and microscopic reversibility');

G = zeros(size(q));
G(q~=0) = 2;
G(idxconstrain) = 1;
G = min(G,G');
% Assigning weights to connections the user wants fixed by microscopic
% reversibility after the weights of connections fixed by physical
% constraints could lead to non-independent constraints if some of the
% physical constraints are on the connections fixed by MR. The weights of
% rates fixed by MR could be set first to avoid this.
G(idxmr) = length(q);
G = max(G,G');

% A = q~=0;
% A(idxconstrain) = false;
% % This may not work if both rates in a connection are constrained, e.g.
% % rate 1->2 and rate 2->1 are constrained
% A = xor(A,A'); %Will not work if both rates in a connection are fixed by physical constraints
% G(A) = 1;
% B = q~=0;
% B(idxmr) = false;
% B = xor(B,B');
% % G(B&~A) = length(q);
% G(B) = length(q);
% % What if the user wants rate r->s to be set by physical constraints but
% % rate s->r to be set by microscopic reversibility? It may be possible to
% % do this, but it could lead to non-independent constraints. It's probably
% % best to ensure the rates set by physical constraints will be included in
% % the minimum spanning tree so they can be set later.

MST = graphminspantree(sparse(G));

% Let's just use the transitions in the lower triangular part of Q to fix for MR
% idxMR = setdiff(find(tril(q)), find(MST));
% [iiMR,jjMR] = find(q & ~MST);
idxMR = find(tril(q) & ~MST);

for kk = 1:numel(idxMR)
    [ii,jj] = ind2sub(size(q),idxMR(kk));
    % Check that the user did not specify a rate to fix for microscopic
    % reversibility that is in the upper triangular
    if ismember(sub2ind(size(q),jj,ii),idxmr)
        tmp = ii;
        ii = jj;
        jj = tmp;
        idxMR(kk) = sub2ind(size(q),ii,jj);
    end
    [~, pth] = graphshortestpath(MST,ii,jj,'Directed',false);
    src = sub2ind(size(q), pth, circshift(pth,[0 -1]));
    tgt = sub2ind(size(q), fliplr(pth), circshift(fliplr(pth),[0 -1]));
    [tmpg,tmpxi] = constrainrate(q,idxall,'loop',src,tgt);
    Gamma=[Gamma;tmpg];
    Xi=[Xi;tmpxi];
    idxConstrain = [idxConstrain; idxMR(kk)];
end

% Now we need to set the rates to comply with constraints
% This may not need to be done here, unless desired, but will be important
% for the fitting algorithm

% We can probably ignore the "extra", non-independent rates when imposing
% constraints by solving the underdetermined system for the minimum norm
% solution (i.e. find [one of the infinitely many sets of] rates satisfying
% the constraints). For example,
% [Q,R] = qr(gamma');
% [~,R1] = qr(gamma',0); %Does Q not need to be partitioned as well? - No,
% I don't think so because the QR decomposition of the transpose of gamma
% will produce in R an upper triangluar matrix, so the first r rows of Q
% (where r = rank(Gamma);) will be the ones corresponding to R1.
% x = Q*[R1'\xi;zeros(length(Q)-length(R1),1)];
% x is now the logarithm of the set of rates that will satisfy the
% constraints

end

