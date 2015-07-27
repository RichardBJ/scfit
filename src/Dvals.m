function [ D, taus, lambda, specA] = Dvals( q, A, F, td)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

qAA = q(A,A);
qFF = q(F,F);
qAF = q(A,F);
qFA = q(F,A);
[taus, lambda, specA] = qmatvals (q);
% specA2 = mat2cell(specA(:,:),length(q),length(q)*ones(1,length(q)));
eQFFt = expm(qFF*td);
D = zeros(length(A),length(A),length(q));
% D2 = cell(length(q),1);%zeros([size(qAA),length(q)]);
for ii=1:length(q)
    D(:,:,ii) = specA(A,F,ii)*eQFFt*qFA;
%     D2{ii} = specA2{ii}(A,F)*eQFFt*qFA;
end

end

