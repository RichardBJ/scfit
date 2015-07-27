function [q, A, F, idxall, idxvary, gamma, xi, nstates, nopenstates, ...
    nshutstates, numconstraints, idxconstrain, fname] = qmatxlsread( varargin )
%QMATXLSREAD Read in single channel model from an Excel file

p = inputParser;
addOptional(p,'filename','',@(x) ischar(x));
parse(p,varargin{:});
if isempty(p.Results.filename) || ~exist(p.Results.filename,'file')
    [file,path]=uigetfile({'*.xls;*.xlsx','Qmatrix Excel file'});
    if file==0
        error('No file selected');
    end
    fname = [path, file];
else
    fname = p.Results.filename;
end

tmp = xlsread(fname,'B1:B3');
nstates=tmp(1);
nopenstates=tmp(2);
nshutstates=tmp(3);
xlRange = sprintf('B4:%c4',char(uint8('B')+nopenstates-1));
A = xlsread(fname,xlRange);
xlRange = sprintf('B5:%c5',char(uint8('B')+nshutstates-1));
F = xlsread(fname,xlRange);
numconstraints = xlsread(fname,'B6:B6');
xlRange = sprintf('C10:%c%d',char(uint8('C')+nstates-1),10+nstates-1);
q = xlsread(fname,xlRange);
q(isnan(q))=0;
xlRange = sprintf('A%d:%c%d',10+nstates+1,char(uint8('A')+nstates+1),10+nstates+1+numconstraints-1);
[num,txt] = xlsread(fname,xlRange);

idxall = find(q);

gamma=[];
xi=[];
idxconstrain=zeros(1,numconstraints);
for ii=1:numconstraints
    switch txt{ii}
        case 'constrain'
            src = sub2ind(size(q),num(ii,1),num(ii,2));
            tgt = sub2ind(size(q),num(ii,3),num(ii,4));
            val = num(ii,5);
            [tmpg, tmpxi] = constrainrate(q,idxall,'constrain',src,tgt,val);
            gamma = [gamma;tmpg];
            xi = [xi;tmpxi];
            idxconstrain(ii)=tgt;
        case 'fix'
            src = sub2ind(size(q),num(ii,1),num(ii,2));
            val = num(ii,3);
            [tmpg, tmpxi] = constrainrate(q,idxall,'fix',src,[],val);
            gamma = [gamma;tmpg];
            xi = [xi;tmpxi];
            idxconstrain(ii)=src;
        case 'loop'
            nloopstates = sum(~isnan(num(ii,:)));
            src = zeros(1,nloopstates);
            tgt = zeros(1,nloopstates);
            for mm=1:nloopstates
                if mm==nloopstates
                    jj=1;
                else
                    jj=mm+1;
                end
                if mm==1
                    kk=nloopstates;
                else
                    kk=mm-1;
                end
                src(mm) = sub2ind(size(q),num(ii,mm),num(ii,jj));
                tgt(mm) = sub2ind(size(q),num(ii,mm),num(ii,kk));
            end
            if ~isempty(setdiff(tgt,idxconstrain))
                idxconstrain(ii) = tgt(find(setdiff(tgt,idxconstrain),1));
            elseif ~isempty(setdiff(src,idxconstrain))
                idxconstrain(ii) = src(find(setdiff(src,idxconstrain),1));
            else
                break;
            end
            [tmpg,tmpxi] = constrainrate(q, idxall, 'loop', src, tgt);
            gamma = [gamma;tmpg];
            xi = [xi;tmpxi];
    end
end

idxvary = setdiff(idxall,idxconstrain);

end

