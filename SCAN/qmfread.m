function [q, A, F, idxall, idxvary, gamma, xi, numstates, numrates, ...
    numconstraints, idxconstrain, root, filename] = qmfread(varargin)
%QMFREAD Read QuB Model File, qmf

p = inputParser;
addOptional(p,'filename','',@(x) ischar(x) || iscellstr(x));
parse(p,varargin{:});
if isempty(p.Results.filename) || ~exist(p.Results.filename,'file')
    [file,path]=uigetfile({'*.qmf','QuB Model File';'*.abf','Axon Binary'});
    if file==0
        root=0;
        filename=0;
        states=0;
        return
    end
    filename = [path file];
else
    filename = p.Results.filename;
end
    

fid = fopen(filename);
% magicStr = fread(fid,[1 12],'*char*1');
magicStr = fscanf(fid,'%12c',1);
if ~strcmp(magicStr,'QUB_(;-)_QFS')
    error('%s is not a valid QuB model file (qmf)',[path file]);
end
root = qtreadnode(fid);
root.child = qtgetchild(fid,root);
for ii=1:numel(root.child)
    root.child(ii).child = qtgetchild(fid,root.child(ii));
    for jj=1:numel(root.child(ii).child)
        root.child(ii).child(jj).child = qtgetchild(fid,root.child(ii).child(jj));
        for kk=1:numel(root.child(ii).child(jj).child)
            root.child(ii).child(jj).child(kk).child = qtgetchild(fid,root.child(ii).child(jj).child(kk));
        end
    end
end
idx = strcmp({root.child.name},'States');
numstates = numel(root.child(idx).child);
q = zeros(numstates);
for ii=1:numstates
    st = root.child(idx).child(ii).child;
    idx_cl = strcmp({st.name},'Class');
    st_class(ii) = st(idx_cl).data;
end
A = find(st_class==1);
F = find(st_class==0);
idx = strcmp({root.child.name},'Rates');
numrates = numel(root.child(idx).child);
for ii=1:numrates
    rt = root.child(idx).child(ii).child;
    idx_st = strcmp({rt.name},'States');
    idx_rt = strcmp({rt.name},'k0');
    states = rt(idx_st).data + 1;
    rates = rt(idx_rt).data;
    q(states(1),states(2)) = rates(1);
    q(states(2),states(1)) = rates(2);
end
idxall = find(q);
idx = strcmp({root.child.name},'Constraints');
numconstraints = numel(root.child(idx).child);
idxconstrain = zeros(1,numconstraints);
gamma=[];
xi=[];
for ii=1:numconstraints
    con = root.child(idx).child(ii);
    switch con.name
        case 'FixRate'
            states = con.data+1;
            src =  sub2ind(size(q),states(1),states(2));
            jdx = strcmp({con.child.name},'HasValue');
            if con.child(jdx).data == 1
                kdx = strcmp({con.child.name},'Value');
                val = con.child(kdx).data;
            else
                val = q(src);
            end
            [tmpg,tmpxi] = constrainrate(q, idxall, 'fix', src, [], val);
            gamma = [gamma;tmpg];
            xi = [xi;tmpxi];
            idxconstrain(ii) = src;
        case 'FixExp'
            %
        case 'ScaleRate'
            states = con.data+1;
            src = sub2ind(size(q),states(1),states(2));
            tgt = sub2ind(size(q),states(3),states(4));
            jdx = strcmp({con.child.name},'HasValue');
            if con.child(jdx).data == 1
                kdx = strcmp({con.child.name},'Value');
                val = 1 ./ con.child(kdx).data;
            else
                val = q(tgt) ./ q(src);
            end
            [tmpg,tmpxi] = constrainrate(q, idxall, 'constrain', src, tgt, val);
            gamma = [gamma;tmpg];
            xi = [xi;tmpxi];
            idxconstrain(ii) = tgt;
        case 'ScaleExp'
            %
        case 'LoopBal'
            states = con.data+1;
            src = zeros(1,numel(states));
            tgt = zeros(1,numel(states));
            for mm=1:numel(states)
                if mm==numel(states)
                    jj=1;
                else
                    jj=mm+1;
                end
                if mm==1
                    kk=numel(states);
                else
                    kk=mm-1;
                end
                src(mm) = sub2ind(size(q),states(mm),states(jj));
                tgt(mm) = sub2ind(size(q),states(mm),states(kk));
            end
            [tmpg,tmpxi] = constrainrate(q, idxall, 'loop', src, tgt);
            gamma = [gamma;tmpg];
            xi = [xi;tmpxi];
%             idxconstrain(ii) = tgt(find(setdiff(tgt,idxconstrain),1));
            idxconstrain(ii) = src(1);
        case 'LoopImbal'
    end
end
idxvary = setdiff(idxall,idxconstrain);
fclose(fid);
end

function child = qtgetchild(fid,node)
    if node.childpos == 0
       child=[];
    else
        fseek(fid,node.childpos,'bof');
        child = qtreadnode(fid);
        nodeIter = child;
        ii=2;
        while nodeIter.siblingpos ~= 0
            fseek(fid,nodeIter.siblingpos,'bof');
            child(ii) = qtreadnode(fid);
            nodeIter = child(ii);
            ii=ii+1;
        end
    end
end

function node = qtreadnode(fid)
    node.flags = fread(fid,3,'uchar');
    node.datatype = fread(fid,1,'uchar');
    node.datasize = fread(fid,1,'uint');
    node.datacount = fread(fid,1,'uint');
    node.datapos = fread(fid,1,'uint');
    node.childpos = fread(fid,1,'uint');
    node.siblingpos = fread(fid,1,'uint');
    fseek(fid,7,'cof');
    node.namelen = fread(fid,1,'uchar');
%     node.name = fread(fid,[1 node.namelen],'*char*1');
    node.name = fscanf(fid,'%c',[1 node.namelen]);
    node.data = qtreaddata(fid,node);
end

function data = qtreaddata(fid,node)
    if node.datatype~=0 && node.datapos~=0
        fseek(fid,node.datapos,'bof');
        switch node.datatype
            case 3
                dtype = '*char*1';
            case 4
                dtype = 'uchar';
            case 5
                dtype = 'schar';
            case 6
                dtype = 'ushort';
            case 7
                dtype = 'short';
            case 8
                dtype = 'uint';
            case 9
                dtype = 'int';
            case 12
                dtype = 'float';
            case 13
                dtype = 'double';
            otherwise
                dtype = sprintf('integer*%d',node.datasize);
        end
        if node.datatype==3
            data = fread(fid,[1 node.datacount],dtype);
        else
            data = fread(fid,node.datacount,dtype);
        end
    else
        data=[];
    end
end