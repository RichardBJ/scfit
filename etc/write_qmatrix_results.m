function write_qmatrix_results( qmffile, idealizedfile, ...
    qfit, A, F, idxall, idxvary, td, tcrit, numdwells, runtime, ll, ...
    eq_Po, opentaus, openareas, shuttaus, shutareas)
%WRITE_QMATRIX_RESULTS Write a comma delimited file with results of
%maximumum likelihood fitting of rate constants to idealized dwells

dateNum = now();
[~,qmfname,~] = fileparts(qmffile);
[~,idlname,~] = fileparts(idealizedfile);
outfilename = sprintf('%s - %s - %s.csv',qmfname,idlname,datestr(dateNum,'yymmdd HH_MM_SS'));

nstates = length(qfit);
assert(nstates == (numel(A)+numel(F)));
nrates = numel(idxall);
% nfree = numel(idxvary);
% nconstrain = nstates-nfree;
[state1,state2] = ind2sub(size(qfit),idxall);

fid = fopen(outfilename,'w');
fprintf(fid,'Results from Q matrix maximum likelihood fitting\n');
fprintf(fid,'Exported on %s\n\n',datestr(dateNum));
fprintf(fid,'Model file: %s\n', qmffile);
fprintf(fid,'Idealized data file: %s\n', idealizedfile);
fprintf(fid,'Dead Time (ms): %f\n',td);
fprintf(fid,'Tcrit (ms): %f\n',tcrit);
fprintf(fid,'Number of Fitted Dwells: %d\n',numdwells);
fprintf(fid,'Runtime (s): %.2f\n\n',runtime);
fprintf(fid,'QMATRIX -- rates in ms-1\n,');
fprintf(fid,'State%d,',1:nstates);
for ii=1:nstates
    fprintf(fid,'\nState%d,',ii);
    fprintf(fid,'%.6e,',qfit(ii,:));
end

fprintf(fid,'\n\nLog Likelihood, %.4e\n\n',ll);

fprintf(fid,'Transition,Rate,Constrain/Free\n');
for ii=1:nrates
    if qfit(idxall(ii))<0.1 || qfit(idxall(ii))>=1e4
        formatSpec = '%d-->%d,%.4e,';
    else
        formatSpec = '%d-->%d,%.3f,';
    end
    fprintf(fid,formatSpec,state1(ii),state2(ii),qfit(idxall(ii)));
    if ismember(idxall(ii),idxvary) %any(idxall(ii)==idxvary)
        fprintf(fid,'Free\n');
    else
        fprintf(fid,'Constrained\n');
    end
end

fprintf(fid,'\n\nState Number,Open/Shut,Equilibrium Occupancy\n');
for ii=1:nstates
    fprintf(fid,'%d,',ii);
    if ismember(ii,A)
        fprintf(fid,'Open,');
    elseif ismember(ii,F)
        fprintf(fid,'Shut,');
    end
    fprintf(fid,'%.4e\n',eq_Po(ii));
end

fprintf(fid,'\n\nOpen Taus (in ms)');
fprintf(fid,',%.4e',opentaus);
fprintf(fid,'\nOpen Areas');
fprintf(fid,',%.4e',openareas);
fprintf(fid,'\nShut Taus (in ms)');
fprintf(fid,',%.4e',shuttaus);
fprintf(fid,'\nShut Areas');
fprintf(fid,',%.4e',shutareas);

fclose(fid);
end