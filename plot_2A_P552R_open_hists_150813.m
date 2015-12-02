[d, s] = dwtread;
[rd, rs] = imposeres(d, s, 0.05, 0.05);
redist_rd = smoothbinpdf(rd, 0.05, Inf);
[opencts, shutcts, xopen, xshut] = plothist(redist_rd, rs);
emdistfit(rd(rs == 1), [0.1 1], [0.3 0.7])
% updated plothist sets 10 bins per decade, or a binwidth of 0.1
AUC_open = 0.1*sum(opencts);
% AUC_open = 0.05*sum(opencts); % 0.05 binwidth from plothist
opencts = opencts / AUC_open;
% for shut times, I think the binwidth in plothist is 0.1, not 0.05
AUC_shut = 0.05*sum(shutcts);
shutcts = shutcts / AUC_shut;
[pthstr, name] = fileparts(fileName);
csvwrite(fullfile(pthstr, [name, '.opens']), [xopen', opencts']);
csvwrite(fullfile(pthstr, [name, '.shuts']), [xshut', shutcts']);
