% Creates the variable `fileName` in the base workspace
[d, s] = dwtread;
[pthstr, name] = fileparts(fileName);

% Log results
logfname = fullfile(pthstr, [name, '.log']);
logfile = fopen(logfname, 'w');
fprintf(logfile, 'Read dwell times from %s\n', name);

tres = 0.05;  % 50 microseconds
fprintf(logfile, 'Imposing %f millisecond open and shut time resolution\n', tres);
[rd, rs] = imposeres(d, s, tres, tres);

sampling_int = 0.025;
fprintf(logfile, 'Smoothing the pdf with sampling interval %f\n', sampling_int);
redist_rd = smoothbinpdf(rd, sampling_int, Inf);

fprintf(logfile, 'Plotting histograms\n-------------------------\n');
openbinwidth = 0.05;
shutbinwidth = 0.1;
fprintf(logfile, 'Open binwidth: %f\nShut binwidth: %f\n', openbinwidth, shutbinwidth);
[opencts, shutcts, xopen, xshut] = plothist(redist_rd, rs, openbinwidth, shutbinwidth);

% fprintf('Fitting open time histogram\n');
% [ opentaus, openareas, ~, Nopens, ~, ~, ~, ~, ~, ~, ~, ~, ...
%     pdfxlog, flog, ~, ~] = emdistfit(rd(rs == 1), [0.1 1], [0.3 0.7]);
fprintf(logfile, 'In a heretical way, we are fitting the redistributed (i.e. smoothed) histogram\n');
[ opentaus, openareas, ~, Nopens, ~, ~, ~, ~, ~, ~, ~, ~, ...
    pdfxlog, flog, ~, ~] = emdistfit(redist_rd(rs == 1), [0.1 1], [0.3 0.7]);
fprintf(logfile, 'Fit open times: %f\n', opentaus);
fprintf(logfile, 'Fit open areas: %f\n', openareas);

AUC_open = openbinwidth*sum(opencts);
% AUC_open = 0.05*sum(opencts); % 0.05 binwidth from plothist
opencts = opencts / AUC_open;
AUC_shut = shutbinwidth*sum(shutcts);
shutcts = shutcts / AUC_shut;

opencts_out = fullfile(pthstr, [name, '.opens']);
fprintf('Saving smooth open counts in %s\n', opencts_out);
csvwrite(opencts_out, [xopen', opencts']);

openfit_out = fullfile(pthstr, [name, '.ofit']);
fprintf('Saving exp fit for plotting to %s\n', openfit_out);
% scale pdf for the effective number of openings
openpdf = Nopens*openbinwidth*flog;
% now scale pdf to match the 'normalize' histogram
openpdf = openpdf / AUC_open;
csvwrite(openfit_out, [pdfxlog', openpdf']);

% Do not save shut times because
% 1 - they are meaningless for 2A(P552R) recordings
% 2 - need to verify bin width in plothist for shut times
% csvwrite(fullfile(pthstr, [name, '.shuts']), [xshut', shutcts']);

fclose(logfile);