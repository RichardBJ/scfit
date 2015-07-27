[files, pth] = uigetfile('*.dwt','multiselect','on');

ncomponents = 2;
opentau = zeros(length(files),ncomponents);
openweight = zeros(length(files),ncomponents);
nropens = zeros(length(files),1);
res = 0.05;

for ii=1:length(files)
    fprintf('Fitting file #%d - %s\n',ii,files{ii});
    [d,s] = dwtread([pth,files{ii}]);
    [rd,rs] = imposeres(d,s,res,res);
    [opentau(ii,:), openweight(ii,:), ~, ~, ~, ~, ~, ~, ~, ...
        ~,~,~,~,~,~,exitflag(ii)] = emdistfit(rd(rs==1),[0.1 10], [1 1]/2);
    title(files{ii},'interpreter','none');
    nropens(ii) = numel(rd(rs==1));
    
    % Print traces to bmp file
    set(gcf,'PaperSize',[11 8.5]);
    set(gcf,'PaperPosition',[0.5 0.5 10 7.5]);
    [~,name,~] = fileparts(files{ii});
%     print('-dbmp','-r300',[pth, name, '.bmp']);
%     BMP format caused very large files to be saved, even when a lower 
%     resolution was used
    print('-djpeg','-r300',[pth, name, '.jpg']);
    
%     close;    
end

[sort_ot,sidx] = sort(opentau,2);
for ii=1:length(files)
    sort_ow(ii,:) = openweight(ii,sidx(ii,:));
end