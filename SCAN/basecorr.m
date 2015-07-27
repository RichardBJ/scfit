function data = basecorr (data, adint, npts, obj)
%BASECORR   subtracts baseline from singel channel data
%

if isempty(npts) npts = 20*40000; end
% npts = length(data);
bdata = zeros(npts,1);
time = adint*(1:npts)';
warning ('off','all');
last = ceil(length(data)/npts);
% h = gaussfir(0.01*sqrt(2),50,1);
% h = gaussfir(0.01,50,1);
h = ones(1,50)./50;
% s.mu=[0; 20];
% s.Sigma=cat(3,1,1);
% s.PComponents=[0.5, 0.5];
for ii=1:last
    openamp=0;
    shutamp=0;
    dataf=0;
    datafx=0;
    ind = [1:npts] + (ii-1)*npts;
    tf = (1:npts)';
    if ind(end)>length(data)
        ind = ind(1):length(data);
        tf = (1:length(ind))';
    end
    bdata(tf) = data(ind);
%     p = polyfit(tf,bdata(tf),1);
%     bdata(tf) = bdata(tf)-polyval(p,tf);
    try
%         gmfit = gmdistribution.fit(bdata(1:length(ind)),2,'Replicates',3);
        if isempty(obj)
            gmfit = gmdistribution.fit(bdata(tf),2);
        else
            gmfit = obj;
        end
        [shutamp,shutstate] = min(gmfit.mu);
        [openamp,openstate] = max(gmfit.mu);
        if openamp-shutamp < 3
%             data(ind) = bdata(1:length(ind))-gaussfilter(bdata(1:length(ind)),0.01);
            try
                gmfit = gmdistribution.fit(bdata(tf),2,'Replicates',3);
                [shutamp,shutstate] = min(gmfit.mu);
                [openamp,openstate] = max(gmfit.mu);
            catch error
                data(ind) = bdata(tf)-filtfilt(h,1,bdata(tf));
%                 data(ind) = bdata(1:length(ind))-filter(h,1,bdata(1:length(ind)));
%                 data(ind) = bdata(tf)-gaussfilter(bdata(tf),0.01);
            end
        end
        if openamp-shutamp < 3
            data(ind) = bdata(tf)-filtfilt(h,1,bdata(tf));
%             data(ind) = bdata(1:length(ind))-filter(h,1,bdata(1:length(ind)));
%             data(ind) = bdata(tf)-gaussfilter(bdata(tf),0.01);
        else
            idx = cluster(gmfit,bdata(tf));
%             dataf = gaussfilter(bdata(idx==shutstate),0.01);
            dataf = filtfilt(h,1,bdata(idx==shutstate));
%             dataf = filter(h,1,bdata(idx==shutstate));
            datafx = interp1(time(idx==shutstate),dataf,time(idx==openstate));
            idxnan = find(isnan(datafx));
            if ~isempty(idxnan) 
                datafx(idxnan) = interp1(time(idx==shutstate),dataf,time(idxnan),'nearest','extrap');
            end
%                 if find(~isnan(datafx),1) ~= 1
%                     datafx(idxnan) = datafx(idxnan(end)+1);
%                 else
%                     datafx(idxnan) = datafx(idxnan(1)-1);
%                 end
%             end
%             plot(time,bdata,time(idx==shutstate),dataf,'.',time(idx==openstate),datafx,'.');
            bdata(idx==shutstate) = bdata(idx==shutstate)-dataf;
            bdata(idx==openstate) = bdata(idx==openstate)-datafx;
            data(ind)=bdata(tf);
        end
    catch error
%         data(ind) = bdata(tf)-gaussfilter(bdata(tf),0.01);
        data(ind) = bdata(tf)-filtfilt(h,1,bdata(tf));
%         data(ind) = bdata(1:length(ind))-filter(h,1,bdata(1:length(ind)));
        ii
    end
    if ~mod(ii,25)
        ii
    end
end

% try
%     gmfit = gmdistribution.fit(data(1:10*npts),2,'Replicates',3);
%     s.mu=gmfit.mu;
%     s.Sigma=gmfit.Sigma;
%     s.PComponents=gmfit.PComponents;
% catch error
%     s.mu=[0, 8];
%     s.Sigma=cat(3,0.1,0.3);
%     s.PComponents=[0.5, 0.5];
% end
% last = ceil(length(data)/npts)
% % last = 1000;
% for ii=1:last
%     ind = [1:npts] + (ii-1)*npts;
%     if ind(end)>length(data)
%         ind = ind(1):length(data);
%     end
%     try 
%         gmfit = gmdistribution.fit(data(ind),2,'Start',s);
%     catch err
%         gmfit = gmdistribution.fit(data(ind),1);
%     end
%     data(ind)=data(ind)-min(gmfit.mu);
%     if ~mod(ii,100)
%         disp(ii);
%     end
% %     ko = cluster (gmfit,abf(ind));
% %     plot(adint*ind(ko==1),abf(ind(ko==1)),'.b',adint*ind(ko==2),abf(ind(ko==2)),'.g');
% %     hold on;
% %     mu(ii,:) = gmfit.mu;
% end

warning('on','all');
