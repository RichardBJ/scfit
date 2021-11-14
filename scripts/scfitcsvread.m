function [durations,amps]=scfitcsvread(fname)
si=0.01;
data=readmatrix(fname,delimiter);
[m,n]=size(data);
if n==1
    disp("1 column, must just be events");
elseif n==2
    disp("2 columns, so could be time vs events or raw and ideal");
    %This bit ...coercing to a time column was completely unecessary!
    % ...given you just need the sample interval!!!
    data=coercetime(data,si);
elseif n==3
    disp("3 columns, must be time, raw, ideal");
    si=data(10,1)-data(9,1);
end
last=0;
dur=0;
loc=0;
durs=zeros(m,1);
states=ones(m,1);
for ii=1:length(data)
    if data(ii,2)==last
        dur=dur+si;
    else
        loc=loc+1;
        durs(loc)=dur;
        states(loc)=last;
        dur= si;
        last= data(ii,2);
    end
end

data=cat(2,durs,states);
data=data(1:loc,:);
durations= data(:,1);
amps = data(:,2);
end


function newdata= coercetime(data,si)
    if istime(data(:,1))==true
        newdata=data;
    else
        disp("make some time using si");
        rows=length(data);
        start=0;
        times=[start:si:start+(rows-1)*si]';
        idl=data(:,2);
        newdata = cat(2,times,idl);
    end
end

function isit = istime(column)
    column=column(1:100);
    diffs=diff(column);
        if min(diffs)>=0
            isit=true;
        else
            isit=false;
        end
end
