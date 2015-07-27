function [ ] = copy_and_rename_results( destination )
%COPY_AND_RENAME_RESULTS Rename fit results with model and data file names
%   Input
%   destination (optional) - destination folder

narginchk(0,1);

if nargin<1
    destination = uigetdir('','Select Destination Folder');
    if destination==0
        return
    end
end
    
[files, path] = uigetfile('*.csv','Select Results Files to Move','MultiSelect','on');
if files==0
    return
end

for ii=1:length(files)
    fid = fopen([path, files{ii}]);
    for jj=1:3
        fgetl(fid);
    end
    str = fgetl(fid);
    str = str(13:end);
    [~,mdl,~] = fileparts(str);
    str = fgetl(fid);
    str = str(22:end);
    [~,datafile,~] = fileparts(str);
    newname = [mdl, ' - ', datafile, '.csv'];
    fclose(fid);
    copyfile([path, files{ii}],[destination, newname]);
end