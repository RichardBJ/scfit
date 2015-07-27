function [ durations, amplitudes, properties, calibration, varargout] = scanread( filename )
%SCANREAD Read transition properties from a SCAN.SCN file
%   Detailed explanation goes here

if nargin<1 || isempty(filename)
    [filename, filePath] = uigetfile({'*.scn','SCAN idealization file'});
    if filename == 0
        durations = 0;
        amplitudes = 0;
        properties = 0;
        calibration = 0;
        return
    end
    filename = [filePath, filename];
    assignin('base','fileName',filename);
end

fid = fopen(filename,'r');

fseek(fid,4,'bof');
dataOffset = fread(fid,1,'int');
dataOffset = dataOffset-1;

nTransitions = fread(fid,1,'int');

fseek(fid,208,'bof');
nfits = fread(fid,1,'int32');
ntmax = fread(fid,1,'int32');
nfmax = fread(fid,1,'int32');

fseek(fid,dataOffset,'bof');
durations = fread(fid,nTransitions,'float');
amplitudes = fread(fid,nTransitions,'int16');
properties = fread(fid,nTransitions,'int8');
minfreq = fread(fid,3,'int32');
maxfreq = fread(fid,3,'int32');
num2read = maxfreq-minfreq+1;
allpoint = fread(fid,num2read(1),'int32');
openpoint = fread(fid,num2read(2),'int32');
shutpoint = fread(fid,num2read(3),'int32');

transtimes = fread(fid,nfits,'real*8');
% transtimes(ii) = t0sav + dfinter*(infirst(ii)-1);
idxfits = fread(fid,nfits,'int32');
baseline = fread(fid,nfits,'int16');
% pos = ftell(fid);


fseek(fid,425,'bof');
calibration = fread(fid,1,'float');
% amplitudes = calibration*amplitudes;

fseek(fid,624,'bof');
imin = fread(fid,1,'int32');
imax = fread(fid,1,'int32');

fseek(fid,469,'bof');
dfinter = fread(fid,1,'real*8');
tlast = fread(fid,1,'real*8');

fseek(fid,501,'bof');
t0sav = fread(fid,1,'float');

fseek(fid,517,'bof');
infit = fread(fid,1,'int32');
infirst = fread(fid,1,'int32');

fseek(fid,393,'bof');
fc = fread(fid,1,'float');

offset_error = 2048;
infit = infit - offset_error;
infirst = infirst - offset_error;
transtimes = transtimes - offset_error*dfinter;
tlast = tlast - offset_error*dfinter;
t0sav = t0sav - offset_error*dfinter;

varargout = {minfreq, maxfreq, allpoint, openpoint, shutpoint, transtimes, ...
    idxfits, baseline, imin, imax, dfinter, tlast, t0sav, infit, infirst, ...
    fc, ntmax, nfmax};

fclose(fid);

end

