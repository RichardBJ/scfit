function [ output_args ] = dwtwrite ( dwells, states, sampling, varargin )
%DWTWRITE Write idealized single channel dwell for importing into QuB
%   dwells
%   states
%   sampling
%   'fileName' (optional)
%   'Start' (optional) - Time of the start of each segment. Default = 0
%   'classCount' (optional) - Default is 2.
%   'amplitudes' (optional) - Default is 0 for class 1 and 1 for class 2
%   'sd' (optional) - Standard deviation of amplitudes

p = inputParser;
p.FunctionName = 'dwtwrite';

default_amps=[0 mean(states(states~=0))];

addRequired(p,'dwells');
addRequired(p,'state');
addRequired(p,'sampling',@isnumeric);
addParameter(p,'filename','',@ischar);
addParameter(p,'start',0,@isnumeric);
addParameter(p,'classcount',2,@isnumeric);
addParameter(p,'amplitudes',default_amps,@isnumeric);
addParameter(p,'sd',[0 0],@isnumeric);

parse(p,dwells,states,sampling,varargin{:});

fileName = p.Results.filename;
start = p.Results.start;
classCount = p.Results.classcount;
amplitudes = p.Results.amplitudes;
sd = p.Results.sd;

if isempty(fileName)
    [fileName, filePath] = uiputfile({'*.dwt','QuB idealization dwell times'});
    if fileName == 0
        msgbox('No file written');
        return
    end
    fileName = [filePath, fileName];
end

nDwells=numel(dwells);

fid = fopen(fileName, 'w');

fprintf(fid,'Segment: 1 Dwells: %d Sampling(ms): %f Start(ms): %f ClassCount: %d %f %f %f %f\n', ...
    nDwells, sampling, start, classCount, amplitudes(1), sd(1),amplitudes(2),sd(2));
% fprintf(fid,'%d\t%f\n',[+(states~=0)';dwells']);
% fprintf(fid,'%f\t%f\n',[states,dwells]');
if all( mod(states,1) == 0)
    fprintf(fid,'%d\t%e\n',[states,dwells]');
else
    fprintf(fid,'%f\t%e\n',[states,dwells]');
end
fclose(fid);
end

