function [ dwells, states, ndwells, amps, startTime, nSegs, segs ] = dwtread( varargin )
%DWTREAD Read idealization data from QuB
%   INPUT: fileName
%         'joinSegs' (optional default is true
%
%   OUTPUT: dwells - duration of sojourns in each state
%           states - identity of each state of dwells
%           ndwells - number of dwells in each segment
%           amps - a 1-by-k cell array where k is the number of segments
%               each cell is a m-by-2 matrix where m is the number of
%               states; the first column is the mean amplitude and the
%               second column is the sd
%           startTime - start time in milliseconds of each segment
%           nSegs - number of segments in idealization
%           segs - the segment identity that corresponds to the data file
%               or source of the idealization

p = inputParser;
addOptional(p,'fileName','',@ischar);
addParameter(p,'joinsegs',true,@islogical);
parse(p,varargin{:});

fileName = p.Results.fileName;
joinSegs = p.Results.joinsegs;

if isempty(fileName)
    [fileName, filePath] = uigetfile({'*.dwt','QuB idealization dwell times'});
    if fileName == 0
        dwells = 0;
        states = 0;
        ndwells = 0;
        amps = 0;
        startTime = 0;
        return
    end
    fileName = [filePath, fileName];
    assignin('base','fileName',fileName);
end

% Note: Opening the file in "text mode" (see help for fopen) resulted in
% problems reading a file with multiple segments.

fid = fopen(fileName);
ii=1;
while ~feof(fid)
    c = textscan(fid, 'Segment: %u Dwells: %u Sampling(ms): %f Start(ms): %f ClassCount: %d',1);
    segs(ii) = c{1};
    ndwells(ii)=c{2};
    startTime(ii)=c{4};
    temp = textscan (fid, '%f %f', c{1,5});
    amps{ii}(:,1) = temp{1};
    amps{ii}(:,2) = temp{2};
%     temp = textscan(fid,'%d %f');
    temp = textscan(fid,'%f %f');
    states{ii} = temp{1};
    dwells{ii} = temp{2};
    ii = ii+1;
end

if size(dwells,2) ~= 1
    dwells=dwells';
    states=states';
end

nSegs = ii-1;
if joinSegs || nSegs==1
    dwells=cell2mat(dwells);
    states=cell2mat(states);
end

fclose (fid);

end

