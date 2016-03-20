function [ dwells, states, ndwells, amps, startTime, nSegs, segs ] = dwtread( varargin )
%DWTREAD Read idealization data from QuB
%   INPUT: fileName
%         'joinSegs' (optional default is true)
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
%
%  ------------------------------------------------------------------------
%  According to QuB DWT documentation (www.qub.buffalo.edu/qubdoc/files/dwt.html), 
%  DWT
% 
% idealized data (event list)
% ASCII text
% segmented (can be discontiguous according to seg. start time)
% segments are separated by a blank line
% each segment begins with a one-line header, which can be incomplete:
%   Segment: %d -- numbered starting with 1
%   Dwells: %d -- number of events in this segment
%   Sampling(ms): %f -- sampling rate of the original recording
%   Start(ms): %f -- time-offset of this segment from beginning of recording
%   ClassCount: %d -- number of distinct amplitude levels
%   {%f %f}* -- {amp, std} of each class
% subsequent lines are dwells of "<class> \t <duration>"
% classes are numbered consecutively starting with 0. Two dwells have the
%   same class if they are indistinguishable (same amplitude)
% durations are given in milliseconds
% Example:
% 
% Segment: 1 Dwells: 17 Sampling(ms): 1 Start(ms): 0 ClassCount: 2 -0.00563935 0.102412 1.02029 0.141354
% 0	134
% 1	3
% 0	3
% 1	6
% 0	6
% 1	14
%
% -------------------------------------------------------------------------
% 
% In practice, however, I have not encountered dwt files that have segments
% separated by blank lines. Rather the segment one-line header appears on
% the line immediately following a dwell.

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
    % Since the segment header line can be incomplete, first read in the
    % line and then parse, checking for each part of the header
    % c = textscan(fid, 'Segment: %u Dwells: %u Sampling(ms): %f Start(ms): %f ClassCount: %d',1);
    
    % Parse segment header line
    % ---------------------------------------------------------------------
    line = fgetl(fid);
    
    % Force the segment header to at least state it's a segment header
    [c, nrmatches] = sscanf(line, 'Segment:%u Dwells:%u');
    if (nrmatches < 2 || isempty(c))
        error('Cannot find "Segment: " in %s\nLine:\n%s', fileName, line);
    end
    segs(ii) = c(1);
    ndwells(ii)=c(2);
    
    % look for other parts of header, which may be missing
    % start time
    c = regexpi(line, 'Start\(ms\):\s+([0-9.]+)', 'tokens');
    if ~isempty(c)
        startTime(ii)=str2double(c{1});
    end
    
    % Class Count
    c = regexpi(line, 'ClassCount:\s+([0-9]+)', 'tokens');
    if ~isempty(c)
%         disp(c{1});
        nr_classes = str2double(c{1});
        c = regexpi(line, 'ClassCount:\s+[0-9]+(.*)$', 'tokens');
        % if there is a match, c is a cell array and each element of the
        % cell array is a cell array of length equal to the number of
        % tokens
        % we expect a single match and we've only got a single token
        if ~isempty(c)
            temp = sscanf(c{1}{1}, '%f %f', [nr_classes, 2]);
            % sscanf returns an array, not a cell array, like textscan
            % if successful, return a nrclassess-by-2 matrix with column
            % one giving the amplitudes and column 2 giving the std of each
            % class
            amps{ii} = temp;
        end
    end
    
    % Read dwells from this segment
    % ---------------------------------------------------------------------
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

