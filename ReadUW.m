function [hdr, seis, masthdr] = ReadUW(filename, varargin)
% [hdr,seis] = ReadUW(filename);
% [hdr,seis] = ReadUW(filename, byteorder, verbosity, datapath, pickpath);
% [hdr,seis,masthdr] = ReadUW(...);
%
% This function reads UW2 format data or pick files (e.g. 05011611163W)
% into Matlab. The old UW interface libraries are NOT required.
%
% INPUTS
%
% filename      String, e.g. '04093023532W'. Accepts a pickfile name, a
%               datafile name, or either one with no trailing letter (e.g.
%               04093023532). Program behavior varies with string spec:
%
%                   Spec                Behavior
%                   [event ID][a-z]     Read data file and pick file
%                   [event ID][D,W]     Read data file only
%                   [event ID]          Read data file and pick file
%
% OPTIONAL INPUTS
%
% byteorder     char (default 'n'). See 'help fopen' for allowed formats.
%
% verbosity     int (default '0'). 
% 
% datapath      String for name of environmental variable giving path to UW
%               data. Default: USER_DATA_PATH (if set), else PATH.
%
% pickpath      String for name of environmental variable with path to pick
%               files. Defaults to PATH.  
%
% OUTPUTS
%
% hdr           UW2 channel header and pickfile structure. Values
%               correspond to channel numbers. 
%
% seis          Seismograms, a cell structure that references by channel
%               number: seis{1} = data from channel 1, etc.  
%
% masthdr       UW2 master header. A structure of master header values. Not
%               especially useful, except for debugging. 
%
% NOTE: Unprocessed UW-format waveform data won't work unless the data file
% contains a valid channel structure.
% 
% =======================================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 2.0, 2016-03-05

%% Parse varargin
byteorder = 'n';
datapath = 'USER_DATA_PATH';
pickpath = 'PATH';
verbose = 0;
if numel(varargin) > 0
    if ~isempty(varargin{1})
        byteorder = varargin{1};
    end
    if numel(varargin) > 1
        if ~isempty(varargin{2})
            verbose = varargin{2};
        end
        if numel(varargin) > 2
            if ~isempty(varargin{3})
                datapath = varargin{3};
            end
            if numel(varargin) > 3
                if ~isempty(varargin{4})
                    pickpath = varargin{4};
                end
            end
        end
    end
end

% Added 2015-08-29 to deal with non-UW paths
if isempty(getenv(datapath)); datapath = 'PATH'; end
if isempty(getenv(pickpath)); pickpath = 'PATH'; end

%% Locate files
ext = (filename(numel(filename)));
if strcmp(ext,'W') || strcmp(ext,'D')
    if verbose
        disp('UW data file specified. Reading data file only.');
    end
    df1 = findfile(filename,datapath);
elseif regexp(ext,'[a-z]')
    if verbose
        disp('UW pick file specified. Reading pick and data files.');
    end
    if strcmp(ext,'d')
        dfname = strcat(filename(1:numel(filename)-1),'D');
        df1 = findfile(dfname,datapath);
        pfname = filename;
        pf1 = findfile(pfname,pickpath);
    else
        dfname = strcat(filename(1:numel(filename)-1),'W');
        df1 = findfile(dfname,datapath);
        pfname = filename;
        pf1 = findfile(pfname,pickpath);
    end
elseif regexp(ext,'[0-9]')
    if verbose
        disp('Generic UW file specified. Searching for data and pick files.');
    end
    dfext = {'D' 'W'};
    df1 = findfile(filename,datapath,dfext);
    pfext = {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' ...
             'k' 'l' 'm' 'n' 'o' 'p' 'q' 'r' 's' 't' ...
             'u' 'v' 'w' 'x' 'y' 'z'};
    pf1 = findfile(filename,pickpath,pfext);
end

% Open data file
fid = fopen(df1,'r',byteorder);

% Initial number of structs
nstructs = 1;

%% Process master header
masthdr.nchan       = fread(fid,1,'int16');
masthdr.lrate       = fread(fid,1,'int32');
masthdr.lmin        = fread(fid,1,'int32');
masthdr.lsec        = fread(fid,1,'int32');
masthdr.length      = fread(fid,1,'int32');
masthdr.tapenum     = fread(fid,1,'int16');
masthdr.eventnum	= fread(fid,1,'int16');
masthdr.flags       = transpose(fread(fid,10,'int16'));
masthdr.extra       = strcat(transpose(fread(fid,10,'*char')));
masthdr.comment     = strcat(deblank(transpose(fread(fid,80,'*char'))));

% Set masthdr time using lmin and lsec
[masthdr.yr, masthdr.mo, masthdr.dy, masthdr.hr, masthdr.mn, masthdr.sc] = ...
    datevec(datenum(masthdr.lmin/1440 + masthdr.lsec/8.64e10) + ...
            datenum(1600,1,1));

% Seek to end of file; get number of structures
fseek(fid, -4, 'eof');
nstructs = fread(fid, nstructs, 'int32');
structs_os = (-12*nstructs)-4;
tc_os = 0;

% Set version of UW seismic data file
if strcmp(masthdr.extra(3),'2')
    uwformat = 2;
else
    uwformat = 1;
end

% Read in UW2 data structures to determine number of channels
if uwformat == 2
    fseek(fid,structs_os,'eof');
    for i1 = 1:1:nstructs
        structtag           = strcat(transpose(fread(fid,4,'*char')));
        nstructs            = fread(fid,1,'int32');
        byteoffset          = fread(fid,1,'int32');
        if strcmp(structtag,'CH2')
            hdr.nchan = nstructs;
        elseif strcmp(structtag,'TC2')
            fpos = ftell(fid);
            fseek(fid, byteoffset,-1);
            chno = zeros(nstructs,1);
            corr = zeros(nstructs,1);
            for n = 1:1:nstructs
                chno(n)     = 1+fread(fid,1,'int32');
                corr(n)     = fread(fid,1,'int32');
            end
            tc_os = -8*nstructs;
            fseek(fid, fpos, -1);
        end
    end
else
    hdr.nchan = masthdr.nchan;
end

if verbose
    disp(['Processing ' num2str(hdr.nchan) ' channels.']);
end

%% Write time corrections to header
hdr.timecorr = zeros(1, hdr.nchan);
if exist('chno','var')
    for n = 1:1:numel(chno)
        hdr.timecorr(chno(n)) = corr(n)/8.64e10;
    end
end

%% Read UW2 channel headers
if uwformat == 2
    fseek(fid,(-56*hdr.nchan)+structs_os+tc_os,'eof');
    fmt = cell(hdr.nchan,1);
    for i1=1:hdr.nchan
        hdr.chlen(i1)       = fread(fid,1,'int32');
        hdr.offset(i1)  	= fread(fid,1,'int32');
        hdr.start_lmin(i1) 	= fread(fid,1,'int32');
        hdr.start_lsec(i1)	= fread(fid,1,'int32');
        
        % Set time for each trace
        [hdr.yr(i1), hdr.mo(i1), hdr.dy(i1), ...
            hdr.hr(i1), hdr.mn(i1), hdr.sc(i1)] = ...
            datevec(datenum((hdr.start_lmin(i1)/1440) + ...
            (hdr.start_lsec(i1)/8.64e10)) + hdr.timecorr(i1) + ...
            datenum(1600,1,1));
        
        hdr.Fs(i1) 		= fread(fid,1,'int32')/1000;
        hdr.expan1(i1)	= fread(fid,1,'int32');
        hdr.lta(i1)		= fread(fid,1,'int16');
        hdr.trig(i1)	= fread(fid,1,'int16');
        hdr.bias(i1)	= fread(fid,1,'int16');
        hdr.fill(i1)	= fread(fid,1,'int16');
        hdr.name{i1}    = strcat(transpose(fread(fid,8,'*char')));

        tmp             = strcat(transpose(fread(fid,4,'*char')));
        for j1=1:numel(tmp)
            if strcmpi(tmp(j1),'f')
                fmt{i1} = 'float32';
            elseif strcmpi(tmp(j1),'l')
                fmt{i1} = 'int32';
            elseif strcmpi(tmp(j1),'s')
                fmt{i1} = 'int16';
            end
        end
        hdr.compflg{i1} = strcat(transpose(fread(fid,4,'*char')));
        hdr.chid{i1}    = strcat(transpose(fread(fid,4,'*char')));
        hdr.expan2{i1}  = strcat(transpose(fread(fid,4,'*char')));
    end
end

%% Read UW channel data
if uwformat == 2
    seis = cell(1,hdr.nchan);
    for i1=1:hdr.nchan
        fseek(fid,hdr.offset(i1),'bof');
        seis{i1} = fread(fid,hdr.chlen(i1),fmt{i1});
    end
end

%% Done data file
fclose(fid);
if verbose
    disp('Done reading data file.');
end

%% Pickfile processing
if exist('pf1','var')
    
    % Does a pickfile exist?
    if exist(pf1,'file') == 2
        if verbose
            disp('Reading pickfile');
        end
        pfid=fopen(pf1,'r',byteorder);
        
        % _______________________________________________________________
        % Acard line
        acard = nextline(pfid,'A');
        
        ycorr = 0;
        if numel(acard) == 75 || numel(acard) == 12
            y2k = 0;
            ycorr = 1900;
        else
            y2k = 2;
        end
        hdr.type = acard(2);
        
        ot = zeros(1,6);
        ot(1) = str2double(acard(3:4+y2k))+ycorr;
        ot(2) = str2double(acard(5+y2k:6+y2k));
        ot(3) = str2double(acard(7+y2k:8+y2k));
        ot(4) = str2double(acard(9+y2k:10+y2k));
        ot(5) = str2double(acard(11+y2k:12+y2k));
        
        if numel(acard) > (14+y2k)
            ot(6)           = str2double(acard(13+y2k:18+y2k));
            hdr.lat         = str2double(acard(19+y2k:21+y2k)) + ...
                                str2double(acard(23+y2k:24+y2k))/60 + ...
                                str2double(acard(25+y2k:26+y2k))/6000;
            latcode         = acard(22+y2k);
            if strcmp(latcode,'S')
                hdr.lat     = -hdr.lat;
            end
            hdr.lon         = str2double(acard(27+y2k:30+y2k)) + ...
                                str2double(acard(32+y2k:33+y2k))/60 + ...
                                str2double(acard(34+y2k:35+y2k))/6000;
            loncode         = acard(31+y2k);
            if strcmp(loncode,'W')
                hdr.lon = -hdr.lon;
            end
            hdr.z           = str2double(acard(36+y2k:41+y2k));
            hdr.fix         = acard(42+y2k);
            hdr.mag         = str2double(acard(43+y2k:46+y2k));
            hdr.numsta      = str2double(acard(47+y2k:49+y2k));
            hdr.numpha      = str2double(acard(51+y2k:53+y2k));
            hdr.gap         = str2double(acard(54+y2k:57+y2k));
            hdr.dmin        = str2double(acard(58+y2k:60+y2k));
            hdr.rms         = str2double(acard(61+y2k:65+y2k));
            hdr.err         = str2double(acard(66+y2k:70+y2k));
            hdr.q           = acard(71+y2k:72+y2k);
            hdr.velmodel    = acard(74+y2k:75+y2k);
        elseif numel(acard) > 12+y2k
            hdr.region      = acard(14+y2k);
        end
        
        % Set origin time at each station
        hdr.ot = zeros(1,hdr.nchan);
        for n = 1:1:hdr.nchan
            v = [hdr.yr(n) hdr.mo(n) hdr.dy(n) hdr.hr(n) hdr.mn(n) ...
                 hdr.sc(n)];
            hdr.ot(n) = 86400*(datenum(ot)-datenum(v));
        end
        
        % _______________________________________________________________
        % Error line
        fseek(pfid,0,'bof');
        eline = nextline(pfid,'E');
        
        if eline ~= -1
            hdr.MeanRMS     = str2double(eline(12:17));
            hdr.SDabout0    = str2double(eline(18:23));
            hdr.SDaboutMean = str2double(eline(24:29));
            hdr.sswres      = str2double(eline(30:37));
            hdr.ndfr        = str2double(eline(38:41));
            hdr.fixxyzt     = eline(42:45);
            hdr.SDx         = str2double(eline(46:50));
            hdr.SDy         = str2double(eline(51:55));
            hdr.SDz         = str2double(eline(56:60));
            hdr.SDt         = str2double(eline(61:65));
            hdr.SDmag       = str2double(eline(66:70));
            hdr.MeanUncert  = str2double(eline(76:79));
        else
            fseek(pfid,0,'bof');
        end
        
        % _______________________________________________________________
        % Alternate magnitude line
        fseek(pfid,0,'bof');
        sline = nextline(pfid,'S');
        
        if sline ~= -1
            warning('Alternate magnitude found, M_d overwritten.');
            hdr.mag         = str2double(sline(1:5));
            hdr.magtype     = sline(6:8);
        else
            fseek(pfid,0,'bof');
        end
        
        % _______________________________________________________________
        % Focal mechanism line
        fseek(pfid,0,'bof');
        m1 = 0;
        mline = nextline(pfid,'M');
        while mline ~= -1
            m1 = m1+1;
            hdr.mech{m1,1}  = str2double(mline(5:7));
            hdr.mech{m1,2}  = str2double(mline(9:10));
            hdr.mech{m1,3}  = str2double(mline(14:16));
            hdr.mech{m1,4}  = str2double(mline(18:19));
            hdr.mech{m1,5}  = str2double(mline(23:25));
            hdr.mech{m1,6}  = str2double(mline(27:28));
            hdr.mech{m1,7}  = str2double(mline(32:34));
            hdr.mech{m1,8}  = str2double(mline(36:37));
            hdr.mech{m1,9}  = str2double(mline(41:43));
            hdr.mech{m1,10} = str2double(mline(45:46));
            hdr.mech{m1,11} = str2double(mline(50:52));
            hdr.mech{m1,12} = str2double(mline(54:55));
            mline = nextline(pfid,'M');
        end
        
        % _______________________________________________________________
        % Pick lines
        fseek(pfid,0,'bof');
        m1 = 0;
        pline = nextline(pfid,'.');
        while pline ~= -1
            m1=m1+1;
            pline = pline(2:numel(pline));
            sta = pline(1:3);
            cmp = pline(5:7);
            
            for j1=1:numel(hdr.name)
                if strcmpi(hdr.name{j1},sta) && ...
                        strcmpi(hdr.compflg{j1},cmp)
                    n = strfind(pline,' ');
                    pline = pline(n(1)+1:end);
                    while ~isempty(pline)
                        [P,n] = textscan(pline, '%s', 1, 'delimiter',')');
                        t = char(P{1});
                        switch t(2)
                            case 'P'
                                C = textscan(t(4:end),'%s %s %f %d %f %f');
                                switch char(C{1})
                                    case 'P'
                                        hdr.P(j1)       = C{3};
                                        hdr.Ppol{j1}    = char(C{2});
                                        hdr.Pqual(j1)   = C{4};
                                        hdr.Punc(j1)    = C{5};
                                        if ~isempty(C{6})
                                            hdr.Perr(j1)= C{6};
                                        end
                                    otherwise
                                        hdr.S(j1)       = C{3};
                                        hdr.Spol{j1}    = char(C{2});
                                        hdr.Squal(j1)   = C{4};
                                        hdr.Sunc(j1)    = C{5};
                                        if ~isempty(C{6})
                                            hdr.Serr(j1)= C{6};
                                        end
                                end
                            case 'D'
                                C = textscan(t(2:end),'%s %f');
                                hdr.D(j1) = C{2};
                            otherwise
                                warning('Unknown pick packet!');
                        end
                        pline = pline(n+1:end);
                    end
                end
            end
            pline = nextline(pfid,'.');
        end
        
        % Done file read
        fclose(pfid);
        if verbose
            disp('Done reading pickfile.');
        end
        
        % _______________________________________________________________
        % Postprocess
        
        % Ensure that all channels have P, S, and duration set
        if isfield(hdr,'P') 
            if numel(hdr.P) < numel(hdr.name)
                hdr.P(numel(hdr.P)+1:numel(hdr.name)) = 0;
            end
        else
            hdr.P = zeros(size(hdr.name));
        end
        
        if isfield(hdr,'S') 
            if numel(hdr.S) < numel(hdr.name)
                hdr.S(numel(hdr.S)+1:numel(hdr.name)) = 0;
            end
        else
            hdr.S = zeros(size(hdr.name));
        end
        
        if isfield(hdr,'D') 
            if numel(hdr.D) < numel(hdr.name)
                hdr.D(numel(hdr.D)+1:numel(hdr.name)) = 0;
            end
        elseif ~isfield(hdr,'D')
            hdr.D = zeros(size(hdr.name));
        end
        
        % If a pick is zero, insert placeholders
        for j1=1:numel(hdr.P)
            if hdr.P(j1) == 0
                hdr.Ppol{j1}    = [];
                hdr.Pqual(j1)   = 4;
                hdr.Punc(j1)    = 1;
                hdr.Perr(j1)    = 10;
            end
            if hdr.S(j1) == 0
                hdr.Spol{j1}    = [];
                hdr.Squal(j1)   = 4;
                hdr.Sunc(j1)    = 1;
                hdr.Serr(j1)    = 10;
            end
        end
    else
        warning('Pickfile not read (invalid path or filename).');
    end
end

% Clean up hdr structure
hdr = rmfield(hdr,'offset');
hdr = rmfield(hdr,'chlen');
hdr = orderfields(hdr);

% =======================================================================
function tmpstring = nextline(pfid,ch)
tmpstring = fgetl(pfid);
while tmpstring(1) ~= ch
    tmpstring = fgetl(pfid);
    if tmpstring == -1
        break
    end
end

% =======================================================================
function fname = findfile(varargin)
filename = varargin{1};
upath = getenv(varargin{2});
ext = {''};
if nargin == 3
    ext = varargin{3};
end

% Dummy check: long format filename, correctly specified
if exist(filename,'file') == 2
    fname = filename;
    return
end

% Try to find data file in current directory
for j1 = 1:numel(ext)
    filename1 = strcat(filename,ext{j1});
    fname = fullfile(pwd,filename1);
    test1 = exist(fname,'file');
    if test1 == 2
        break
    end
end

% Try to find data file in other path directories if this fails
while test1 ~= 2
    if ~isempty(upath)
        P = textscan(upath, '%s', 'delimiter', pathsep);
        i = 1;
        while i < numel(P)
            for j1 = 1:numel(ext)
                fname = fullfile(P{i},strcat(filename,ext{j1}));
                test1 = exist(fname,'file');
                if test1 == 2
                    break
                end
            end
            i = i+1;
        end
    else
        if strcmp(ext{1},'a')
            warning(['No file matching ' filename '[a-z] found in path!']);
        elseif strcmp(ext{1},'G')
            warning(['No file matching ' filename '[G,W] found in path!']);
        end
        break
    end
end
