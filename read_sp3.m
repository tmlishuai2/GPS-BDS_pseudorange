function sp3 = read_sp3(files, folder, satsys, start_time, end_time)
% READ_SP3 reads SP3 files.
%
% SYNTAX:
%	[sp3, hdr] = read_sp3(sp3_file);
%
% INPUT:
%	files - SP3 files
%   folder -
%   satsys =
%
% OUTPUT:
%   sp3 - matrix of the satellite position, velocity and clock data
%	hdr - Struct of the sp3 file header
%   sp3(:,1:2) = gpsweek, sow]
%   sp3(:,3)   = sat number
%   sp3(:,4:7) = [x,y,z,dtr]
%   sp3(:,8:11)= [x_dot,y_dot,z_dot,dtr_dot]
%
% EXAMPLE:
%	sp3 = read_sp3 ('d:\0422\rinex3\gbm18935.sp3');
%
% See also READ_RINEX_N, READ_SP3, FIXY2K, CAL2GPST.
%
% Copyright 2002-2015  zhengdong.bai@gmail.com

% validate number of input arguments
narginchk(1,5);

if nargin<2, folder = []; end
if nargin<3 || isempty(satsys), satsys = 'GCRESJI'; end
if nargin<4, start_time  = []; end
if nargin<5, end_time    = []; end

if ischar(files), files = cellstr(files); end

sp3 = [];
n = length(files);
if(n == 0), return; end

for i=1:length(files)
    sp3_file = fullfile(folder, files{i});
    [sp3i,hdri] = read_sp3_file(sp3_file);  
    if ~isempty(hdri) && ~isempty(sp3i)
        sp3 = [sp3; sp3i]; %#ok<AGROW>
    end
end

ii = all(isnan(sp3(:,1:11)),2);
sp3(ii,:) = [];

prn = char(sp3(:,3)/100);
gps = []; bds = []; glo = []; gal = []; sbs = []; qzs = []; irn = [];

if(find(satsys=='G')>0), ii = (prn(:,1) == 'G'); gps = sp3(ii,:); end
if(find(satsys=='C')>0), ii = (prn(:,1) == 'C'); bds = sp3(ii,:); end
if(find(satsys=='R')>0), ii = (prn(:,1) == 'R'); glo = sp3(ii,:); end
if(find(satsys=='E')>0), ii = (prn(:,1) == 'E'); gal = sp3(ii,:); end
if(find(satsys=='S')>0), ii = (prn(:,1) == 'S'); sbs = sp3(ii,:); end
if(find(satsys=='J')>0), ii = (prn(:,1) == 'J'); qzs = sp3(ii,:); end
if(find(satsys=='I')>0), ii = (prn(:,1) == 'I'); irn = sp3(ii,:); end

sp3 = [gps; bds; glo; gal; sbs; qzs; irn];

% remove the observation data before the first epoch or after the last epoch
if ~isempty(start_time)
    ts = gpst2sec(sp3(:,1:2));
    ii = ts >= gpst2sec(cal2gpst(start_time));
    sp3 = sp3(ii,:);
end

if ~isempty(end_time)
    ts = gpst2sec(sp3(:,1:2));
    ii = ts <= gpst2sec(cal2gpst(end_time));
    sp3 = sp3(ii,:);
end

if ~isempty(sp3)
    [~,ii] = unique(sp3(:,1:3),'rows');
    sp3 = sp3(ii,:);
end

end

%-------------------------------------------------------------------------------
function [sp3, hdr] = read_sp3_file(sp3_file)

% validate number of input arguments
narginchk(1,1);

% % see if this file has already been formatted for MATLAB
mat_file = sprintf('%s.mat',lower(strtrim(sp3_file)));
if exist(mat_file, 'file')
    mat = load(mat_file, 'hdr', 'sp3');
    hdr = mat.hdr;
    sp3 = mat.sp3;
    return
end

% read the whole file to a temporary cell array
[fid, message] = fopen(sp3_file,'rt');
if fid == -1
    error ('Open %s failed: %s.\n', sp3_file, message);
else
    buf = textscan(fid,'%s','delimiter','\n','whitespace','');
    buf = buf{1};
end
fclose(fid);

% check the 1st line for version symbol
line  = [buf{1}, blanks(256)];
if ~strcmp(line(1:2), '#a') && ~strcmp(line(1:2), '#c')
    error ('%s is NOT a SP3 file', sp3_file);
end

% read observation file data from the buf
n = length(buf);
sp3 = nan(n,11); % preallocate memory
k = 0;

for i = 2:n
    line  = [buf{i}, blanks(256)];
    if strcmpi(line(1:3), 'EOF'), break; end
    symbol = line(1);

    switch symbol
        case '*'
            [date,m] = sscanf(line(2:31),'%f%f%f%f%f%f');
            if(m==6)
                gpst = cal2gpst(date');
            end
        case 'P'
            sat = prn2sat(line(2:4));
            [pos,m] = sscanf(line(5:60),'%f%f%f%f');
            if(~isempty(sat) && m==4)
                k = k + 1;
                sp3(k,1:3) = [gpst,sat];
                sp3(k,4:4+m-1) = pos;
            end
        case'V'
            sat = prn2sat(line(2:4));
            [vel,m] = sscanf (line(5:60),'%f%f%f%f');
            if(sat == sp3(k,3) && m==4)
                sp3(k,8:8+m-1) = vel;
            end
        otherwise
    end
end

% remove the extra preallocated memory
sp3(k+1:end,:) = [];

% sort and remove duplicates
if ~isempty(sp3)
    [~,ii] = sortrows(sp3,1:3);
    sp3 = sp3(ii,:);
    hdr.gpst = [sp3(1,1:2), sp3(end,1:2)];
else
    hdr.gpst = [NaN, NaN];
end

% save to a mat file for later use
save(mat_file, 'hdr', 'sp3');

end
