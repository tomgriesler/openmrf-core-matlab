function [twix_obj, study_info, PULSEQ] = pulseq_read_meas_siemens(study, flagRemoveOS)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% select study via selection dialog box
[~, ~, pulseq_path] = pulseq_get_user_definitions([], 0);
if nargin<1
    study = [];
end
if isempty(study)    
    study_path = [pulseq_path '/Rawdata'];
    if isempty(study_path)
        study_path = userpath();
    end
    [study_name, study_path] = uigetfile('*.*', 'Select a file', study_path );
    study_path = strrep(study_path, '\', '/');
    study = [study_path study_name];
else
    study = strrep(study, '\', '/');
    [study_path, study_name, ext] = fileparts(study);
    study_path = [study_path '/'];
    study_name = [study_name ext];
    clear ext;
end

%% optional: remove oversampling
if nargin<2
    flagRemoveOS = 0;
end
if isempty(flagRemoveOS)
    flagRemoveOS = 0;
end

%% check for file existence and emit clear exception upon missing file
checkexists(study, file_qualifier='scanner raw data');

%% Read Siemens meas file from VB/VD MRI raw data
if flagRemoveOS==0
    twix_obj = mapVBVD(study, 'ignoreSeg');
elseif flagRemoveOS==1
    twix_obj = mapVBVD(study, 'ignoreSeg', 'removeOS');
end
if numel(twix_obj)==1
    twix_obj = twix_obj{1};
end
if numel(twix_obj)>1
    twix_obj = twix_obj{2};
end

%% read header information
external_name = twix_obj.hdr.Meas.tFree;             % important! this encodes your pulseq user and time stamp
md5_hash      = twix_obj.hdr.Meas.tSequenceVariant;  % md5 hash can be used for cross check
meas_name     = twix_obj.hdr.Meas.tProtocolName;     % only the name in the MRI protocol list

% get time stamp of measurement; works only for our Skyra
mem_uid    = twix_obj.hdr.Meas.ExamMemoryUID;
[~,temp]   = sort(mem_uid=='_', 'descend');
meas_date  = mem_uid( temp(3)+1 : temp(4)-1 );
meas_clock = mem_uid( temp(4)+1 : temp(5)-1 );
clear temp;

%% import backup of pulseq workspace
[~, external_name] = fileparts(external_name);
[~,temp]      = sort(external_name=='_', 'descend');
seq_id        = external_name(1:11);
scan_id       = seq_id;
scan_id(7)    = [];
scan_id       = int64(str2num(scan_id));
pulseq_user   = external_name( 13 : temp(3)-1 );
seq_name      = external_name( temp(3)+1 : end );
if scan_id==0
    error('scan ID not found!');
end
clear temp;
backup_path = [pulseq_path() '/Pulseq_Workspace/' pulseq_user '/' seq_id(1:6) '/' seq_id];
mat_file    = dir(fullfile(backup_path, '*.mat'));
backup_path = [backup_path '/' mat_file.name];
checkexists(backup_path, file_qualifier='pulseq backup');
try
    load(backup_path);
catch
    error('no pulseq backup found!')
end

% old version: rename pulseq_workspace -> PULSEQ
if exist('pulseq_workspace', 'var')
    PULSEQ = pulseq_workspace;
    clear pulseq_workspace;
end

% compare md5 hash
if isfield(PULSEQ, 'md5_hash')
    if md5_hash ~= PULSEQ.md5_hash
        warning(' >>>>> MD5 Hashs do not agree! <<<<<')
    end
else
    warning('old pulseq version: MD5 Hash can not be compared!')
end

%% output
study_info.hdr           = twix_obj.hdr;
study_info.seq_name      = seq_name;
study_info.seq_id        = seq_id;
study_info.scan_id       = scan_id;
study_info.md5_hash      = md5_hash;
study_info.pulseq_user   = pulseq_user;
study_info.study_name    = study_name;
study_info.study_path    = study_path;
study_info.backup_file   = backup_path;
study_info.external_name = external_name;
study_info.meas_name     = meas_name;
study_info.meas_date     = meas_date;
study_info.meas_clock    = meas_clock;
study_info.f0            = twix_obj.hdr.Meas.lFrequency;
study_info.time_stamps   = twix_obj.image.timestamp(:);
if ~exist('PULSEQ','var')
   PULSEQ = 0; 
end

end

