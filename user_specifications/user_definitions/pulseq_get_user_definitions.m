function [pulseq_user, pulseq_lab, pulseq_path, pulseq_scanner, system] = pulseq_get_user_definitions(pulseq_scanner, flag_output)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

    if nargin==0
        pulseq_scanner = [];
        flag_output = 0;
    end
    if nargin==1
        flag_output = 0;
    end

    %% get path of user .csv file
    path_def = which('pulseq_get_user_definitions');
    path_def = strrep(path_def(1:end-29), '\', '/');
    path_def = [path_def 'pulseq_user_definitions.csv'];

    %% read user definitions from file
    if exist(path_def, 'file') > 0
        temp_defs   = readcell(path_def);
        pulseq_user = temp_defs{1,2};
        pulseq_lab  = temp_defs{2,2};
        pulseq_path = temp_defs{4,2};
        if isempty(pulseq_scanner)
            pulseq_scanner = temp_defs{3,2};
        end
        clear temp_defs;
    
    %% create new .csv file
    else

        % set pulseq username
        pulseq_user = inputdlg('Enter your pulseq username:', 'Pulseq User', [1 50], {'anonymous'});
        pulseq_user = pulseq_user{1};

        % set MRI lab
        labList = {'University of Wuerzburg, Department of Physics, EP5, Wuerzburg, Germany', ...
                   'University of Michigan, Department of Biomedical Engineering, Ann Arbor, MI, USA', ...
                   'University Medical Center Freiburg, Department of Radiology, Medical Physics, Freiburg, Germany', ...
                   'University Hospital Jena, Institute for Diagnostic and Interventional Radiology, Medical Physics Group, Jena, Germany', ...
                   'University of Erlangen-Nuernberg, Department of Neuroradiology at Uniklinikum Erlangen, Erlangen, Germany', ...
                   'University of Bordeaux, IHU Liryc, Bordeaux, France', ...
                   'University Hospital Wuerzburg, Department of Diagnostic and Interventional Radiology, Wuerzburg, Germany', ...
                   'University Hospital Wuerzburg, Department of Neuroradiology, Wuerzburg, Germany', ...
                   'University Hospital Wuerzburg, Department of Radiotherapy and Radiation Oncology, Wuerzburg, Germany', ...
                   'University Hospital Wuerzburg, Comprehensive Heart Failure Center, Wuerzburg, Germany', ...
                   'Other (enter manually)'};        
        [idx, ~] = listdlg('ListString', labList, ...
                           'ListSize', [600, 200], ...
                           'SelectionMode', 'Single', ...
                           'PromptString', 'Select you Lab or choose "Other" to enter manually', ...
                           'Initialvalue', 1, ...
                           'Name', 'Make choice');        
        if idx == numel(labList)
            pulseq_lab = inputdlg('Enter the name of your MRI lab:', 'Custom Lab Name', [1 60]);
            pulseq_lab = pulseq_lab{1};
        else
            pulseq_lab = labList{idx};
        end
        clear labList idx;

        % set pulseq backup path
        pulseq_path = uigetdir('/', 'select Path of Pulseq Worskpace data');
        pulseq_path = strrep(pulseq_path, '\', '/');

        % set default MRI scanner
        scanner_path = pulseq_get_path('pulseq_find_scanner_csv');
        scanner_path = dir(fullfile(scanner_path, '*.csv'));
        for j=1:numel(scanner_path)
            scannerList{j,1} = strrep(scanner_path(j).name(1:end-4), '_', ' ');
        end
        [idx, ~] = listdlg('ListString', scannerList, ...
                   'ListSize', [300, 250], ...
                   'SelectionMode', 'Single', ...
                   'PromptString', 'Select you default MRI scanner', ...
                   'Initialvalue', 1, ...
                   'Name', 'Make choice'); 
        pulseq_scanner = scanner_path(idx).name(1:end-4);
        clear scanner_path scannerList idx;

        % save definitions in .csv file
        writelines({['user;',    pulseq_user]; ...
                    ['lab;',     pulseq_lab]; ...
                    ['scanner;', pulseq_scanner]; ...
                    ['path;',    pulseq_path]}, ...
                    path_def);
    end

    %% correct path for unix or windows
    if isunix
        if pulseq_path(2) == ':'
            pulseq_path = strrep(['/' pulseq_path(1) pulseq_path(3:end)], '\','/');
        end
    end
    if ispc
        if pulseq_path(1) == '/'
            pulseq_path = strrep([pulseq_path(2) ':' pulseq_path(3:end)], '\','/');
        end
    end

    %% init pulseq system struct
    scanner_path = pulseq_get_path('pulseq_find_scanner_csv');
    scanner_path = fullfile(scanner_path, [pulseq_scanner '.csv']);
    temp_tab     = readtable(scanner_path, setvartype(detectImportOptions(pulseq_scanner), 'char'));    
    system       = mr.opts( 'B0',                   str2num(temp_tab(1,2).Var2{1}), ...
                            'gamma',                str2num(temp_tab(2,2).Var2{1}), ...
                            'maxGrad',              str2num(temp_tab(3,2).Var2{1}), ...
                            'gradUnit',             temp_tab(4,2).Var2{1}, ...
                            'maxSlew',              str2num(temp_tab(5,2).Var2{1}), ...
                            'slewUnit',             temp_tab(6,2).Var2{1}, ...
                            'maxB1',                mr.convert(str2num(temp_tab(7,2).Var2{1}), 'uT'), ...
                            'rfDeadTime',           str2num(temp_tab(8,2).Var2{1}), ...
                            'rfRingdownTime',       str2num(temp_tab(9,2).Var2{1}), ...
                            'adcDeadTime',          str2num(temp_tab(10,2).Var2{1}), ...
                            'adcRasterTime',        str2num(temp_tab(11,2).Var2{1}), ...
                            'rfRasterTime',         str2num(temp_tab(12,2).Var2{1}), ...
                            'gradRasterTime',       str2num(temp_tab(13,2).Var2{1}), ...
                            'blockDurationRaster',  str2num(temp_tab(14,2).Var2{1}) );    
    system.ascfile = temp_tab(15,2).Var2{1};
    if strcmp(system.ascfile, '')
        system.ascfile = [];
    end

    %% output results
    if flag_output==1
        disp(['     ... pulseq lab:       ' pulseq_lab]);
        disp(['     ... pulseq username:  ' pulseq_user]);
        disp(['     ... pulseq path:      ' pulseq_path]);
        disp(['     ... system specs:     ' pulseq_scanner]);
    end

end