function save_study_results(study_name, res, default_dir)
    if ~isempty(study_name)
        [~, savename] = fileparts(study_name);
        savename = sprintf('%s.mat', savename);
    else
        savename = '.mat';
    end
    if nargin < 3
        [filename, pathname] = uiputfile('*.mat', 'Save results as', savename);
    else 
        [filename, pathname] = uiputfile('*.mat', 'Save results as', fullfile(default_dir, savename));
    end
    if nargin==1
        fprintf('no variables to save')
        return 
    end

    % If user didn't cancel, save to the chosen location
    if isequal(filename, 0) || isequal(pathname, 0)
        disp('Save canceled by user.');
    else
        fullpath = fullfile(pathname, filename);
        save(fullpath, 'res', '-v7.3');
        fprintf('Saved variables to: %s\n', fullpath);
    end
end