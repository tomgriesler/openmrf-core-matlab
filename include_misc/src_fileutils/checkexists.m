function checkexists(filepath, options)
    arguments
        filepath char
        options.file_qualifier char = ''
    end
    if ~isempty(options.file_qualifier)
        qualifier = [' ' options.file_qualifier];
    else
        qualifier = options.file_qualifier;
    end
    % Check that a file exists at filepath
    if isfolder(filepath)
        error('OpenMRF:IsADirectoryError', ...
              'Expected a%s file but found a directory at ''%s''', ...
              qualifier, filepath);
    elseif ~isfile(filepath)
        error('OpenMRF:FileNotFoundError', ...
              'Required%s file at ''%s'' does not exist!', ...
              qualifier, filepath);
    end 
end