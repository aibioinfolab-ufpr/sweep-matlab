function save_file(obj, name, ext)
% Saves the informed file with the informed name to the folder output. This
% function will be used to intermediary files required to manage the
% computers memory and the final outputs.

% Args:
%       obj: Object to save
%       name: File name to save to disk. Intermediary files must be named
%             with the extention .sweep

% Mariane Goncalves Kulik (mgkulik) - 2018-nov-13
% UFPR Bioinformatics team - http://www.bioinfo.ufpr.br/

global USER_PATH;

comp_path = char(fullfile(USER_PATH, char(strcat(name, ext))));

%try
    if strcmp(ext, '.csv')
       dlmwrite(comp_path, obj, 'delimiter','\t');
    else
        S.(char(name)) = obj;
        save(comp_path, '-struct', 'S', '-v7.3');
    end
%catch ME
%    message = strcat('Unable to save the file', {' '}, char(name), char(ext), {' '}, 'to disk.', {' '}, ...
%                   'Please check the folder output permissions. Error: ', {' '}, ME.identifier);
%    generate_log(message, 2);
%    error(ME.identifier, char(message));
%end

end