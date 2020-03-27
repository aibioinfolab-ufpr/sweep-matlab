function config = load_config
% Function reads the config file and returns its values

% Returs:
%       config: a struct with all the setted configs from the file

% Mariane Goncalves Kulik (mgkulik) - 2018-dec-10
% UFPR Bioinformatics team - http://www.bioinfo.ufpr.br/

path = 'config.ini';
fileID = fopen (path,'r');
if fileID == -1
    message = 'Cannot open config.ini file. Please make sure its available in the source code folder.';
    generate_log(message, 2);
    error(message);
end
out = textscan(fileID, '%s', 'Delimiter', '\n');
out = [out{:}];

try
    iscomm = cellfun(@(x) strfind(x, '%'), out, 'un', 0);
    idx_comm = find(not(cellfun('isempty',iscomm)));
    idx_emp = find(cell2mat(cellfun(@(x) isempty(x), out, 'un', 0)));
    not_comm = out;
    not_comm([idx_comm;idx_emp]) = [];
    
    data = cellfun(@(x) strrep(x, ' ', ''), not_comm, 'un', 0);
    conf_data = cellfun(@(x) regexp(x, '=', 'split'), data, 'un', 0);
    conf_data = reshape(cellstr([conf_data{:}]'), 2, length(data))';
    conf_values = cellfun(@(x) regexp(x,',', 'split'), conf_data(:,2), 'un', 0);
    config = cell2struct(conf_values, conf_data(:,1));
catch
    message = strcat('The config.ini isn´t properly formatted. Please check our ', ...
                'documentation to set it correctly');
    generate_log(message, 2);
    error(message);
end

end