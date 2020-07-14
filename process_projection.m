function projMat = process_projection(proj, mask, exec_name, fasta_type)
% Function checks if a projection was chosen or if its size was informed
% and generates a new random projection matrix, than saves it to file with
% the informed name.

% Args:
%       proj: The projection matrix path or its size
%       mask: The mask which will be used in the HDV transformation allows
%           the row proj matrix definition;
%       exec_name: The name which will be used to save the projection to
%           disk.

% Returns:
%       projMat: The loaded random projection matrix

% Mariane Goncalves Kulik (mgkulik) - 2018-nov-16
% UFPR Bioinformatics team - http://www.bioinfo.ufpr.br/

% Generates the random projection matrix if only its size be informed and 
% dump it to disk

if strcmp(fasta_type,'AA')
    hdv_size = 20^mask(1)*20^mask(3);
    name = strcat(exec_name, 'AA_', num2str(proj));
else
    hdv_size = 4^mask(1)*4^mask(3);
    name = strcat(exec_name, 'NT_', num2str(proj));
end

if isnumeric(proj)
    message = strcat('Generating the projection matrix size', {' '}, num2str(hdv_size),...
                     'x', num2str(proj), '. Please wait...');
    generate_log(message, 0);
    projMat = orthbase(hdv_size, proj);
    
    message = strcat('Dumping the projection', {' '}, name, {' '}, 'to disk. Please wait...');
    generate_log(message, 0);
    save_file(projMat, name, '.sweepproj');
else
    % Loads the informed projection matrix
    try
        message = 'Loading the projection matrix. Please wait...';
        generate_log(message, 0);
        projMat = load('-mat', proj);
        field_name = char(fieldnames(projMat));
        projMat = projMat.(field_name);
    catch
        message = 'Unable to load the projection matrix. Please check if the informed path is correct. If you are working with gene/transcript sequences make sure you downloaded the default projection matrix from: https://sourceforge.net/projects/spacedwordsprojection/';
        generate_log(message, 2);
        error(message);
    end
    
    projSize = size(projMat);
    if (strcmp(fasta_type,'AA') && projSize(1)~= hdv_size)
        message = 'The informed projection matrix must have 160,000 (20^4) values for AA sequences.';
        generate_log(message, 2);
        error(message);
    elseif (strcmp(fasta_type,'NT') && projSize(1)~= hdv_size)
        message = 'The informed projection matrix must have 65,536 (4^8) values for NT sequences.';
        generate_log(message, 2);
        error(message);
    end
end