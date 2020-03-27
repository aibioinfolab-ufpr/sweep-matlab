function [sldv, headers, hdv] = execute_sweep(multifasta, mask, method, proj_mat, exec_name, keep_hdv, fasta_type)
% Loads the fasta files, if not in memory, and executes the the function to
% generate the HDV and project it to the LDV.

% Args:
%       multifasta: The multifasta in memory or its prefix to load from
%                   disk
%       mask: The mask to use in the slide window
%       method: The HDV composition (Binary, Prime or Count)
%       proj_mat: The orthonormal projection matrix to execute the LDV
%       exec_name: Job name will be used to save the HDV to disk
%       keep_hdv: Flag to indicate if the .mat HDV file must be created
%       fasta_type: Info of fasta type AA or NT

% Returns:
%       SWeeP LDV: The lower dimension vector for each sequence;
%       headers: Each sequence header

% Mariane Goncalves Kulik (mgkulik) - 2018-nov-16
% UFPR Bioinformatics team - http://www.bioinfo.ufpr.br/

if isstruct(multifasta)

    % Executes with the transformation function to generate the LDV with
    % its ordered headers.
    [hdv, headers, sldv] = fas2sweep(multifasta, fasta_type, method, proj_mat, keep_hdv, mask, exec_name);
    
else
	% Otherwise, loads the multifasta to memory.
    try
        xfas = fastaread(multifasta);
        [hdv, headers, sldv] = fas2sweep(xfas, fasta_type, method, proj_mat, keep_hdv, mask, exec_name);
    catch ME
        switch ME.identifier
            case 'bioinfo:fastaread:FastaNotValid'
                message = strcat('Unable to load the multifasta from the folder output.', {' '}, ...
                                'Please check if the informed name is correct.');
                generate_log(message, 2);
                error(message);
            otherwise
                message = strcat('Please contact the developer. Error:', {' '}, ME.identifier);
                generate_log(message, 2);
                error(message);
        end
    end
end