function new_fasta = process_fastas(fasta_path, varargin)
% Reads the folder with fastas or the fasta file from disk and collects
% its names and basic information/or loads it if already informed

% ARGS:
%       fasta_path: The folder name which contain all the target fastas or
%       the file name, when only one fasta file is informed.

% OPTIONAL ARGS:
%       varargin{1} - seq_type: AA (Default) or NT
%       varargin{2} - gen_dists: Values 0 or 1 to inform if the phylogenetic tree
%               must be generated.
%       varargin{3} - exec_name: The name used to dump the multifasta to disk

% RETURNS:
%       new_fasta: Multifasta file loaded or created.

% Mariane Goncalves Kulik (mgkulik) - 2018-nov-13
% UFPR Bioinformatics team - http://www.bioinfo.ufpr.br/

% Process optional arguments
if nargin > 1; seq_type = varargin{1};
else; seq_type = 'AA'; end

if nargin > 2; gen_dists = varargin{2};
else; gen_dists = 1; end

if nargin > 3; exec_name = varargin{3};
else; exec_name = 'sweep'; end

if isdeployed==0
    process_config(seq_type);
    global DEV_MODE;
    if isempty(DEV_MODE)
        setGlobalDev(1);
        DEV_MODE = getGlobalDev;
    end
else
    setGlobalDev(0);
    DEV_MODE = getGlobalDev;
end

% Takes the informed path and separates its components to access the folder
global MAIN_PATH;

[filepath,name,ext] = fileparts(fasta_path);

if isempty(ext)
    % If the extention is empty, the folder was informed and all the file
    % names are collected
    cd(fasta_path)
    fasta_ext = cellstr(['*.fa   '; '*.faa  '; '*.fasta'; '*.fas  '; '*.fna  ']);
    file_info = cellfun(@(x) dir(x), fasta_ext, 'un', 0);
    file_info = [file_info{:}];
    
    N = length(file_info);
    if N==0
        message = 'The informed path is empty or fasta files aren`t available. Please inform the correct path with protein fasta files.';
        cd(MAIN_PATH);
        generate_log(message, 2);
        error(message);
    elseif N==1
        message = 'The folder must have at least 3 organisms/strains multifasta files, otherwise phylogeny analysis will not be possible.';
        cd(MAIN_PATH);
        generate_log(message, 2);
        error(message);
    end
    
    % Dealing with empty fasta files
    ids_files = cell2mat({file_info(1:end).bytes}')>0;
    if sum(ids_files)<length(file_info)
        message = 'Some of your fasta files are empty. They will be ignored.';
        generate_log(message, 0);
    end
    file_list = {file_info(ids_files).name}';
    
    %The files size were also collected to help memory management
    %file_sizes = cell2mat({file_info(1:end).bytes}');
    %exp_size = round((sum(file_sizes)*GROWTH_RATE)/1024/1024, 2);
    
    % Returns to the source-code path
    cd(MAIN_PATH);
	   
    % Checks if the Parallel Computing toolbox is available, so we're able
    % to delete the pool when the process finishes
    parallel = 0;
    if isToolboxAvailable('Parallel Computing Toolbox')
        parallel = 1;
    end
    
    message = strcat('Total of ', {' '}, num2str(length(file_list)), {' '}, 'multifastas found.', {' '}, ...
                    'Processing the multifasta. Please wait...');
    generate_log(message, 0);
    new_fasta = generate_multifasta(file_list, fasta_path, exec_name, seq_type);
    check_fasta(new_fasta(1).Sequence, seq_type);
    
    message = 'Multifasta saved to disk. Please check the folder output.';
    generate_log(message, 0);
    
    % Deleting the parallel pool
    %if parallel == 1
    %    try
    %        poolobj = gcp('nocreate');
    %        delete(poolobj);
    %    catch
    %        message = 'Problems to stop the Parallel Pool. Please stop it manually.';
    %        generate_log(message,1);
    %        warning(message);
    %    end
    %end
    
else
    % Otherwise collects the same info from the informed file
    try
        new_fasta = fastaread(fasta_path);
    catch
        message = strcat('Unable to load the multifasta from the informed path. Please inform the correct complete file path.');
        generate_log(message, 2);
        error(message);
    end
end

check_fasta(new_fasta(1).Sequence, seq_type);

if (length(new_fasta) == 1 && gen_dists)
    message = strcat('Only 1 sequence was informed and the phylogenetic tree cannot be generated.', {' '}, ...
        'Please set the option gen_dists to 0');
    generate_log(message, 2);
    error(char(message));
end

end

%% ########################################################################
%                          AUXILIAR FUNCTIONS
% #########################################################################
function new_fasta = generate_multifasta(file_list, path, fasta_name, seq_type)
% Reads each fasta file and concatenate all its sequences appending 
% 5 asteriscs between each protein. The file name will be used as the new
% sequence header and all the original headers will be discarded.

% Function uses the toolbox "Parallel Computing Toolbox" if available.

% Args:
%       file_list: The collected list with the fasta file names;
%       path: The complete path. Just to assure the correct path is
%               accessed;
%       fasta_name: The name which will be used to save the fasta file;
%       seq_type: AA or NT. It's required to assure the right ammount of #
%               will be placed after each gene (5-AA and 15-NT)

% Returns:
%       new_fasta: The new multifasta with one organism by row

% Mariane Goncalves Kulik (mgkulik) - 2018-nov-06
% UFPR Bioinformatics team - http://www.bioinfo.ufpr.br/
    
global MAIN_PATH;
global DEV_MODE;
cd(path);

if strcmp(seq_type, 'AA')
    spacer = repmat('*', 1, 5);
else
    spacer = repmat('*', 1, 15);
end

% Parfor is 40% faster than cellfun when parallel computing is available.
parfor i=1:length(file_list)
    fas = fastaread(char(file_list(i)));
    [~,file_name] = fileparts(char(file_list(i)));
    new_fasta(i).Header = strcat(num2str(i), '_', regexprep(char(file_name), '[^a-zA-Z0-9]','_'));
    seqs = {fas(1:end).Sequence};
    new_seqs = cellfun(@(x) [x spacer], seqs,'un',0)';
    new_fasta(i).Sequence = [new_seqs{:}];
end
cd(MAIN_PATH);

if ~DEV_MODE
    save_fasta(new_fasta, fasta_name);
end

end


function comp_path = save_fasta(fasta, name)
% Saves the multifasta to disk.

global USER_PATH;
comp_path = char(fullfile(USER_PATH, char(strcat(name,'.fasta'))));
try
    message = 'Saving the multifasta to disk. Any file with the same Job Name will be replaced.';
    generate_log(message, 0);
    if exist(comp_path, 'file')
        delete(comp_path);
    end
    fastawrite(comp_path, fasta);
catch
    message = 'Unable to save the multifasta to disk. Please check the folder output permissions.';
    generate_log(message, 2);
    error(message);
end

end