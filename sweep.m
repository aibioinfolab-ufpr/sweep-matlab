function [sldv, headers, exec_data, new_fasta, proj_mat, tree, dists, hdv] = sweep(exec_data, varargin)
%function sweep(exec_data, varargin)
% SWeeP main function. Executes all the steps required to generate the
% SWeeP package products.

% ARGS (8 variables from console):
%       fasta_path: The complete path to the fastas to execute the
%               pre-processing tasks or the complete fastatree path + file name
%               or a single multifasta complete path.
%       type: The type of the informed fasta AA or NT. Later in the process
%               we'll assure the corret type was informed, but it must be
%               informed to validate the default mask.
%       gen_dists: Values 0 or 1 to inform if the phylogenetic tree
%               must be generated. Default is 1.
%       dist_type: Selected method to calculate the tree distances. Default
%               is Euclidean.
%       k_hdv: Values 0 or 1 to inform if the LDV matrix should be dumped
%               to disk. Default is 0.
%       method: The HDV transformation method. Default is Binary;
%       proj_info: The projection size or a projection complete path, with 
%               file name. This projection has extention .sweepproj and 
%               was save from a previous SWeeP execution.
%               Its reuse will assure identical trees over the same
%               fasta data. Default is 600.
%       mask: The mask to the slide window, required by the HDV
%               transformation. Suggested default 11011 for AA and 
%               111100000001111 to NT.

% OUTPUTS (dumped to disk):

%   CSV format
%       SWeeP LDV (Lower Dimension Vectors - SWEEP) - The projected matrix 
%                   generated by processing the SROM and the HDV. This is 
%                   the source to calculate the tree distances.

%   FASTA format
%       SWeeP MultiFasta file(s) - Miltifasta transformed file.


%   MATLAB format
%       SWeeP HDV (High Definition Vectors - SHDV ) - Sparse matrix with 
%                   the original HDV projections.
%       SWeeP Random Orthonormal Matrix (SROM) - Orthonormal random
%                   projections used to reduce dimensionality with minimum 
%                   information loss. Saved only in case of a new one is 
%                   requested;
%       SWeeP Tree - Phylogenetic tree generated through processing the
%                   SLDV using Neighbor Join clustering method.


% Roberto Tadeu Raittz (rraittz) - 2018-nov-26
% Mariane Goncalves Kulik (mgkulik) - 2018-nov-26
% UFPR Bioinformatics team - http://www.bioinfo.ufpr.br/

if isdeployed==0
    global DEV_MODE;
    if isempty(DEV_MODE)
        setGlobalDev(1);
        DEV_MODE = getGlobalDev;
    end
else
    setGlobalDev(0);
    DEV_MODE = getGlobalDev;
end

init_time = clock;

message = 'SWeeP started. You will be updated each step status.';
generate_log(message, 0);
generate_log(' ', 3);

% Check if the values were informed using the UI. If not, they are checked
% and added to a struct, so the application can manage both entries the
% same way.
%if ~isstruct(exec_data)
exec_data = validate_entry(exec_data, varargin);
%end

config_data = process_config(exec_data.type);

% Execute additional data transformations, applied to both entry ways.
exec_data = ui_transform(exec_data);   

% Reads the fasta(s) to generate the multifasta formated file
tic;
new_fasta = process_fastas(exec_data.fasta, exec_data.type, exec_data.gen_dists, exec_data.exec_name);
sec_time = toc;
tot_duration = num2hms(sec_time);
message = strcat('Total time to generate and save/load the multifasta:', {' '}, char(tot_duration));
generate_log(message, 0);
generate_log(' ', 3);


% Manages the random projection matrix (Generate or load it)
tic;
proj_mat = process_projection(exec_data.proj, exec_data.mask, exec_data.exec_name, exec_data.type);
sec_time = toc;
tot_duration = num2hms(sec_time);
message = strcat('Time to generate and save/load the Random projection matrix:', {' '}, char(tot_duration));
generate_log(message, 0);
generate_log(' ', 3);


% Generates the HDV and the projected matrix LDV
tic;
[sldv, headers, hdv] = execute_sweep(new_fasta, exec_data.mask, exec_data.method, proj_mat, exec_data.exec_name, exec_data.k_hdv, exec_data.type);
sec_time = toc;
tot_duration = num2hms(sec_time);
message = strcat('Time to generate the HDV and LDV and dump to disk:', {' '}, char(tot_duration));
generate_log(message, 0);
generate_log(' ', 3);


% Generates the tree object and saves to disk
tree_msg = ' ';
tree_msg2 = ' ';
if exec_data.gen_dists
    tic;
    [tree, dists] = phylomat(sldv, headers, exec_data.dist_type, exec_data.exec_name, exec_data.tree_type, exec_data.type);
    sec_time = toc;
    tot_duration = num2hms(sec_time);
    tree_msg =  '+Tree';
    tree_msg2 = strcat('You can load the .tree file generated in the folder output ', {' '}, ...
                 'with the viewers (dendroscope.org) or (itol.embl.de).');
    message = strcat('Time to generate the tree file:', {' '}, char(tot_duration));
    generate_log(message, 0);
    generate_log(' ', 3);
else
    tree = "";
    dists = "";
end

sec_time = etime(clock, init_time);
tot_duration = num2hms(sec_time);

message = strcat('SWeeP execution:', {' '}, exec_data.exec_name, {' '}, 'finished.', {' '}, ...
                 'Elapsed time: ', {' '} , tot_duration, ' (SWeeP', tree_msg,')', {' '}, tree_msg2);
generate_log(message, 0);

end

function struct_data = ui_transform(struct_data)
% Additional tranformation required to guarantee the informed values are in
% the expected format to be used by the functions

% Removes the characters which cause trouble to a file name. It'll be used
% as name for all the output files.
if isfield(struct_data, 'exec_name')
    exec_name = replace(replace(struct_data.exec_name, ' ', '_'), '.', '_');
else
    %If not informed, use actual date/time as entry
    now_vect = datevec(datestr(now,'yyyy-mm-dd HH:MM'));
    exec_name = strcat('sweep_', num2str(now_vect(1)),'_', num2str(now_vect(2)),...
                        '_', num2str(now_vect(3)),'_',num2str(now_vect(4)),...
                        '_',num2str(now_vect(5)), '_files');
end
struct_data.exec_name = exec_name;


% Adjust the mask to a vector format expected by fas2sweep
if strcmp(struct_data.type, 'AA')
    if strcmp(struct_data.mask,'11011'); mask = [2 1 2];
    elseif strcmp(struct_data.mask, '10111'); mask = [1 1 3];
    elseif strcmp(struct_data.mask, '11101'); mask = [3 1 1];
    elseif strcmp(struct_data.mask, '101'); mask = [1 1 1];
    elseif strcmp(struct_data.mask, '1101'); mask = [2 1 1];
    elseif strcmp(struct_data.mask, '1011'); mask = [1 1 2];
    elseif strcmp(struct_data.mask, '1001'); mask = [1 2 1];
    else
        message = 'Invalid mask. Please check the readme.txt and inform an allowed AA mask.';
        generate_log(message, 2);
        error(message);
    end
else
    if strcmp(struct_data.mask,'111110000011111'); mask = [5 5 5];
    elseif strcmp(struct_data.mask,'111100000001111'); mask = [4 7 4];
    elseif strcmp(struct_data.mask, '111100001111'); mask = [4 4 4];
    elseif strcmp(struct_data.mask, '111000000111'); mask = [3 6 3];
    elseif strcmp(struct_data.mask, '111000111'); mask = [3 3 3];
    elseif strcmp(struct_data.mask, '110011'); mask = [2 2 2];
    elseif strcmp(struct_data.mask, '101'); mask = [1 1 1];
    else
        message = 'Invalid mask. Please check the readme.txt and inform an allowed NT mask.';
        generate_log(message, 2);
        error(message);
    end
end
struct_data.mask = mask;

% Transforms the method text in numeric value, expected by fas2sweep
if strcmpi(struct_data.method, 'binary'); method = 0;
elseif strcmpi(struct_data.method, 'prime position')||strcmpi(struct_data.method, 'prime'); method = 1;
elseif strcmpi(struct_data.method, 'count'); method = 2; 
end
struct_data.method = method;

if strcmp(struct_data.proj, 'DEFAULT')
    global PROJ_PATH;
    struct_data.proj = PROJ_PATH;
end
    
end

function time_string = num2hms(time)
% Functions transforms the tic/toc ouptut to print in the log file
time_string = datestr(datenum(0,0,0,0,0,time),'HH:MM:SS');
end