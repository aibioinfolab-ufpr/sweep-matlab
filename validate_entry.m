function entry_data = validate_entry(fasta_path, varargin)
% Loading function, used only by the command line version to check the
% entry parameters.

% Args:
%       fasta_path: The complete path to the fastas to execute the
%               pre-processing tasks or the complete fasta path + file name
%               or a single multifasta complete path.
%       type: The type of the informed entry must be AA, NT or Genome. 
%               Later in the process we'll assure the corret type was 
%               informed, but we must validate the default mask.
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

% Returns:
%       entry_data: A struct with 8 fields, informed or replaced by its
%       default.

% Mariane Goncalves Kulik (mgkulik) - 2018-nov-26
% Mariane Goncalves Kulik (mgkulik) - 2019-feb-15
% UFPR Bioinformatics team - http://www.bioinfo.ufpr.br/

if nargin>1
    param = varargin{:};
end

% Fasta parameter
[filepath_fas,name_fas,ext_fas] = fileparts(fasta_path);
if ~isempty(ext_fas)
    if (strcmp(ext_fas, '.fa')&&strcmp(ext_fas, '.faa')&&strcmp(ext_fas, '.fasta')&&strcmp(ext_fas, '.fas')&&strcmp(ext_fas, '.fna'))
        message = 'Invalid fasta file. Please check.';
        generate_log(message, 2);
        error(message);
    elseif isempty(filepath_fas)
        fasta_path = fullfile(pwd, char(strcat(name_fas,ext_fas)));
    else
        if ~exist(filepath_fas, 'dir')
            check = strfind(fasta_path(1), '/');
            if (check > 0)
                fasta_path = fasta_path(2:end);
            end
        end
    end
else
    if ~exist(fasta_path, 'dir')
        check = strfind(fasta_path(1), '/');
        if (check > 0)
            fasta_path = fasta_path(2:end);
        end
        if ~exist(fasta_path, 'dir')
            message = strcat('Invalid fasta folder (', fasta_path, '). Please check.');
            generate_log(message, 2);
            error(message);
        end
    end
end
entry_data.fasta = fasta_path;


% Check sequence type
type = vararginnamed(param, 'SeqType', 'AA');
%if ~(strcmpi(type, 'AA') || strcmpi(type, 'NT') || strcmpi(type, 'Genome'))
if ~(strcmpi(type, 'AA') || strcmpi(type, 'NT'))
    message = 'Please inform the correct entry type: AA or NT.';
    generate_log(message, 2);
    error(message);
end
entry_data.type = type;


% Check if the user wants to generate the phylogenetic tree
gen_dists = vararginnamed(param, 'GenerateTree', 1);
gen_dists = cast(gen_dists, 'logical');
%whos gen_dists;
%disp(gen_dists==0);
%disp(gen_dists==1);
%disp(strcmpi(gen_dists, '0'));
%disp(strcmpi(gen_dists, '1'));
if ~((gen_dists==0)||(gen_dists==1))
    message = 'Please inform 0 or 1 to generate the Phylogenetic Tree.';
    generate_log(message, 2);
    error(message);
end
entry_data.gen_dists = gen_dists;


% Check which distance method should be used
dist_type = vararginnamed(param, 'DistanceMethod', 'euclidean');
if (strcmpi(dist_type, 'euclidean')...
        &&strcmpi(dist_type,'spearman')&&strcmpi(dist_type, 'cosine'))
    message = 'Invalid distance method. Please inform: Euclidean OR Spearman OR Cosine.';
    generate_log(message, 2);
    error(message);
end
entry_data.dist_type = dist_type;


% Check which tree type should be used
tree_type = vararginnamed(param, 'ClusteringType', 1);
tree_type = cast(tree_type, 'int8');
if (~(tree_type==1||tree_type==2))
    message = 'Invalid clustering method. Please inform: 1 (Ward) or 2 (SeqNeighborJoin).';
    generate_log(message, 2);
    error(message);
end
entry_data.tree_type = tree_type;

k_hdv = vararginnamed(param, 'SaveHDV', 0);
k_hdv = cast(k_hdv, 'logical');
% Check if the user wants to save the HDV file
if ~((k_hdv==0)||(k_hdv==1))
    message = 'Please inform 0 or 1 to save the binary matrix to disk (default 0).';
    generate_log(message, 2);
    error(message);
end
entry_data.k_hdv = k_hdv;


% Check informed method
method = vararginnamed(param, 'SWeePMethod', 'binary');
if ~(strcmpi(method,'binary')||strcmpi(method, 'prime')||strcmpi(method, 'count'))
    message = 'Invalid method. Please inform: Binary OR Prime OR Count.';
    generate_log(message, 2);
    error(message);
end
entry_data.method = method;


% Check if projection is integer or, if file path, if the extension is
% correct and if the complete path or just the file name was informed.
proj = vararginnamed(param, 'ProjMatrix', 'DEFAULT');
if ~(isnumeric(proj)||strcmpi(proj, 'DEFAULT'))
    [filepath,name,ext] = fileparts(proj);
    if ~isempty(ext)
        if ~strcmp(ext, '.sweepproj')
            message = 'Please inform a valid .sweepproj file.';
            generate_log(message, 2);
            error(message);
        else
            if ~exist(filepath, 'dir')
                check = strfind(filepath(1), '/');
                if (check > 0)
                    filepath = filepath(2:end);
                end
                if ~exist(filepath, 'dir')
                    proj = fullfile(pwd, char(strcat(name, ext)));
                else
                    proj = fullfile(filepath, char(strcat(name, ext)));
                end
            end
        end
    else
        message = 'Please inform a valid .sweepproj file.';
        generate_log(message, 2);
        error(message);
    end
end
entry_data.proj = proj;


% Check mask
if strcmpi(type, 'AA')
    num_mask = str2num((vararginnamed(param, 'SWeepMask', '11011')));
elseif strcmpi(type, 'NT')
    num_mask = str2num((vararginnamed(param, 'SWeepMask', '111100000001111')));
%else
%    num_mask = str2num((vararginnamed(param, 'SWeepMask', '111110000011111')));
end
if ~isnumeric(num_mask)
    message = 'Invalid mask format. Please follow the example: 11011, 1 for keep and 0 for GAP.';
    generate_log(message, 2);
    error(message);
else
    mask = num2str(num_mask);
end
entry_data.mask = mask;

% Generates the structure with the data values
%field_names = {'fasta', 'type', 'gen_dists', 'dist_type', 'tree_type', 'k_hdv', 'method', 'proj', 'mask'};
%entry_data = cell2struct([cellstr(fasta_path) cellstr(type) num2cell(gen_dists) cellstr(dist_type) num2cell(tree_type) num2cell(k_hdv) cellstr(method) proj cellstr(mask)], field_names, 2);