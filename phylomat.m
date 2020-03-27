function [tree_values, dists] = phylomat(ldv, headers, varargin)
% Calculates the parwise distances using the informed distance type and
% generates the Neighbor Join clustering file to open in external tools.

% ARGS:
%       ldv: The low dimension vector for each sequence
%       headers: Each sequence header

% OPTIONAL ARGS: 
%       varargin{1} - t_dist: Selected distance method (Euclidean - 
%               Default, Spearman or Cosine)
%       varargin{2} - tree_type: Selected clustering method (Ward -
%               Default, NeighborJoin)
%       varargin{3} - exec_name: The execution name will be used to save
%               the tree file to disk if not in DEV_MODE
%       varargin{4} - seq_type: AA (Default) or NT

% RESULTS:
%       tree_values: The tree information formatted to open in external
%               visualization tools.

% Roberto Tadeu Raittz (rraittz) - 2018-nov-18
% Mariane Goncalves Kulik (mgkulik) - 2018-nov-18
% UFPR Bioinformatics team - http://www.bioinfo.ufpr.br/

% UPDATES:
%   mgkulik - 2019-apr-09 - Ward parameter addition

if nargin > 2; t_dist = varargin{1};
else; t_dist = 'euclidean'; end

if nargin > 3; tree_type = varargin{2};
else; tree_type = 1; end

if nargin > 4; exec_name = varargin{3};
else; exec_name = 'sweep'; end

if nargin > 5; seq_type = varargin{4};
else; seq_type = 'AA'; end

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

global USER_PATH;

message = 'Starting the tree generation. Please wait...';
generate_log(message, 0);

% Generates the tree object
dists = pdist(ldv,t_dist);

if tree_type==1
    tree_values = seqneighjoin(dists,'equivar',headers);
else
    tree_values = linkage(dists,'ward');
end
message = strcat('Dumping ', {' '}, exec_name, '.tree to disk. It will be available at the folder output.');
generate_log(message, 0);

% Dumps the tree to disk
if ~DEV_MODE
    try
        exec_name = char(exec_name);
        save_file(squareform(dists), strcat(exec_name, '_DISTS'), '.csv');
        comp_path = fullfile(USER_PATH, strcat(exec_name,'.tree'));
        if tree_type==1
            phytreewrite(comp_path, tree_values);
        else
            gen_tree = phytree(tree_values,headers);
            phytreewrite(comp_path, gen_tree);
        end
    catch
        message = strcat('Unable to save the file ', exec_name, '.tree to disk. ', ...
                        'Please check the folder output permissions.');
        generate_log(message, 2);
        error(message);
    end
end