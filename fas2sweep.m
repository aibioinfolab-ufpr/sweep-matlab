function [finalMat, headers, projMat] = fas2sweep(xfas, varargin)
% Function catalogs the presence of the mask using the slide window by
% sequence.

% ARGS:
%   xfas = Multifasta file with the list of genes of one organism or the
%          list of genes of each organism separated by ** (2 asteriscs) in
%          a single line

% OPTIONAL ARGS: 
%   varargin{1} - fasta_type: Info of fasta type AA or NT
%	varargin{2} - Defines the HDV composition:
%        0: Binary matrix, to indicate absence and presence;
%        1: Prime number of the absence and presence positions;
%        2: Occurrence counts for the present values;
%   varargin{3} - Defines the function output:
%        DEFAULT: The sparse original matrix. 
%                 #### The maximum of sequences sugested to a 16GB computer
%                 is 3.000 sequences. More then this value may cause slowness. ####
%        matProj: The output will be the Orthonormal projection of the
%                 original matrix. #### Requires an Orthonormal matrix of 
%                 random basis to execute the projection.
%   varargin{4} - k_hdv: Flag indicating if the user expects to access
%                 the HDV file later. If not requested, those objects will be
%                 discarded.
%   varargin{5} - mask - The mask with 3 positions containing the first 
%                 part, the mismatch and mismatch and the last part. Allowed 
%                 values: [1 1 3] [2 0 2] [2 1 2] [2 0 3] [3 0 1] [3 0 2] [3 1 1]
%                 ([4 7 4] - NT only)
%   varargin{6} - exec_name: Job name will be used to save the files to disk
%               if not in DEV_MODE.
%   varargin{7} - Defines the order of the entries. #### Requires a vector
%                of the size of the fasta with the new ordering values.

% RETURNS:
%       headers: Sequencial headers
%       finalMat: HDV sparse matrix
%       projMat: SWeeP Projected Matrix


% Roberto Tadeu Raittz (rraittz) - 2018-nov-06
% Mariane Goncalves Kulik (mgkulik) - 2018-nov-06
% UFPR Bioinformatics team - http://www.bioinfo.ufpr.br/

% mgkulik (Update): 2019-feb-15

% Starting variables
message = 'Starting the HDV generation and transformation to SWEEP. Please wait...';
generate_log(message, 0);

withPos = 0;                % Default matrix values are binary
proj = 0;                   % The default is not execute the ortonormal projection. It avoid errors when the projection matrix isn't informed
mask = [2 1 2];             % The default mask is 2 AA, 1 GAP, 2 AA
Ps = primon(1);             % Sets the first prime number as 1
k_hdv = 0;                  % By default the HDV are discarded.
fasta_type = 'AA';          % Default operation is AA, but NT can be used
defSize = 20;               % Default value to AA (20). Will be changed to NT (4)
exec_name = 'sweep';        % Name used only if not informed

% Extracts the sequences from the fasta file
fas_cell = struct2cell(xfas);
headers = fas_cell(1,:)';
seqs = fas_cell(2,:);

% Checks the optional variable values
if nargin>1
    fasta_type = varargin{1};
    if strcmp(fasta_type,'NT')
        defSize = 4;
    end
end

if nargin>2
    if varargin{2} == 1                                                                                                                                                                                                           
        withPos = 1;
        % Gets maximum sequence size to pre-calculate the primes
        max_seqSize = max(cellfun(@length,seqs));
        Ps = primon(max_seqSize);
    elseif varargin{2} == 2
        withPos = 2;
    end
end

if nargin>3
    proj = 1; 
    matProj = varargin{3};
end

if nargin>4; k_hdv = varargin{4};end
if nargin>5
    mask = varargin{5};
    if length(mask) ~= 3
        message = strcat('Mask must be an array with 3 integer values:', {' '}, ...
                        '1: Fist portion size, 2: Gap size, 3: Last portion size.');
        generate_log(message, 2);
        error(message);
    end 
end

if nargin>6; exec_name = varargin{6}; end
if nargin>7
    iord = varargin{7};
    seqs = seqs(iord);
    headers = headers(iord);
end

% Checking if all sequences are bigger than de mask size
seq_size = cellfun(@length, seqs);
headers_small = char(headers(seq_size < sum(mask)));
if ~isempty(headers_small)
    headers_small = sprintf([headers_small repmat('; ', length(cellstr(headers_small)),1)]');
    message = strcat('There are sequences smaller than the mask size.', {' '}, ...
                        'Please remove them from the fasta file:', {' '}, ...
                        char(headers_small));
    generate_log(message, 2);
    error(char(message));
end

% Executing the sliding mask count by chunk. As the struggle here is memory
% usage, we manage the execution by slices.
chunks = calculate_chunks(length(xfas));
idx = generate_chunk(chunks, length(xfas));

if k_hdv
    finalMat = cell(chunks, 1);
end


if proj==1
    projMat = zeros(length(xfas), size(matProj,2));
end

if isdeployed==0
    process_config(fasta_type);
    global DEV_MODE;
    if isempty(DEV_MODE)
        setGlobalDev(1);
        DEV_MODE = getGlobalDev;
    end
else
    setGlobalDev(0);
    DEV_MODE = getGlobalDev;
end

%% Executes the transformations itself
for i = 1:chunks
    parcM = seqs(idx(i,1):idx(i,2));
    W160k = cell2mat(cellfun(@(x) matinline(mask2vec(x, Ps, withPos, mask, defSize)),parcM, 'un', 0)');
    
    % Saves each HDV to a cell to later save it to disk.
    finalMat{i} = sparse(W160k);
    
    % Executes the ortonormal projection with the provided matrix
    rows_size = defSize^mask(1)*defSize^mask(3);
    if (proj == 1)
        try
            projMat(idx(i,1):idx(i,2), :) = full(W160k*matProj);
        catch
            message = char(strcat('Projection matrix dimensions must agree with the informed mask.', {' '}, ...
                        'The informed mask requires the projection of size ', {' '}, ...
                        num2str(rows_size) ,{' '}, 'by M. Please check.'));
            generate_log(message, 2);
            error(message);
        end
    end
end


%% Saves the files to disk if not in DEV_MODE
if ~DEV_MODE
    clear matProj;
    % Saving the SWEEP to disk
    try
        name = strcat(exec_name, '_SWEEP');
        save_file(projMat, name, '.csv');
        message = strcat('Low dimension vector (SWEEP) file', {' '}, name, '.csv dumped to disk.', {' '}, ...
                             'Available at the output folder.');
        generate_log(message, 0);
    catch ME
        message = strcat('Unable to dump the low dimension vector (SWEEP) file to disk.', {' '}, ...
                             'Please check the `output` folder permissions. Error: ', {' '}, ME.identifier);
        generate_log(message, 2);
        error(message);
    end
end

% Execute the next step only when required.
if k_hdv||DEV_MODE
    finalMat = cell2mat(finalMat');
    finalMat = reshape(finalMat, length(seqs), rows_size);
end

% Saves the HDV single file to disk, if required and not in DEV_MODE
if k_hdv
    try
        name = strcat(exec_name, '_HDV');
        save_file(finalMat, name, '.mat');
        message = strcat('File', {' '}, name, '.mat dumped to disk.', {' '}, ...
                     'You`ll find it at the `output` folder.');
        generate_log(message, 0);
    catch
        message = strcat('Unable to dump the high dimension vector (HDV) file to disk.', {' '}, ...
                     'Please check the `output` folder permissions.');
        generate_log(message, 2);
        error(message);
    end
end


end

%% ########################################################################
%                          AUXILIAR FUNCTIONS
% #########################################################################
function chunks = calculate_chunks(fas_len)
% Calculates the number of parts the input fasta must be divided to assure
% it fits memory. Default value is 4, but depending on the ammount lines it
% can be increased.

global MAX_SEQS;

chunks = ceil(fas_len/MAX_SEQS);

end

function mret = matinline(M)
% Reshapes the matrix into a 1xN array

% Args:
%       M: The matrix to reshape

% Returns:
%       mret: The 1xN array

% Roberto Tadeu Raittz (rraittz) - 2018-nov-09
% UFPR Bioinformatics team - http://www.bioinfo.ufpr.br/

mret = M(find(M==M))';

end