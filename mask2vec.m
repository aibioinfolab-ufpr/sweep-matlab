function [mret, xy] = mask2vec(xseq, Ps, withPos, mask, defSize)
% Executes all the steps to check for each sequence the absence or presence
% of the amino acids group, using the sliding mask.

% Example:
% SEQ:     MKGASLLGATLL      mask = [2 1 2]
% Slide: 1 MK-AS
%        2  KG-SL
%        3   GA-LL
%        4    AS-LG
%        5     SL-GA
%        6      LL-AT
%        7       LG-TL
%        8        GA-LL

% Each of the occurrences above will be registered, but GA-LL occurs at
% lines 3 and 8. All the non-registered occurrences, e.g. AA-AA will be
% cataloged as missing (zero). The variable withPos will define how the
% founded values will be registered in the matrix.

% Args:
%       xseq: Sequence to analyse;
%       Ps: List of primes. It will be used only with withPos = 1;
%       withPos: Type of expected register for presence in the matrix
%              (0: binary; 1: prime of the position; 2: occurrences count);
%       mask: The mask to be used in the slide window

% Returns:
%       mret: Sparse vector with information with max of 160 thousand
%             positions of absence and presence of amino acid groups
%             determined by the mask over the slide window;
%       xy: The matrix position indexes for each of the presence detected

% Roberto Tadeu Raittz (rraittz) - 2018-nov-08
% Mariane Goncalves Kulik (mgkulik) - 2018-nov-08
% UFPR Bioinformatics team - http://www.bioinfo.ufpr.br/

t = mask(1) + mask(2) + mask(3);

% Slices the sequence through slide window with the mask total size and
% defines its position in the final matrix if present.
slice = seq2list(xseq, t);
xy = [aa2idx(slice(:,1:mask(1)), defSize) aa2idx(slice(:,t-(mask(3)-1):t), defSize)];
inot = ~prod(1-double(xy<=0),2);
xy(inot,:) = [];

lines = defSize^mask(1);
inds = ij2inds(xy,lines);

% Assures no index is bigger than the matrix size and initializes the 
% output matrix
cols = defSize^mask(3);
inds(inds>(lines*cols)) = [];
M = zeros(lines,cols,'double');

% Insert the prime of the position for the present values. This method will
% allow reverse the matrix to the original sequence
if withPos==1
    [vals, inds] = unifysameinds(inds, Ps(1:length(inds)), @(x) prod(x.^(1/length(x))));
    M(inds) = vals;

% Insert the counts for all the present values
elseif withPos==2
    accum_counts = accumarray(inds, 1, [cols*lines 1]);
    M(1:cols*lines) = accum_counts;

% Insert the value one to the present values. Even if more than one mask 
% value was detected in the sequence, the number 1 will be inserted to
% generate a binary matrix
else
    M(inds) = 1;
end

% Returns the sparse matrix
mret = sparse(M);

%% Auxiliary functions
function mret = seq2list(sdna, size)
% Generates a list of amino acids through the slide window with the total
% size of the informed mask

% Args:
%       sdna: The sequence to use
%       size: The size of the window to slide

% Returns:
%       mret: A cell list with the slided amino acid sequence

% Mariane Gonçalves Kulik (mgkulik) - 2018-nov-09
% UFPR Bioinformatics team - http://www.bioinfo.ufpr.br/

n = length(sdna);

xcol = uint32(1:(n-size+1))';
inds = repmat(xcol,1,size) + repmat(uint32(1:size),(n-size+1),1)-1;
mret = sdna(inds);

if size==1
    mret = mret';
end


function mret = aa2idx(xseq, defSize)
% Uses the 20 amino acid numeric representation to calculate the position
% of the combined AA in the matrix.

% Example:
% AAgroup - Positions
% AA   -    1
% AR   -    21
% AN   -    41
% AC   -    61
% AD   -    81
% AQ   -    101

% Args:
%      xseq: Amino acid group to transform. Are considered groups each
%            portion of the mask, gap excluded.
%      defSize: The values considered to AA (20) and NT (4). Required to
%            define the final matrix size.

% Returns:
%      mret: The AA numeric index adjusted to fit the matrix

% Roberto Tadeu Raittz (rraittz) - 2018-nov-09
% Mariane Goncalves Kulik (mgkulik) - 2018-nov-09
% UFPR Bioinformatics team - http://www.bioinfo.ufpr.br/

[n m] = size(xseq);

% Execute matrix transformations to ensure equal values get the same aa
% numeric index
if defSize == 20
    vls = double(aa2int(upper(xseq)))-1;
else
    vls = double(nt2int(upper(xseq)))-1;
end
pot = repmat(double(0:(m-1)),n,1);
mret = sum((repmat(defSize,n,m).^pot).*vls,2)+1;

% Ensures no indexing error
inot = find(~prod(1-double(vls<0 | vls>defSize-1),2));
mret(inot) = -1;



function mret = ij2inds(ijs, tmcol)
% Converts the subscript to the matrix index using as subscript the amino
% acid numeric information

% Args:
%       ijs: Values calculated from the amino acid numeric representation
%       tmcol: Number of rows calculated from to the fist portion of the
%       mask

% Returns:
%       mret: The matrix index where the present values will be placed

% Mariane Gonçalves Kulik (mgkulik) - 2018-nov-09
% UFPR Bioinformatics team - http://www.bioinfo.ufpr.br/

mret = (ijs(:,2)-1)*tmcol + (ijs(:,1));


function [vals, inds] = unifysameinds(ind0, x, func)
% Gets a vector of indexes and a vector of values with different sizes,
% equalizes them and group the ones with same id to apply the informed 
% function.

% Args:
%       ind0: The ids vector (1xN)
%       x: The values vector (1xN). Usually it will have a different
%               size from the ids vector.
%       func: The function to apply in the x values

% Returns:
%       vals: The values transformed by the function;
%       inds: The new ids, with the duplicates removed.

% Roberto Tadeu Raittz (rraittz) - 2018-nov-09
% UFPR Bioinformatics team - http://www.bioinfo.ufpr.br

itn = vectincol(x)';
ids = vectincol(ind0)';

[idsord, ind] = sort(ids);
itnord = itn(ind);
[uidx, ii] = unique(idsord);

%list = [(1+ii-[1;diff(ii)])' ii'];
list = [(1+ii-[1;diff(ii)]) ii];
list(1,1) = 1;

iun = (cellfun(@(x) x(1):x(end),mat2celllines(list),'UniformOutput',false));
vals = cell2mat(cellfun(@(x) func(itnord(x)),iun,'UniformOutput',false));
inds = uidx;


function mret = mat2celllines(mat)
% Group all the columns informed in a matrix in one cell column.

% Args:
%       mat: The matrix to transform

% Returns:
%       mret: The Nx1 cell array

% Mariane Gonçalves Kulik (mgkulik) - 2018-nov-09
% UFPR Bioinformatics team - http://www.bioinfo.ufpr.br

mret = mat2cell(mat,ones(1,length(mat(:,1))),length(mat(1,:)));