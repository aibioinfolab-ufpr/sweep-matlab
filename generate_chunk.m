function idxMat = generate_chunk(chunks, file_len)
% Calculates the indexes according to original file size and ammount of
% chunks informed to particionate the fasta files in memory.

% Args:
%       chunks: Ammount of chunks to break the file
%       file_len: Total of rows of the file

% Returns:
%       idxMat: Matrix with the start:end id for each chunk

% Mariane Goncalves Kulik (mgkulik) - 2018-nov-06
% UFPR Bioinformatics team - http://www.bioinfo.ufpr.br/

z = floor(file_len/chunks);
idxMat = repmat(z, chunks, 2) .* [((1:chunks)-1)' (1:chunks)'] + [ones(chunks,1) zeros(chunks,1)];
idxMat(end) = file_len;