function check_fasta(seq, type)
% Checks if the informed fasta contains aminoacids.
% Returns an error if it contains nucleotides.

notNT = [6,7,10,11,14,15];
int_valsAA = unique(aa2int(seq));
isNT = isempty(intersect(notNT, int_valsAA));
if strcmpi(type, 'AA')
    text_val = 'NTs';
else
    text_val = 'AAs';
end

if (isNT && strcmpi(type, 'AA')||~isNT && strcmpi(type, 'NT'))
    message = strcat('The "fasta type" parameter is setted to', {' '}, type, ',', ...
        {' '}, 'but the fasta has', {' '}, text_val, '. Please review the parameters and the fasta file.');
    generate_log(message, 2);
    error(char(message));
end