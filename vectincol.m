function vret = vectincol(vec)
% Certifies the vector informed is shaped as 1xN

% Args:
%       vec: The vector to reshape

% Returns:
%        vret: The vector shaped as 1xN

% Mariane Gonçalves Kulik (mgkulik) - 2018-nov-09
% UFPR Bioinformatics team - http://www.bioinfo.ufpr.br

vret = reshape(vec, [1,length(vec)]);