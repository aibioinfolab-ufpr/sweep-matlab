function mret = primon(n)
% Returns all the first n primes, starting by 2.

% Args:
%       n: The expected number of primes

% Returns:
%       mret: The list of n first primes

% Roberto Tadeu Raittz (rraittz) - 2018-nov-06
% Mariane Goncalves Kulik (mgkulik) - 2018-nov-06
% UFPR Bioinformatics team - http://www.bioinfo.ufpr.br/

m = n;
nmax = 100;
ps = find(isprime(1:nmax));

while length(ps) < m
    nmax = nmax*2;
    ps = find(isprime(1:nmax));
end

mret = ps(1:m);