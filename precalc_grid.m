%%% OUTPUT
%%%    m, N, h
%%%    functions TH, PH1, PH2
%%%    unknowns

% grid parameters
m = 1000;
N = 2 * m;
X2 = 2;
h = X2 / N;

% indices in the matrix of a linear system
TH = @(i) i + 1;   % i = 0..N
PH1 = @(i) (N + 1) + i + 1;  % i = 0..m
PH2 = @(i) (N + 1) + i + 2;  % i = m..N

unknowns = (N + 1) + (N + 2);
