%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculation the solution of BVP with Fresnel
%% conjugation conditions at the interface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Required: n1 != n2 (so that G is defined)
%%% INPUT
%%%    K1, K2, B1, B2, C1, C2
%%%    tilde_gamma1, tilde_gamma2, thetab1, thetab2
%%%    G
%%% OUTPUT
%%%    pair, theta, phi1, phi2

precalc
assert(!EQUAL_N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculations by finite difference method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

precalc_grid

A = sparse(unknowns, unknowns);
rhs = zeros(unknowns, 1);

% initial guess for Newton's method
thetaold = thetab1 + (thetab2 - thetab1) * (0:N)' / N;

% Newton's method

eps_iter = 1e-7;  % epsilon in the stopping condition
num_iter = 0;

disp("\nNewton's method");
while (1)
  disp(++num_iter);

  A(TH(0), TH(0)) = 1;
  rhs(TH(0)) = thetab1;

  A(PH1(0), PH1(0)) = B1 / h + tilde_gamma1 + h / 2 * K1;
  A(PH1(0), PH1(1)) = - B1 / h;
  A(PH1(0), TH(0)) = - h / 2 * K1 * (4 * thetaold(TH(0)) ^ 3);
  rhs(PH1(0)) = tilde_gamma1 * thetab1 ^ 4 - h / 2 * K1 * (3 * thetaold(TH(0)) ^ 4);

  for i = 1 : m - 1
    A(TH(i), TH(i - 1)) = - C1 / h ^ 2;
    A(TH(i), TH(i)) = 2 * C1 / h ^ 2 + K1 * (4 * thetaold(TH(i)) ^ 3);
    A(TH(i), TH(i + 1)) = - C1 / h ^ 2;
    A(TH(i), PH1(i)) = - K1;
    rhs(TH(i)) = K1 * (3 * thetaold(TH(i)) ^ 4);

    A(PH1(i), PH1(i - 1)) = - B1 / h ^ 2;
    A(PH1(i), PH1(i)) = 2 * B1 / h ^ 2 + K1;
    A(PH1(i), PH1(i + 1)) = - B1 / h ^ 2;  
    A(PH1(i), TH(i)) = - K1 * (4 * thetaold(TH(i)) ^ 3);
    rhs(PH1(i)) = - K1 * (3 * thetaold(TH(i)) ^ 4);
  endfor

  A(TH(m), TH(m - 1)) = - C1 / h;
  A(TH(m), TH(m)) = C1 / h + C2 / h + h / 2 * K1 * (4 * thetaold(TH(m)) ^ 3) ...
    + h / 2 * K2 * (4 * thetaold(TH(m)) ^ 3);
  A(TH(m), TH(m + 1)) = - C2 / h;
  A(TH(m), PH1(m)) = - h / 2 * K1;
  A(TH(m), PH2(m)) = - h / 2 * K2;
  rhs(TH(m)) = h / 2 * K1 * (3 * thetaold(TH(m)) ^ 4) ...
    + h / 2 * K2 * (3 * thetaold(TH(m)) ^ 4);

  A(PH1(m), PH1(m - 1)) = - B1 / h;
  A(PH1(m), PH1(m)) = G + B1 / h + h / 2 * K1;
  A(PH1(m), PH2(m)) = - G;
  A(PH1(m), TH(m)) = - h / 2 * K1 * (4 * thetaold(TH(m)) ^ 3);
  rhs(PH1(m)) = - h / 2 * K1 * (3 * thetaold(TH(m)) ^ 4);

  A(PH2(m), PH1(m)) = - G;
  A(PH2(m), PH2(m)) = G + B2 / h + h / 2 * K2;
  A(PH2(m), PH2(m + 1)) = - B2 / h;
  A(PH2(m), TH(m)) = - h / 2 * K2 * (4 * thetaold(TH(m)) ^ 3);
  rhs(PH2(m)) = - h / 2 * K2 * (3 * thetaold(TH(m)) ^ 4);

  for i = m + 1 : N - 1
    A(TH(i), TH(i - 1)) = - C2 / h ^ 2;
    A(TH(i), TH(i)) = 2 * C2 / h ^ 2 + K2 * (4 * thetaold(TH(i)) ^ 3);
    A(TH(i), TH(i + 1)) = - C2 / h ^ 2;
    A(TH(i), PH2(i)) = - K2;
    rhs(TH(i)) = K2 * (3 * thetaold(TH(i)) ^ 4);

    A(PH2(i), PH2(i - 1)) = - B2 / h ^ 2;
    A(PH2(i), PH2(i)) = 2 * B2 / h ^ 2 + K2;
    A(PH2(i), PH2(i + 1)) = - B2 / h ^ 2;
    A(PH2(i), TH(i)) = - K2 * (4 * thetaold(TH(i)) ^ 3);
    rhs(PH2(i)) = - K2 * (3 * thetaold(TH(i)) ^ 4);
  endfor

  A(TH(N), TH(N)) = 1;
  rhs(TH(N)) = thetab2;
                                                                              
  A(PH2(N), PH2(N - 1)) = - B2 / h;
  A(PH2(N), PH2(N)) = B2 / h + tilde_gamma2 + h / 2 * K2;
  A(PH2(N), TH(N)) = - h / 2 * K2 * (4 * thetaold(TH(N)) ^ 3);
  rhs(PH2(N)) = tilde_gamma2 * thetab2 ^ 4 - h / 2 * K2 * (3 * thetaold(TH(N)) ^ 4);

  pair = A \ rhs;
  
  theta = pair(TH(0) : TH(N));

  diff_theta = norm(theta - thetaold)
  if (diff_theta < eps_iter)
    break;
  endif

  thetaold = theta;

endwhile

phi1 = pair(PH1(0) : PH1(m));
phi2 = pair(PH2(m) : PH2(N));
