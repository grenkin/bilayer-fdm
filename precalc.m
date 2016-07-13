%%% INPUT
%%%    layer_data, boundary_data
%%% OUTPUT
%%%    K1, K2, B1, B2, C1, C2
%%%    tilde_gamma1, tilde_gamma2, thetab1, thetab2
%%%    G (if n1 != n2)
%%%    EQUAL_N

% refractive indices for boundaries
for j = 1:2
  boundary_data(j).n = layer_data(j).n;
endfor

EQUAL_N = (layer_data(1).n == layer_data(2).n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions for calculating normalized model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function retval = gamma_func(emiss)
  retval = emiss / ( 2 * (2 - emiss) );
endfunction

function retval = calc_layer(l)
  alpha = 1 / (3 - l.A * l.omega);
  D = l.n ^ 2 * alpha;
  retval.B = D / l.L;
  retval.C = l.Nc * l.n ^ 2 / l.L;
  retval.K = l.L * l.n ^ 2 * (1 - l.omega);
endfunction

function retval = calc_boundary(b)
  emiss = 1 - b.R;
  gam = gamma_func(emiss);
  retval.tilde_gamma = b.n ^ 2 * gam;
  retval.thetab = b.thetab;
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% model parameters calculating
for j = 1:2
  layer(j) = calc_layer(layer_data(j));
  boundary(j) = calc_boundary(boundary_data(j));
endfor

% define scalar variables for the convenience when solving a linear system
layer_params_list = {"K", "B", "C"};
boundary_params_list = {"tilde_gamma", "thetab"};

for j = 1:2
  for iter = layer_params_list
    eval([iter{} int2str(j) " = " num2str(layer(j).(iter{}))]);    
  endfor
  for iter = boundary_params_list
    eval([iter{} int2str(j) " = " num2str(boundary(j).(iter{}))]);    
  endfor
endfor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions for calculating coefficient G
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function retval = psi_ij(mu, nij)
  radicand = 1 - nij ^ 2 * (1 - mu ^ 2);
  if (radicand > 0)
    retval = sqrt(radicand);
  else
    retval = 0;
  endif
endfunction

function retval = R_ij(mu, nij)
  if (mu == 0)
    retval = 1;
  else
    frac1 = (psi_ij(mu, nij) - nij * mu) / (psi_ij(mu, nij) + nij * mu);
    frac2 = (nij * psi_ij(mu, nij) - mu) / (nij * psi_ij(mu, nij) + mu);
    retval = 0.5 * (frac1 ^ 2 + frac2 ^ 2);
  endif
endfunction

function retfunc = integrand1(n1, n2)
  n12 = n1 / n2;
  retfunc = @(mu_vect) arrayfun(@(mu)(mu * (1 - R_ij(mu, n12))), mu_vect);
endfunction

function retfunc = integrand2(n1, n2)
  n12 = n1 / n2;
  n21 = n2 / n1;
  retfunc = @(mu_vect) arrayfun(@(mu)(mu ^ 2 * (R_ij(mu, n12) + R_ij(mu, n21))), mu_vect);
endfunction

function retval = calcG(n1, n2)
  int1 = quadcc(integrand1(n1, n2), 0, 1);
  int2 = quadcc(integrand2(n1, n2), 0, 1);
  retval = (n1 ^ 2 * int1) / (3 * int2);
endfunction

% calculating coefficient G
if (!EQUAL_N)
  G = calcG(layer_data(1).n, layer_data(2).n)
endif
