%%% INPUT
%%%    filename_theta, filename_phi1, filename_phi2
%%%    layer_data
%%%    pair
%%%    m, N, h
%%%    functions TH, PH1, PH2
%%% OUTPUT
%%%    output files for theta, phi1, phi2

% normalized coordinate of i-th node
X = @(i) i * h;   % i = 0..N

L1 = layer_data(1).L;
L2 = layer_data(2).L;

ftheta = fopen(filename_theta, "wt");
for i = 0 : N
  if (i < m)
    tau = L1 * X(i);
  else
    tau = L1 + (X(i) - 1) * L2;
  endif
  fprintf(ftheta, "%.12f %.12f\n", tau, pair(TH(i)));  
endfor
fclose(ftheta);

fphi1 = fopen(filename_phi1, "wt");
for i = 0 : m
  tau = L1 * X(i);
  fprintf(fphi1, "%.12f %.12f\n", tau, pair(PH1(i)));
endfor
fclose(fphi1);

fphi2 = fopen(filename_phi2, "wt");
for i = m : N
  tau = L1 + (X(i) - 1) * L2;
  fprintf(fphi2, "%.12f %.12f\n", tau, pair(PH2(i)));
endfor
fclose(fphi2);
