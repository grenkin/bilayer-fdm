pre0

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two layers with various refractive indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data for layer 1
layer_data(1) = struct(
  "L",     1,
  "omega", 0.2,
  "n",     1.8,
  "Nc",    0.001,
  "A",       0);
% data for the left boundary
boundary_data(1) = struct(
  "R",      0.3,
  "thetab",   1);

% data for layer 2
layer_data(2) = struct(
  "L",       1,
  "omega", 0.9,
  "n",       1,
  "Nc",   0.0001,
  "A",       0);
% data for the right boundary
boundary_data(2) = struct(
  "R",      0.2,
  "thetab", 0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%

n1_range = 1.1 : 0.1 : 2;
cnt = 0;
for n1_value = n1_range
  layer_data(1).n = n1_value;
  calc_fresnel
  theta_fresnel = theta;
  phi1_fresnel = phi1;
  phi2_fresnel = phi2;
  calc_none
  theta_none = theta;
  phi1_none = phi1;
  phi2_none = phi2;
   
  rms(++cnt) = sqrt(trapz((0 : N)' * h, (theta_fresnel .- theta_none) .^ 2));
  % it works for equal L1 and L2 
endfor

frms = fopen("output_rms/rms.dat", "wt");
for i = 1 : cnt
  fprintf(frms, "%.12f %.12f\n", n1_range(i), rms(i));
endfor
fclose(frms);

figure
plot(n1_range, rms, "b", "linewidth", 2);
xlabel("n1");
ylabel("rms");
saveas(1, "output_rms/rms.png");
