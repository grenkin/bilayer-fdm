pre0

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two layers
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

EQUAL_N = (layer_data(1).n == layer_data(2).n);

if (!EQUAL_N)
  calc_fresnel
  filename_theta = "output_bilayer/theta_fresnel.dat";
  filename_phi1 = "output_bilayer/phi1_fresnel.dat";
  filename_phi2 = "output_bilayer/phi2_fresnel.dat";
  print_theta_phi
  theta_fresnel = theta;
  phi1_fresnel = phi1;
  phi2_fresnel = phi2;
endif

calc_none
filename_theta = "output_bilayer/theta_none.dat";
filename_phi1 = "output_bilayer/phi1_none.dat";
filename_phi2 = "output_bilayer/phi2_none.dat";
print_theta_phi
theta_none = theta;
phi1_none = phi1;
phi2_none = phi2;

% make plots

if (!EQUAL_N)
  filename_plot_theta = "output_bilayer/theta_all.png";
  filename_plot_phi = "output_bilayer/phi_all.png";
  plot_fresnel_none  
else
  filename_plot_theta = "output_bilayer/theta_none.png";
  filename_plot_phi = "output_bilayer/phi_none.png";
  plot_none
endif
