pre0

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Single layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%

LEN = 1;  % layer thickness

% data for layer 1
layer_data(1) = struct(
  "L",     LEN / 2,
  "omega", 0.9,
  "n",     1,
  "Nc",    1/600,
  "A",     1);
% data for the left boundary
boundary_data(1) = struct(
  "R",      0.3,
  "thetab", 1);
% data for layer 2
layer_data(2) = layer_data(1);
% data for the right boundary
boundary_data(2) = struct(
  "R",      0.3,
  "thetab", 0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%

calc_none
filename_theta = "output_single/theta.dat";
filename_phi1 = "output_single/phi1.dat";
filename_phi2 = "output_single/phi2.dat";
print_theta_phi
theta_none = theta;
phi1_none = phi1;
phi2_none = phi2;

% make plots

filename_plot_theta = "output_single/theta.png";
filename_plot_phi = "output_single/phi.png";
plot_none