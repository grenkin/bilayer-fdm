%%% INPUT
%%%    filename_plot_theta, filename_plot_phi
%%%    layer_data
%%%    m, N, h
%%%    theta_none, phi1_none, phi2_none
%%% OUTPUT
%%%    plots for theta, phi

figure
x_total = [(0 : m - 1) * h * layer_data(1).L, layer_data(1).L + ((m : N) * h - 1) * layer_data(2).L];
plot(x_total, theta_none, "color", [0 .6 0], "linewidth", 2);
xlabel("tau");
ylabel("theta");
saveas(1, filename_plot_theta);  

figure
x_1 = (0 : m) * h * layer_data(1).L;
x_2 = layer_data(1).L + ((m : N) * h - 1) * layer_data(2).L;
plot(x_1, phi1_none, "color", [0 .6 0], "linewidth", 2, ...
  x_2, phi2_none, "color", [0 .6 0], "linewidth", 2);
xlabel("tau");
ylabel("phi");
saveas(2, filename_plot_phi);
