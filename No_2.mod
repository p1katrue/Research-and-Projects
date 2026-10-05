% No_2.mod
%
% DESCRIPTION: Question 2 - Comparison of Exogenous (#1) vs Endogenous (#2) Foreign Block
%
%----------------------------------------------------------------
% 1. Defining variables (in LOG DEVIATION) and parameters
%----------------------------------------------------------------
var c_hat i_hat pi_hat y_hat S_hat q_hat pi_H_hat pi_F_hat mc_hat 
    a_hat e_c_hat ystar_hat pistar_hat istar_hat phi_hat 
    eps_z_hat eps_H_hat eps_ystar_hat eps_pistar_hat eps_istar_hat 
    CA_hat;

varexo eps_m_shock eta_phi_shock eta_z_shock eta_H_shock 
       eta_ystar_shock eta_pistar_shock eta_istar_shock;

parameters sig alpha eta bet theta_H chi varphi h gam delta_H 
           rho_i psi_pi psi_y psi_dy psi_e 
           rho_phi rho_z rho_H 
           bet_star rho_i_star kappa_star psi_pi_star psi_y_star 
           rho_eps_ystar rho_eps_pistar;

%----------------------------------------------------------------
% 2. Calibration (of parameters)
%----------------------------------------------------------------
sig = 1;
alpha = 0.30;
eta = 0.85;
bet = 0.99;
theta_H = 0.75;
chi = 0.01;
varphi = 1.26;
h = 0.25;
gam = 0.15;
delta_H = 0.30;

rho_i = 0.80;
psi_pi = 1.90;
psi_y = 0.05;
psi_dy = 0.55;
psi_e = 0; 

rho_phi = 0.70;
rho_z = 0.70;
rho_H = 0.70;

bet_star = 0.99;
rho_i_star = 0.70;
kappa_star = 0.05;
psi_pi_star = 1.50;
psi_y_star = 0.125;
rho_eps_ystar = 0.70;
rho_eps_pistar = 0.70;

%----------------------------------------------------------------
% 3. Model Block
%----------------------------------------------------------------
model (linear);
    c_hat = (h/(1+h))*c_hat(-1) + (1/(1+h))*c_hat(+1) - (1/sig)*((1-h)/(1+h))*(i_hat - pi_hat(+1)) + gam*eps_z_hat;
    y_hat = (1-alpha)*c_hat + alpha*ystar_hat + alpha*eta*(2-alpha)*S_hat;
    q_hat = (1-alpha)*S_hat;
    S_hat - S_hat(-1) = pi_F_hat - pi_H_hat;
    pi_H_hat - delta_H*pi_H_hat(-1) = bet*(pi_H_hat(+1) - delta_H*pi_H_hat) + ((1-theta_H)*(1-theta_H*bet)/theta_H)*mc_hat + eps_H_hat;
    mc_hat = varphi*y_hat + alpha*S_hat + sig*c_hat;
    pi_hat = pi_H_hat + alpha*(S_hat - S_hat(-1));
    i_hat - istar_hat = e_c_hat(+1) - chi*(a_hat + phi_hat(+1));
    y_hat - c_hat = a_hat - (1/bet)*a_hat(-1) + (alpha/(1-alpha))*q_hat;
    pi_F_hat = e_c_hat + pistar_hat;
    i_hat = rho_i*i_hat(-1) + psi_pi*pi_hat + psi_y*y_hat + psi_dy*(y_hat - y_hat(-1)) + psi_e*e_c_hat - eps_m_shock;
    phi_hat = rho_phi*phi_hat(-1) - eta_phi_shock;
    eps_z_hat = rho_z*eps_z_hat(-1) + eta_z_shock;
    eps_H_hat = rho_H*eps_H_hat(-1) + eta_H_shock;
    
    % Foreign Block Equations 
    ystar_hat = ystar_hat(+1) - (istar_hat - pistar_hat(+1)) + eps_ystar_hat;
    pistar_hat = bet_star*pistar_hat(+1) + kappa_star*ystar_hat + eps_pistar_hat;
    istar_hat = rho_i_star*istar_hat(-1) + (1-rho_i_star)*(psi_pi_star*pistar_hat + psi_y_star*ystar_hat) + eps_istar_hat;
    eps_ystar_hat = rho_eps_ystar*eps_ystar_hat(-1) + eta_ystar_shock;
    eps_pistar_hat = rho_eps_pistar*eps_pistar_hat(-1) + eta_pistar_shock;
    eps_istar_hat = eta_istar_shock;
    
    CA_hat = y_hat - c_hat;
end;

%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------
shocks;
    var eta_ystar_shock = 1^2; 
end;

check;

stoch_simul(irf = 12, order = 1, nograph, irf_shocks=(eta_ystar_shock)) 
    y_hat c_hat i_hat pi_hat e_c_hat q_hat CA_hat;

irf_y_q2 = oo_.irfs.y_hat_eta_ystar_shock;
irf_c_q2 = oo_.irfs.c_hat_eta_ystar_shock;
irf_i_q2 = oo_.irfs.i_hat_eta_ystar_shock;
irf_pi_q2 = oo_.irfs.pi_hat_eta_ystar_shock;
irf_ec_q2 = oo_.irfs.e_c_hat_eta_ystar_shock;
irf_q_q2 = oo_.irfs.q_hat_eta_ystar_shock;
irf_e_q2 = cumsum(oo_.irfs.e_c_hat_eta_ystar_shock);
irf_ca_q2 = oo_.irfs.CA_hat_eta_ystar_shock;

irf_y_q1c   = [0.1997, 0.1169, 0.0648, 0.0336, 0.0152, 0.0042, -0.0026, -0.0068, -0.0094, -0.0110, -0.0120, -0.0124];
irf_c_q1c   = [0.0178, 0.0239, 0.0301, 0.0372, 0.0439, 0.0491, 0.0526, 0.0545, 0.0552, 0.0550, 0.0542, 0.0529];
irf_i_q1c   = [0.0234, 0.0193, 0.0149, 0.0107, 0.0071, 0.0043, 0.0023, 0.0009, -0.0001, -0.0008, -0.0012, -0.0015];
irf_pi_q1c  = [-0.0508, 0.0212, 0.0131, 0.0075, 0.0042, 0.0024, 0.0014, 0.0009, 0.0006, 0.0004, 0.0003, 0.0002];
irf_ec_q1c  = [-0.2329, 0.0260, 0.0236, 0.0203, 0.0168, 0.0136, 0.0110, 0.0090, 0.0075, 0.0064, 0.0056, 0.0049];
irf_q_q1c   = [-0.1821, -0.1773, -0.1668, -0.1539, -0.1413, -0.1301, -0.1205, -0.1124, -0.1055, -0.0995, -0.0942, -0.0895];
irf_e_q1c   = [-0.2329, -0.2069, -0.1832, -0.1629, -0.1461, -0.1325, -0.1216, -0.1126, -0.1051, -0.0987, -0.0931, -0.0881];
irf_ca_q1c  = [0.1819, 0.0930, 0.0347, -0.0036, -0.0286, -0.0449, -0.0551, -0.0613, -0.0646, -0.0660, -0.0661, -0.0654];

%----------------------------------------------------------------
% 5. Plotting 
%----------------------------------------------------------------
figure('Name', 'Question 2: Exogenous vs Endogenous Comparison', 'Color', 'w');
t_axis = 0:(options_.irf - 1); 

% 1. Output
subplot(4,2,1);
plot(t_axis, irf_y_q1c, 'k-', 'LineWidth', 1.5); hold on;
plot(t_axis, irf_y_q2, 'k--', 'LineWidth', 1.5); hold off;
title('Output'); xlabel('quarter'); ylabel('% deviation'); yline(0, 'r--'); xlim([0 11]);
legend('#1', '#2', 'Location', 'southeast');

% 2. Consumption
subplot(4,2,2);
plot(t_axis, irf_c_q1c, 'k-', 'LineWidth', 1.5); hold on;
plot(t_axis, irf_c_q2, 'k--', 'LineWidth', 1.5); hold off;
title('Consumption'); xlabel('quarter'); ylabel('% deviation'); yline(0, 'r--'); xlim([0 11]);
legend('#1', '#2', 'Location', 'southeast');

% 3. Nominal Interest Rate
subplot(4,2,3);
plot(t_axis, irf_i_q1c, 'k-', 'LineWidth', 1.5); hold on;
plot(t_axis, irf_i_q2, 'k--', 'LineWidth', 1.5); hold off;
title('Nominal int. rate'); xlabel('quarter'); ylabel('deviation (%)'); yline(0, 'r--'); xlim([0 11]);
legend('#1', '#2', 'Location', 'southeast');

% 4. CPI Inflation
subplot(4,2,4);
plot(t_axis, irf_pi_q1c, 'k-', 'LineWidth', 1.5); hold on;
plot(t_axis, irf_pi_q2, 'k--', 'LineWidth', 1.5); hold off;
title('CPI Inflation'); xlabel('quarter'); ylabel('deviation (%)'); yline(0, 'r--'); xlim([0 11]);
legend('#1', '#2', 'Location', 'southeast');

% 5. Nominal Depreciation Rate
subplot(4,2,5);
plot(t_axis, irf_ec_q1c, 'k-', 'LineWidth', 1.5); hold on;
plot(t_axis, irf_ec_q2, 'k--', 'LineWidth', 1.5); hold off;
title('Nom. XR dep. rate'); xlabel('quarter'); ylabel('deviation (%)'); yline(0, 'r--'); xlim([0 11]);
legend('#1', '#2', 'Location', 'southeast');

% 6. Real Exchange Rate
subplot(4,2,6);
plot(t_axis, irf_q_q1c, 'k-', 'LineWidth', 1.5); hold on;
plot(t_axis, irf_q_q2, 'k--', 'LineWidth', 1.5); hold off;
title('Real XR'); xlabel('quarter'); ylabel('% deviation'); yline(0, 'r--'); xlim([0 11]);
legend('#1', '#2', 'Location', 'southeast');

% 7. Nominal Exchange Rate Level
subplot(4,2,7);
plot(t_axis, irf_e_q1c, 'k-', 'LineWidth', 1.5); hold on;
plot(t_axis, irf_e_q2, 'k--', 'LineWidth', 1.5); hold off;
title('Nom. XR (level)'); xlabel('quarter'); ylabel('% deviation'); yline(0, 'r--'); xlim([0 11]);
legend('#1', '#2', 'Location', 'southeast');

% 8. Current Account
subplot(4,2,8);
plot(t_axis, irf_ca_q1c, 'k-', 'LineWidth', 1.5); hold on;
plot(t_axis, irf_ca_q2, 'k--', 'LineWidth', 1.5); hold off;
title('Current account'); xlabel('quarter'); ylabel('% deviation'); yline(0, 'r--'); xlim([0 11]);
legend('#1', '#2', 'Location', 'southeast');

sgtitle('Q2: Comparison between Exogenous (#1) and Endogenous (#2) Foreign Block');