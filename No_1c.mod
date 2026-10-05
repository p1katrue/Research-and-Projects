% Q1c.mod
%
% DESCRIPTION: dynare.mod file to solve and simulate the SOE model
%             Question 1(c): Foreign-output shock
%
%----------------------------------------------------------------
% 1. Defining variables (in LOG DEVIATION) and parameters
%----------------------------------------------------------------
var c_hat i_hat pi_hat y_hat S_hat q_hat pi_H_hat pi_F_hat mc_hat 
    a_hat e_c_hat ystar_hat pistar_hat istar_hat phi_hat 
    eps_z_hat eps_H_hat CA_hat;
varexo eps_m_shock eta_phi_shock eta_z_shock eta_H_shock 
       eta_ystar_shock eta_pistar_shock eta_istar_shock;
parameters sig alpha eta bet theta_H chi varphi h gam delta_H 
           rho_i psi_pi psi_y psi_dy psi_e 
           rho_phi rho_z rho_H rho_ystar rho_pistar rho_istar;
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
rho_ystar = 0.70;
rho_pistar = 0.70;
rho_istar = 0.70;
%----------------------------------------------------------------
% 3. Model (Log-Linearized Model)
%----------------------------------------------------------------
model (linear);
    % (1) Consumption Euler-equation (the IS equation)
    c_hat = (h/(1+h))*c_hat(-1) + (1/(1+h))*c_hat(+1) - (1/sig)*((1-h)/(1+h))*(i_hat - pi_hat(+1)) + gam*eps_z_hat;
    
    % (2) Goods-market clearing condition
    y_hat = (1-alpha)*c_hat + alpha*ystar_hat + alpha*eta*(2-alpha)*S_hat;
    
    % (3) Terms of trade and real exchange rate
    q_hat = (1-alpha)*S_hat;
    
    % (4) Changes of the terms of trade
    S_hat - S_hat(-1) = pi_F_hat - pi_H_hat;
    
    % (5) Domestic-price inflation (Phillips curve)
    pi_H_hat - delta_H*pi_H_hat(-1) = bet*(pi_H_hat(+1) - delta_H*pi_H_hat) + ((1-theta_H)*(1-theta_H*bet)/theta_H)*mc_hat + eps_H_hat;
    
    % (6) The real marginal cost
    mc_hat = varphi*y_hat + alpha*S_hat + sig*c_hat;
    
    % (7) The wedge between CPI- and PPI-inflation
    pi_hat = pi_H_hat + alpha*(S_hat - S_hat(-1));
    
    % (8) The uncovered interest-parity (UIP) condition
    i_hat - istar_hat = e_c_hat(+1) - chi*(a_hat + phi_hat(+1));
    
    % (9) The net-foreign-assets accumulation
    y_hat - c_hat = a_hat - (1/bet)*a_hat(-1) + (alpha/(1-alpha))*q_hat;
    
    % (10) Imported-good inflation (Law of one price)
    pi_F_hat = e_c_hat + pistar_hat;
    
    % (11) Monetary-policy (Taylor) rule
    i_hat = rho_i*i_hat(-1) + psi_pi*pi_hat + psi_y*y_hat + psi_dy*(y_hat - y_hat(-1)) + psi_e*e_c_hat - eps_m_shock;
    
    % (12) Evolution of risk premium
    phi_hat = rho_phi*phi_hat(-1) - eta_phi_shock;
    
    % (13) Evolution of preference shock
    eps_z_hat = rho_z*eps_z_hat(-1) + eta_z_shock;
    
    % (14) Evolution of cost-push shock
    eps_H_hat = rho_H*eps_H_hat(-1) + eta_H_shock;
    
    % (15.A) Foreign output
    ystar_hat = rho_ystar*ystar_hat(-1) + eta_ystar_shock;
    
    % (16.A) Foreign inflation
    pistar_hat = rho_pistar*pistar_hat(-1) + eta_pistar_shock;
    
    % (17.A) Foreign interest rate
    istar_hat = rho_istar*istar_hat(-1) + eta_istar_shock;
    
    % Current Account
    CA_hat = y_hat - c_hat;
end;
%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------
shocks;
    var eta_ystar_shock = 1^2;
end;
check;
stoch_simul(irf = 12, order = 1, nograph, irf_shocks=(eta_ystar_shock)) y_hat c_hat i_hat pi_hat e_c_hat q_hat CA_hat;
e_level_irf = cumsum(oo_.irfs.e_c_hat_eta_ystar_shock);
%----------------------------------------------------------------
% 5. Plotting 
%----------------------------------------------------------------
figure('Name', 'Question 1(c): Foreign-Output Shock', 'Color', 'w');
t_axis = 0:(options_.irf - 1); 

subplot(4,2,1);
plot(t_axis, oo_.irfs.y_hat_eta_ystar_shock, 'LineWidth', 1.5);
title('Output');
xlabel('quarter');
ylabel('% deviation');
yline(0, 'r--');
xlim([0 11]);

subplot(4,2,2);
plot(t_axis, oo_.irfs.c_hat_eta_ystar_shock, 'LineWidth', 1.5);
title('Consumption');
xlabel('quarter');
ylabel('% deviation');
yline(0, 'r--');
xlim([0 11]);

subplot(4,2,3);
plot(t_axis, oo_.irfs.i_hat_eta_ystar_shock, 'LineWidth', 1.5);
title('Nominal int. rate');
xlabel('quarter');
ylabel('deviation (%)');
yline(0, 'r--');
xlim([0 11]);

subplot(4,2,4);
plot(t_axis, oo_.irfs.pi_hat_eta_ystar_shock, 'LineWidth', 1.5);
title('CPI Inflation');
xlabel('quarter');
ylabel('deviation (%)');
yline(0, 'r--');
xlim([0 11]);

subplot(4,2,5);
plot(t_axis, oo_.irfs.e_c_hat_eta_ystar_shock, 'LineWidth', 1.5);
title('Nom. XR dep. rate');
xlabel('quarter');
ylabel('deviation (%)');
yline(0, 'r--');
xlim([0 11]);

subplot(4,2,6);
plot(t_axis, oo_.irfs.q_hat_eta_ystar_shock, 'LineWidth', 1.5);
title('Real XR');
xlabel('quarter');
ylabel('% deviation');
yline(0, 'r--');
xlim([0 11]);

subplot(4,2,7);
plot(t_axis, e_level_irf, 'LineWidth', 1.5);
title('Nom. XR (level)');
xlabel('quarter');
ylabel('% deviation');
yline(0, 'r--');
xlim([0 11]);

subplot(4,2,8);
plot(t_axis, oo_.irfs.CA_hat_eta_ystar_shock, 'LineWidth', 1.5);
title('Current account');
xlabel('quarter');
ylabel('% deviation');
yline(0, 'r--');
xlim([0 11]);

sgtitle('IRFs to a positive 1% Foreign-Output Shock');


