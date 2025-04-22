%% ASSIGNMENT 1 A - CANTILEVER BEAM

clear all
close all
clc

%% Natural frequencies and mode shapes

L = 1.2; % [m]
h = 0.008; % [m]
b = 0.04; % [m]
rho = 2700; % [kg/m^3]
E = 68e9; % [Pa]
m = rho*b*h; % [kg/m]
J = b*h^3/12;

f = linspace(0.01,200,1e6);
omega = 2*pi*f;
gamma_quarta = m*omega.^2/(E*J);
gamma = sqrt(sqrt(gamma_quarta));

H = @(omega) [1                         0                         1                         0;...
              0                         1                         0                         1;...
              -cos(sqrt(sqrt(m*omega.^2/(E*J)))*L)  -sin(sqrt(sqrt(m*omega.^2/(E*J)))*L)  cosh(sqrt(sqrt(m*omega.^2/(E*J)))*L)  sinh(sqrt(sqrt(m*omega.^2/(E*J)))*L);...
              sin(sqrt(sqrt(m*omega.^2/(E*J)))*L)   -cos(sqrt(sqrt(m*omega.^2/(E*J)))*L)  sinh(sqrt(sqrt(m*omega.^2/(E*J)))*L)  cosh(sqrt(sqrt(m*omega.^2/(E*J)))*L)];

for i=1:length(omega)
    dets(i) = det(H(omega(i)));
end

figure
box on
semilogy(f,abs(dets),'-b')
hold on
grid on
xlabel('f [Hz]')
%xlim([0 8])

i_nat =[];
for i=2:length(dets)-1
    if abs(dets(i)) < abs(dets(i-1)) && abs(dets(i)) < abs(dets(i+1))
        i_nat(end+1) = i;
    end
end

fprintf('Natural frequencies: \n%f\n%f\n%f\n%f\n',f(i_nat))
fprintf('Omega: \n%f\n%f\n%f\n%f\n',omega(i_nat))
plot(f(i_nat),abs(dets(i_nat)),'or')

% Reduced system

for i_mode=1:length(i_nat)
    fprintf('MODE %i \n',i_mode)
    omega_i=omega(i_nat(i_mode));
    Hi=H(omega_i);
    Hi_hat=Hi(2:4,2:4);
    Ei_hat=Hi(2:4,1);
    Ci_hat=[1; -Hi_hat\Ei_hat]
    
    C_hat(:,i_mode)=Ci_hat;
end

% Mode shapes computation

x = linspace(0,L,1000);
dx = x(2);

for i_mode=1:length(i_nat)
    omega_i=omega(i_nat(i_mode));
        gamma_i = sqrt(sqrt(m*omega_i^2/(E*J)));

    phi(i_mode,:)= C_hat(1,i_mode)*cos(gamma_i*x)  + C_hat(2,i_mode)*sin(gamma_i*x) + C_hat(3,i_mode)*cosh(gamma_i*x)  + C_hat(4,i_mode)*sinh(gamma_i*x);
end


% Plot of the modes

figure, grid on, box on, hold on
red = [1, 0, 0];
blue = [0, 0, 1];
colors_p = [linspace(red(1),blue(1),size(phi,1))', linspace(red(2),blue(2),size(phi,1))', linspace(red(3),blue(3),size(phi,1))'];
plot([0 L],[0 0],'--k','HandleVisibility','off')

for ii = 1:size(phi,1)
    plot(x,phi(ii,:),'LineWidth',2,'Color',[colors_p(ii,:)])
    %plot([L-0.5 L+0.5 L+0.5 L-0.5 L-0.5],[-2.6 -2.6 -3.4 -3.4 -2.6],'LineWidth',2,'Color',[colors_p(ii,:)],'HandleVisibility','off')
end
%ylim([-5 2])
legend({'Mode 1','Mode 2','Mode 3','Mode 4'},'Location','NorthOutside')

%%  FRFs

n = size(phi);
x_j = x(167);
x_k = x(1000);
Omega = 2*pi*f;
psi = 0.01;
s= 1i*Omega;

omega_1 = omega(i_nat(1));
gamma_1 = sqrt(sqrt(m*omega_1^2/(E*J)));

phi_funz_1 = @(xv) C_hat(1,1)*cos(gamma_1*xv)  + C_hat(2,1)*sin(gamma_1*xv) + C_hat(3,1)*cosh(gamma_1*xv)  + C_hat(4,1)*sinh(gamma_1*xv);
phi_vals_1 = phi_funz_1(x);


omega_2 = omega(i_nat(2));
gamma_2 = sqrt(sqrt(m*omega_2^2/(E*J)));

phi_funz_2 = @(xv) C_hat(1,2)*cos(gamma_2*xv)  + C_hat(2,2)*sin(gamma_2*xv) + C_hat(3,2)*cosh(gamma_2*xv)  + C_hat(4,2)*sinh(gamma_2*xv);
phi_vals_2 = phi_funz_2(x);



omega_3 = omega(i_nat(3));
gamma_3 = sqrt(sqrt(m*omega_3^2/(E*J)));

phi_funz_3 = @(xv) C_hat(1,3)*cos(gamma_3*xv)  + C_hat(2,3)*sin(gamma_3*xv) + C_hat(3,3)*cosh(gamma_3*xv)  + C_hat(4,3)*sinh(gamma_3*xv);
phi_vals_3 = phi_funz_3(x);

omega_4 = omega(i_nat(4));
gamma_4 = sqrt(sqrt(m*omega_4^2/(E*J)));

phi_funz_4 = @(xv) C_hat(1,4)*cos(gamma_4*xv)  + C_hat(2,4)*sin(gamma_4*xv) + C_hat(3,4)*cosh(gamma_4*xv)  + C_hat(4,4)*sinh(gamma_4*xv);
phi_vals_4 = phi_funz_4(x);


% normaliz1=max(abs([phi_vals_1]));
% normaliz2=max(abs([phi_vals_2]));
% normaliz3=max(abs([phi_vals_3]));
% normaliz4=max(abs([phi_vals_4]));
% 
% phi_vals_1=phi_vals_1./normaliz1;
% phi_vals_2=phi_vals_2./normaliz2;
% phi_vals_3=phi_vals_3./normaliz3;
% phi_vals_4=phi_vals_4./normaliz4;

integrando_1 = m*phi_vals_1.^2;
mm_1= trapz(x,integrando_1);

integrando_2 = m*phi_vals_2.^2;
mm_2= trapz(x,integrando_2);

integrando_3 = m*phi_vals_3.^2;
mm_3= trapz(x,integrando_3);

integrando_4 = m*phi_vals_4.^2;
mm_4= trapz(x,integrando_4);

mm = [mm_1 mm_2 mm_3 mm_4];



G_jk_1 = @(Omega) (phi(1,167)*phi(1,1000)/mm(1))./(-Omega.^2+1i*2*psi*omega_1*Omega+omega_1.^2);
G_jk_v_1 = G_jk_1(Omega);

G_jk_2 = @(Omega) (phi(2,167)*phi(2,1000)/mm(2))./(-Omega.^2+1i*2*psi*omega_2*Omega+omega_2.^2);
G_jk_v_2 = G_jk_2(Omega);

G_jk_3 = @(Omega) (phi(3,167)*phi(3,1000)/mm(3))./(-Omega.^2+1i*2*psi*omega_3*Omega+omega_3.^2);
G_jk_v_3 = G_jk_3(Omega);

G_jk_4 = @(Omega) (phi(4,167)*phi(4,1000)/mm(4))./(-Omega.^2+1i*2*psi*omega_4*Omega+omega_4.^2);
G_jk_v_4 = G_jk_4(Omega);

G_jk_a = G_jk_v_1 + G_jk_v_2 + G_jk_v_3 + G_jk_v_4;

magnitudea = abs(G_jk_a);
phasea = angle(G_jk_a)*(180/pi);

figure
subplot(2,1,1)
semilogy(f,magnitudea,'LineWidth',3)
title('FRF')
xlabel('Frequency [Hz]')
ylabel('Magnitude [m/N]')
ylim([1e-6 1e-1])
grid on
legend('x_j = 0.2 x_k = 1.2')
subplot(2,1,2)
plot(f,phasea,'LineWidth',3)
xlabel ('Frequency [Hz]')
ylabel('Phase [°]')
grid on

%% Other combinations

n = size(phi);
x_j = x(330);
x_k = x(1000);
Omega = 2*pi*f;
psi = 0.01;
s= 1i*Omega;

omega_1 = omega(i_nat(1));
gamma_1 = sqrt(sqrt(m*omega_1^2/(E*J)));

phi_funz_1 = @(xv) C_hat(1,1)*cos(gamma_1*xv)  + C_hat(2,1)*sin(gamma_1*xv) + C_hat(3,1)*cosh(gamma_1*xv)  + C_hat(4,1)*sinh(gamma_1*xv);
phi_vals_1 = phi_funz_1(x);

integrando_1 = m*phi_vals_1.^2;
mm_1= trapz(x,integrando_1);

omega_2 = omega(i_nat(2));
gamma_2 = sqrt(sqrt(m*omega_2^2/(E*J)));

phi_funz_2 = @(xv) C_hat(1,2)*cos(gamma_2*xv)  + C_hat(2,2)*sin(gamma_2*xv) + C_hat(3,2)*cosh(gamma_2*xv)  + C_hat(4,2)*sinh(gamma_2*xv);
phi_vals_2 = phi_funz_2(x);

integrando_2 = m*phi_vals_2.^2;
mm_2= trapz(x,integrando_2);

omega_3 = omega(i_nat(3));
gamma_3 = sqrt(sqrt(m*omega_3^2/(E*J)));

phi_funz_3 = @(xv) C_hat(1,3)*cos(gamma_3*xv)  + C_hat(2,3)*sin(gamma_3*xv) + C_hat(3,3)*cosh(gamma_3*xv)  + C_hat(4,3)*sinh(gamma_3*xv);
phi_vals_3 = phi_funz_3(x);

integrando_3 = m*phi_vals_3.^2;
mm_3= trapz(x,integrando_3);

omega_4 = omega(i_nat(4));
gamma_4 = sqrt(sqrt(m*omega_4^2/(E*J)));

phi_funz_4 = @(xv) C_hat(1,4)*cos(gamma_4*xv)  + C_hat(2,4)*sin(gamma_4*xv) + C_hat(3,4)*cosh(gamma_4*xv)  + C_hat(4,4)*sinh(gamma_4*xv);
phi_vals_4 = phi_funz_4(x);

integrando_4 = m*phi_vals_4.^2;
mm_4= trapz(x,integrando_4);

mm = [mm_1 mm_2 mm_3 mm_4];

G_jk_1 = @(Omega) (phi(1,330)*phi(1,1000)/mm(1))./(-Omega.^2+1i*2*psi*omega_1*Omega+omega_1.^2);
G_jk_v_1 = G_jk_1(Omega);

G_jk_2 = @(Omega) (phi(2,330)*phi(2,1000)/mm(2))./(-Omega.^2+1i*2*psi*omega_2*Omega+omega_2.^2);
G_jk_v_2 = G_jk_2(Omega);

G_jk_3 = @(Omega) (phi(3,330)*phi(3,1000)/mm(3))./(-Omega.^2+1i*2*psi*omega_3*Omega+omega_3.^2);
G_jk_v_3 = G_jk_3(Omega);

G_jk_4 = @(Omega) (phi(4,330)*phi(4,1000)/mm(4))./(-Omega.^2+1i*2*psi*omega_4*Omega+omega_4.^2);
G_jk_v_4 = G_jk_4(Omega);

G_jk_b = G_jk_v_1 + G_jk_v_2 + G_jk_v_3 + G_jk_v_4;

magnitudeb = abs(G_jk_b);
phaseb = angle(G_jk_b)*(180/pi);

figure
subplot(2,1,1)
semilogy(f,magnitudeb,'LineWidth',3)
title('FRF')
xlabel('Frequency [Hz]')
ylabel('Magnitude [m/N]')
ylim([1e-8 1e-1])
legend('x_j = 0.4 x_k = 1.2')
grid on
subplot(2,1,2)
plot(f,phaseb,'LineWidth',3)
xlabel ('Frequency [Hz]')
ylabel('Phase [°]')
grid on

%% FRFs numerically computed

f_min_1 = 4;
f_max_1 = 5;
i_min1=19952;
i_max1=24952;
freq1 = linspace(f_min_1,f_max_1,100);
G_jk_EXP = [G_jk_a; G_jk_b];
% for i=1:length(f)
%     if f(i)<f_max_1
%         i=i+1;
%     else
%     end
%         i_max1=i;
% end

params1 = [omega_1 psi 0.0565 0];
G_a_EXP_1 = interp1(f, G_jk_EXP(1,:), freq1, 'spline');
%G_a_EXP_1 = G_jk_EXP(f_min_1:f_max_1, 1);
%G_a_NUM_1_fun = @(params1,f)...
      %params1(3)./(-(f.^2)+1i*2*params1(2)*params1(1)*f+params1(1)^2) + params1(4); 
% params1(3) = A  params1(2) = psi  params1(1) = natural frequency  params1(4)= Rh (residuals of higher modes)
% G_a_NUM_1 = G_a_NUM_1_fun(params1,freq1);
% diff = G_a_EXP_1 - G_a_NUM_1;
% eps = sum(sum( real(diff).^2 + imag(diff).^2 ));
% [om1, psi1, A1, Rh1] = lsqnonlin(eps,params1);
% 
% fun = @(params) cost_function_FRF(params, Omega, G_exp)

function err = cost_function_FRF(params, freq, G_exp)

    om1 = params(1);
    psi1= params(2);
    A1= params(3);
    Rh1=params(4);
    
    % FRF del singolo modo
    G_model = A1./ (-(2*pi*freq).^2 + 1i*2*psi1*om1*2*pi*freq + (om1)^2)+Rh1;
    diff = G_exp - G_model;
    % Errore: parte reale e immaginaria
    %err = [real(G_model - G_exp); imag(G_model - G_exp)];
    err =sum(real(diff).^2 + imag(diff).^2);
end
err = @(params) cost_function_FRF(params, freq1, G_a_EXP_1);
% Limiti inferiori e superiori
lb = [zeros(1,4)];
ub = [Inf(1,4)];

% Ottimizzazione
opts = optimoptions('lsqnonlin','Display','iter','MaxFunctionEvaluations',5000);

x_opta_1 = lsqnonlin(err, params1, lb, ub, opts);
G_a_NUM_1 = x_opta_1(3)./ (-(2*pi*freq1).^2 + 1i*2*x_opta_1(2)*x_opta_1(1)*2*pi*freq1 + x_opta_1(1)^2)+x_opta_1(4);

magnitudea_NUM1 = abs(G_a_NUM_1);
phasea_NUM1 = angle(G_a_NUM_1)*(180/pi);

figure
subplot(2,1,1)
semilogy(f,magnitudea,'LineWidth',3)
hold on 
semilogy(freq1,magnitudea_NUM1,'or')
title('FRF')
xlabel('Frequency [Hz]')
ylabel('Magnitude [m/N]')
ylim([1e-6 1e-1])
xlim([3 6])
grid on
legend('x_j = 0.2 x_k = 1.2')

subplot(2,1,2)
plot(f,phasea,'LineWidth',3)
hold on 
plot(freq1,phasea_NUM1,'or')
xlim([3 6])
xlabel ('Frequency [Hz]')
ylabel('Phase [°]')
grid on

% Other FRF (b)

paramsb1 = [omega_1 psi 0.65 0];
G_b_EXP_1 = interp1(f, G_jk_EXP(2,:), freq1, 'spline');
err = @(params) cost_function_FRF(paramsb1,freq1,G_b_EXP_1);
x_optb_1 = lsqnonlin(err, paramsb1, lb, ub, opts);
G_b_NUM_1 = x_optb_1(3)./ (-(2*pi*freq1).^2 + 1i*2*x_optb_1(2)*x_optb_1(1)*2*pi*freq1 + x_optb_1(1)^2)+x_optb_1(4);

magnitudeb_NUM1 = abs(G_b_NUM_1);
phaseb_NUM1 = angle(G_b_NUM_1)*(180/pi);

figure
subplot(2,1,1)
semilogy(f,magnitudeb,'LineWidth',3)
hold on 
semilogy(freq1,magnitudeb_NUM1,'or')
title('FRF')
xlabel('Frequency [Hz]')
ylabel('Magnitude [m/N]')
ylim([1e-6 1e-1])
xlim([3 6])
grid on
legend('x_j = 0.4 x_k = 1.2')

subplot(2,1,2)
plot(f,phaseb,'LineWidth',3)
hold on 
plot(freq1,phaseb_NUM1,'or')
xlim([3 6])
xlabel ('Frequency [Hz]')
ylabel('Phase [°]')
grid on

%% Computation for the other modes

% Second mode 
freq2 = linspace(26,30,200);
paramsa2 = [omega_2 psi 0.85 0 0];
G_a_EXP_2 = interp1(f, G_jk_EXP(1,:), freq2, 'spline');

function err = cost_function_FRF_seismic(params, freq, G_exp)

    om1 = params(1);
    psi1= params(2);
    A1= params(3);
    Rh1=params(4);
    Rk1=params(5);
    
    % FRF del singolo modo
    G_model = A1./ (-(2*pi*freq).^2 + 1i*2*psi1*om1*2*pi*freq + (om1)^2)+Rh1 + Rk1./((2*pi*freq).^2);
    diff = G_exp - G_model;
    % Errore: parte reale e immaginaria
    %err = [real(G_model - G_exp); imag(G_model - G_exp)];
    err =sum(real(diff).^2 + imag(diff).^2);
end

err = @(params) cost_function_FRF_seismic(paramsa2,freq2,G_a_EXP_2);
x_opta_2 = lsqnonlin(err, paramsa2, lb, ub, opts);
G_a_NUM_2 = x_opta_2(3)./ (-(2*pi*freq2).^2 + 1i*2*x_opta_2(2)*x_opta_2(1)*2*pi*freq2 + x_opta_2(1)^2)+x_opta_2(4)+x_opta_2(5)./((2*pi*freq2).^2);

magnitudea_NUM2 = abs(G_a_NUM_2);
phasea_NUM2 = angle(G_a_NUM_2)*(180/pi);

paramsb2 = [omega_2 psi 2.3 0 0];
G_b_EXP_2 = interp1(f, G_jk_EXP(2,:), freq2, 'spline');
err = @(params) cost_function_FRF_seismic(paramsb2,freq2,G_b_EXP_2);
x_optb_2 = lsqnonlin(err, paramsb2, lb, ub, opts);
G_b_NUM_2 = x_optb_2(3)./ (-(2*pi*freq2).^2 + 1i*2*x_optb_2(2)*x_optb_2(1)*2*pi*freq2 + x_optb_2(1)^2)+x_optb_2(4)+x_optb_2(5)./((2*pi*freq2).^2);

magnitudeb_NUM2 = abs(G_b_NUM_2);
phaseb_NUM2 = angle(G_b_NUM_2)*(180/pi);

% Third mode

freq3 = linspace(76,82,200);
paramsa3 = [omega_3 psi 1.87 0 0];
G_a_EXP_3 = interp1(f, G_jk_EXP(1,:), freq3, 'spline');

err = @(params) cost_function_FRF_seismic(paramsa3,freq3,G_a_EXP_3);
x_opta_3 = lsqnonlin(err, paramsa3, lb, ub, opts);
G_a_NUM_3 = x_opta_3(3)./ (-(2*pi*freq3).^2 + 1i*2*x_opta_3(2)*x_opta_3(1)*2*pi*freq3 + x_opta_3(1)^2)+x_opta_3(4)+x_opta_3(5)./((2*pi*freq3).^2);

magnitudea_NUM3 = abs(G_a_NUM_3);
phasea_NUM3 = angle(G_a_NUM_3)*(180/pi);

paramsb3 = [omega_3 psi 2.8 0 0];
G_b_EXP_3 = interp1(f, G_jk_EXP(2,:), freq3, 'spline');
err = @(params) cost_function_FRF_seismic(paramsb3,freq3,G_b_EXP_3);
x_optb_3 = lsqnonlin(err, paramsb3, lb, ub, opts);
G_b_NUM_3 = x_optb_3(3)./ (-(2*pi*freq3).^2 + 1i*2*x_optb_3(2)*x_optb_3(1)*2*pi*freq3 + x_optb_3(1)^2)+x_optb_3(4)+x_optb_3(5)./((2*pi*freq3).^2);

magnitudeb_NUM3 = abs(G_b_NUM_3);
phaseb_NUM3 = angle(G_b_NUM_3)*(180/pi);

% Fourth mode

freq4 = linspace(150,160,100);
paramsa4 = [omega_4 psi 2.72 0 ];
G_a_EXP_4 = interp1(f, G_jk_EXP(1,:), freq4, 'spline');

function err = cost_function_FRF_static(params, freq, G_exp)

    om1 = params(1);
    psi1= params(2);
    A1= params(3);
    Rk1=params(4);
    
    % FRF del singolo modo
    G_model = A1./ (-(2*pi*freq).^2 + 1i*2*psi1*om1*2*pi*freq + (om1)^2) + Rk1./((2*pi*freq).^2);
    diff = G_exp - G_model;
    % Errore: parte reale e immaginaria
    %err = [real(G_model - G_exp); imag(G_model - G_exp)];
    err =sum(real(diff).^2 + imag(diff).^2);
end

err = @(params) cost_function_FRF_static(paramsa4,freq4,G_a_EXP_4);
x_opta_4 = lsqnonlin(err, paramsa4, lb, ub, opts);
G_a_NUM_4 = x_opta_4(3)./ (-(2*pi*freq4).^2 + 1i*2*x_opta_4(2)*x_opta_4(1)*2*pi*freq4 + x_opta_4(1)^2)+x_opta_4(4)./((2*pi*freq4).^2);

magnitudea_NUM4 = abs(G_a_NUM_4);
phasea_NUM4 = angle(G_a_NUM_4)*(180/pi);

paramsb4 = [omega_4 psi 0.85 0];
G_b_EXP_4 = interp1(f, G_jk_EXP(2,:), freq4, 'spline');
err = @(params) cost_function_FRF_static(paramsb4,freq4,G_b_EXP_4);
x_optb_4 = lsqnonlin(err, paramsb4, lb, ub, opts);
G_b_NUM_4 = x_optb_4(3)./ (-(2*pi*freq4).^2 + 1i*2*x_optb_4(2)*x_optb_4(1)*2*pi*freq4 + x_optb_4(1)^2)+x_optb_4(4)./((2*pi*freq4).^2);

magnitudeb_NUM4 = abs(G_b_NUM_4);
phaseb_NUM4 = angle(G_b_NUM_4)*(180/pi);

figure
subplot(2,2,1)
semilogy(f,magnitudea,'LineWidth',3)
hold on 
semilogy(freq1,magnitudea_NUM1,'or')
hold on
semilogy(freq2,magnitudea_NUM2,'or')
hold on
semilogy(freq3,magnitudea_NUM3,'or')
hold on
semilogy(freq4,magnitudea_NUM4,'or')
title('FRF')
xlabel('Frequency [Hz]')
ylabel('Magnitude [m/N]')
grid on
legend('$G^{EXP}_{jk} x_j = 0.2 x_k = 1.2$','$G^{NUM}_{jk}$','Interpreter','latex')

subplot(2,2,2)
semilogy(f,magnitudeb,'LineWidth',3)
hold on 
semilogy(freq1,magnitudeb_NUM1,'or')
hold on
semilogy(freq2,magnitudeb_NUM2,'or')
hold on
semilogy(freq3,magnitudeb_NUM3,'or')
hold on
semilogy(freq4,magnitudeb_NUM4,'or')
title('FRF')
xlabel('Frequency [Hz]')
ylabel('Magnitude [m/N]')
grid on
legend('$G^{EXP}_{jk}  x_j = 0.4 x_k = 1.2$','$G^{NUM}_{jk}$','Interpreter','latex')


subplot(2,2,3)
plot(f,phasea,'LineWidth',3)
hold on 
plot(freq1,phasea_NUM1,'or')
hold on
plot(freq2,(180+phasea_NUM2),'or')
hold on
plot(freq3,phasea_NUM3,'or')
hold on
plot(freq4,(180+phasea_NUM4),'or')
xlabel ('Frequency [Hz]')
ylabel('Phase [°]')
grid on

subplot(2,2,4)
plot(f,phaseb,'LineWidth',3)
hold on 
plot(freq1,phaseb_NUM1,'or')
hold on
plot(freq2,(180+phaseb_NUM2),'or')
hold on
plot(freq3,phaseb_NUM3,'or')
hold on
plot(freq4,(180+phaseb_NUM4),'or')
xlabel ('Frequency [Hz]')
ylabel('Phase [°]')
grid on

%% COMPARISON

ck_1 = omega_1*psi*2*mm_1;

figure, grid on, box on, hold on

plot([0 L],[0 0],'--k','HandleVisibility','off')

plot(x,phi(1,:),'LineWidth',2)
hold on
plot(x(167),-imag(max(G_a_NUM_1)*(omega_1*ck_1)/phi(1,1000)),'ob','LineWidth',2)
hold on 
plot(x(330),-imag(max(G_b_NUM_1)*(omega_1*ck_1)/phi(1,1000)),'ob','LineWidth',2)