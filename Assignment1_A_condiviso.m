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
% G_jk_v = 0;
s= 1i*Omega;
% mm = zeros(n);

% for i_mode=1:4
%     omega_i=omega(i_nat(i_mode));
%         gamma_i = sqrt(sqrt(m*omega_i^2/(E*J)));
% 
%     phi_funz = @(xv) C_hat(1,i_mode)*cos(gamma_i*xv)  + C_hat(2,i_mode)*sin(gamma_i*xv) + C_hat(3,i_mode)*cosh(gamma_i*xv)  + C_hat(4,i_mode)*sinh(gamma_i*xv);
%     phi_vals(i,:) = phi_funz(x);
%     integrando = m*phi_vals.^2;
%     mm(i)= trapz(x,integrando);
% end

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

% for i=1:4
%     G_jk = @(Omega) G_jk_v+(phi(i,167)*phi(i,1000)/mm(i))./(-Omega.^2+1i*2*psi*omega(i)*Omega+omega(i).^2);
%     G_jk_v = G_jk(s);
% end

G_jk_1 = @(Omega) (phi(1,167)*phi(1,1000)/mm(1))./(-Omega.^2+1i*2*psi*omega(1)*Omega+omega(1).^2);
G_jk_v_1 = G_jk_1(Omega);

G_jk_2 = @(Omega) (phi(2,167)*phi(2,1000)/mm(2))./(-Omega.^2+1i*2*psi*omega(2)*Omega+omega(2).^2);
G_jk_v_2 = G_jk_2(Omega);

G_jk_3 = @(Omega) (phi(3,167)*phi(3,1000)/mm(3))./(-Omega.^2+1i*2*psi*omega(3)*Omega+omega(3).^2);
G_jk_v_3 = G_jk_3(Omega);

G_jk_4 = @(Omega) (phi(4,167)*phi(4,1000)/mm(4))./(-Omega.^2+1i*2*psi*omega(4)*Omega+omega(4).^2);
G_jk_v_4 = G_jk_4(Omega);

G_jk_v = G_jk_v_1 + G_jk_v_2 + G_jk_v_3 + G_jk_v_4;


magnitude = abs(G_jk_v);
phase = angle(G_jk_v)*(180/pi);

figure
semilogx(Omega,magnitude)
grid on