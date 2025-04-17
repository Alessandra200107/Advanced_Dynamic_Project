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
mm = zeros(n);

for i=1:n
    m_funz = @(xv) m*phi(i,xv).^2;
    mm(i)= trapz(x,m_funz);
end

for i=1:n
    G_jk = @(Omega) G_jk+(phi(i,167)*phi(i,1000)/mm(i))/(-Omega.^2+1i*2*psi*omega(i)*Omega+omega(i).^2);
end