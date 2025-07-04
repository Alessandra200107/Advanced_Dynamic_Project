%% ASSIGNMENT 2 - FEM

clear all
close all
clc

% load input file and assemble structure
[file_i,xy,nnod,sizee,idb,ndof,incidenze,l,gamma,m,EA,EJ,posiz,nbeam,pr]=loadstructure;

% draw structure
dis_stru(posiz,l,gamma,xy,pr,idb,ndof);

%%

% assemble mass and stiffness matrix
[M,K] = assem(incidenze,l,m,EA,EJ,gamma,idb);

M_FF = M(1:ndof,1:ndof);
M_FC = M(1:ndof,ndof+1:end);
M_CF = M(ndof+1:end,1:ndof);
M_CC = M(ndof+1:end,ndof+1:end);

K_FF = K(1:ndof,1:ndof);
K_FC = K(1:ndof,ndof+1:end);
K_CF = K(ndof+1:end,1:ndof);
K_CC = K(ndof+1:end,ndof+1:end);

[modes,omega2] = eig(M_FF\K_FF);
omega = sqrt(diag(omega2));
[omega,i_omega] = sort(omega);
freq0 = omega/2/pi;
modes = modes(:,i_omega);

scale_factor = 1.5;

figure()
diseg2(modes(:,1),scale_factor,incidenze,l,gamma,posiz,idb,xy);
title('Mode 1: freq[Hz] = 6.4821')
xlabel('x [m]')
ylabel('y [m]')

figure()
diseg2(modes(:,2),scale_factor,incidenze,l,gamma,posiz,idb,xy);
title('Mode 2: freq[Hz] = 6.7823')
xlabel('x [m]')
ylabel('y [m]')

figure()
diseg2(modes(:,3),scale_factor,incidenze,l,gamma,posiz,idb,xy);
title('Mode 3: freq[Hz] = 14.5331')
xlabel('x [m]')
ylabel('y [m]')

%% DAMPING

ab = [0.1 2e-4];
C = ab(1)*M + ab(2)*K;

C_FF = C(1:ndof,1:ndof);
C_FC = C(1:ndof,ndof+1:end);
C_CF = C(ndof+1:end,1:ndof);
C_CC = C(ndof+1:end,ndof+1:end);


%% FRF between the force applied in A and the vertical displacement in A

F0 = zeros(ndof,1);
F0((idb(21,2))) = 1;
om = (0:0.01:20)*2*pi;
for ii = 1:length(om)
    A_force = -om(ii)^2*M_FF + 1i*om(ii)*C_FF + K_FF;
    Xf(:,ii) = A_force\F0;
    Xfp(:,ii) = 1i*om(ii)*Xf(:,ii);
    Xfpp(:,ii) = -om(ii)^2*Xf(:,ii); % to obtain the acceleration spectrum
end

FRF_dispA = Xf(idb(21,2),:);

figure()
subplot(2,1,1)
semilogy(om/(2*pi),abs(FRF_dispA),'LineWidth',3)
title('FRF with respect to the vertical displacement of point A')
xlabel('Frequency [Hz]')
ylabel('Magnitude [m/N]')
grid on
subplot(2,1,2)
plot(om/(2*pi),(angle(FRF_dispA)*180/pi),'LineWidth',3)
xlabel ('Frequency [Hz]')
ylabel('Phase [°]')
grid on

FRF_dispB = Xf(idb(7,2),:);

figure()
subplot(2,1,1)
semilogy(om/(2*pi),abs(FRF_dispB),'LineWidth',3)
title('FRF with respect to the vertical displacement of point B')
xlabel('Frequency [Hz]')
ylabel('Magnitude [m/N]')
grid on
subplot(2,1,2)
plot(om/(2*pi),(angle(FRF_dispB)*180/pi),'LineWidth',3)
xlabel ('Frequency [Hz]')
ylabel('Phase [°]')
grid on

%% MODAL APPROACH

ii = 1:2;
Phi = modes(:,ii);
Mmod = Phi'*M_FF*Phi;
Kmod = Phi'*K_FF*Phi;
Cmod = Phi'*C_FF*Phi;
Fmod = Phi'*F0;
for ii = 1:length(om)
    xx_mod(:,ii) = (-om(ii)^2*Mmod + 1i*om(ii)*Cmod +Kmod)\Fmod;
end
xx_m = Phi*xx_mod;
FRF_modB = xx_m(idb(7,2),:);

figure()
subplot(2,1,1)
semilogy(om/(2*pi),abs(FRF_dispB),'LineWidth',3)
title('FRF with respect to the vertical displacement of point B')
xlabel('Frequency [Hz]')
ylabel('Magnitude [m/N]')
grid on
hold on
semilogy(om/(2*pi),abs(FRF_modB),'LineWidth',3)
legend('FEM','Modal')
subplot(2,1,2)
plot(om/(2*pi),(angle(FRF_dispB)*180/pi),'LineWidth',3)
xlabel ('Frequency [Hz]')
ylabel('Phase [°]')
grid on
hold on
plot(om/(2*pi),(angle(FRF_modB)*180/pi),'LineWidth',3)

%% STATIC RESPONSE DUE TO THE WEIGHT

g = 9.81; 
acc = zeros(ndof,1);
acc_c = zeros(nnod*3-ndof,1);
vert = idb(:,2);

% for ii = 1:length(vert)
%     if idb(ii,2)<=170
%         vert_dof(ii) = idb(ii,2);
%     end
% end

vert_c = [vert(1); vert(41)];
% vert(1) = [];
% vert(40) = [];
acc(vert) = -g;
acc_c(vert_c-ndof) = -g;
Fg = M*acc;
Fg_FF = Fg(1:ndof);
Xg = K_FF\Fg_FF;
% Fg_FC = M_CF*acc;
% Xg_FC = K_CF\Fg_FC;
% Xg = [Xg_FF; Xg_FC];

figure()
diseg2(Xg,150,incidenze,l,gamma,posiz,idb,xy);
title('Static deflecton due to the weight')

max_def = max(abs((Xg))); % maximum vertical displacement in m
fprintf('The maximum vertical displacement is %.6f mm\n', max_def*1000);