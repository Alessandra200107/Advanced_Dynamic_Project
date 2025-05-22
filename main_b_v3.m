clear
clc
close all
set(0,'defaultfigurecolor',[1 1 1])
set(0, 'DefaultFigureWindowStyle', 'docked')

%% data
load("Data.mat")
figure()
subplot(3, 1, 1)

plot(freq, log(abs(frf)))
ylim([-10, 0])
legend

subplot(3, 1, 2)
plot(freq, angle(frf))
subplot(3, 1, 3)
plot(freq, cohe)

N = size(frf, 2);
dw = (freq(2)-freq(1))*2*pi;

xj = 15:15:180;
xk = zeros(N, 1);

%% Modal parameters identification, simplified procedure
%Identificazione campi di interesse frequenze proprie

ei = 1; %Eigenfrequency of interest


if ei == 1
    range = dw:dw:875*2*pi;
    sensors =  [1:3, 5:8, 9, 11, 12];
 
elseif ei == 2
    range = 1000*2*pi:dw:2000*2*pi;
    % sensors =  1:12;
     sensors = [1:3, 5:6, 8:10, 12];
elseif ei == 4
        range = 3800*2*pi:dw:5000*2*pi;
        sensors = [1:6, 8:11];
end

frf = [frf(:, sensors)];
xj = [xj(:, sensors)];
N = size(frf, 2);
M = length(range);
omega = 2*pi*freq;
[omega0, csi_exp, Ajk, Gr_exp] = modal_parameters_simplified_4(frf, range, omega);

Re_Rh = zeros(1, N);
Im_Rh = zeros(1, N);
Re_Rl = zeros(1, N);
Im_Rl = zeros(1, N);

%% Blocco smorzamento

syms x omega_s
% x rappresenta le variabili i cui errori vanno minimizzati:  wi, csii, Ai, Rl,
% Rh
Gr_num_og = @(x, omega_s, ii) x(ii)./(-omega_s.^2 + 2*1i*csi_exp(1)*omega0(1).*omega_s +omega0(1).^2) + (x(N+ii)  +1i*x(2*N+ii)) + (x(3*N+ii) + 1i*x(4*N+ii))./(omega_s.^2);

x_index = 1:N;

Gr_num = @(x) Gr_num_og(x, range, x_index);


 % 
 % Gr_num = @(x) ones(1, M).*x(2.*ones(N,1)+x_index')./(-range.^2 + 2*1i*x(2)*x(1).*range +x(1).^2) ...
 %     + ones(1, M).*x((N+2).*ones(N,1)+x_index') ...
 %     + ones(1, M).*x((2*N+2).*ones(N,1)+ x_index')./(range.^2);

diff = @(x) real(Gr_exp-Gr_num(x)).^2 + imag(Gr_exp-Gr_num(x)).^2;

% err = @(x) sum(sum(  (Gr_exp - Gr_num(x, range, x_index')).*conj(Gr_exp-Gr_num(x, range, x_index'))));
err = @(x) sum(sum(  diff(x)));

x0 = [Ajk;
    Re_Rh';
    Im_Rh';
    Re_Rl';
    Im_Rl';];

figure()
plot(range, abs(Gr_num(x0)))
title('x0')

options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt', 'StepTolerance', 10e-20);

soluzione = lsqnonlin(@(x) err(x), x0, [], []);
% Gr_num(x0, range, x_index');
err(soluzione)

Gr_sol = Gr_num(soluzione);

for ii = 1:size(Gr_sol, 1)
    % subplot(size(Gr_sol, 1), 1, ii)
    figure()
    plot(range, abs(Gr_exp(ii, :)), 'r', range, abs(Gr_sol(ii, :)), 'ob')
    title(['xj = ', num2str(xj(ii)), ' xk = ', num2str(xk(ii))])
end

% for ii = 1:size(Gr_sol, 3)
%     title(['xj = ', num2str(xj(ii)), ' xk = ', num2str(xk(ii))])
%     subplot(3, 1, ii)
%     plot(range, abs(Gr_exp(ii, :)), 'r', range, abs(Gr_sol(ii, :)), 'ob')
% end

% %% 3 - deformed shape plot
% mod_sperimentale = soluzione(1:N);
% mod_sperimentale = mod_sperimentale./(max(abs(mod_sperimentale)));
% 
% theta = deg2rad([xj -xj]);
% R_base = 10;
% rho = abs([R_base+mod_sperimentale; R_base+mod_sperimentale]);
% 
% % theta_vect 
% 
% figure()
% 
% 
% 
% polarplot(theta, rho, '*')
% hold on
% title(['Deformed shape n°: ', num2str(ei), ' f0 = ', num2str(omega0(1)/(2*pi)), ' damping ratio: ', num2str(csi_exp(1))])
% 


%% 3 - deformed shape plot with interpolation 
mod_sperimentale = soluzione(1:N);
mod_sperimentale = mod_sperimentale ./ max(abs(mod_sperimentale)); % Normalizzazione

% Coordinate polari
theta = deg2rad([xj -xj]); % Angoli originali
R_base = 10; % Valore base del raggio
rho = abs([R_base + mod_sperimentale; R_base + mod_sperimentale]);

% Interpolazione spline
theta_interp = linspace(min(theta), max(theta), 1000); % Maggiore risoluzione angolare
rho_interp = interp1(theta, rho, theta_interp, 'spline'); % Interpolazione spline

% Circonferenza di riferimento
theta_circle = linspace(0, 2*pi, 1000); % Angoli per la circonferenza
rho_circle = R_base * ones(size(theta_circle)); % Raggio costante

% Plot
figure()
polarplot(theta, rho, '*', 'DisplayName', 'Punti originali') % Punti originali
hold on
polarplot(theta_interp, rho_interp, '-', 'DisplayName', 'Interpolazione spline','LineWidth',1.5) % Curva interpolata
polarplot(theta_circle, rho_circle, '--k', 'DisplayName', 'Circonferenza di riferimento') % Circonferenza di riferimento
title(['Deformed shape n°: ', num2str(ei), ' f0 = ', num2str(omega0(1)/(2*pi)), ...
    ' Hz, damping ratio: ', num2str(csi_exp(1))])
legend('Location', 'best')
hold off
