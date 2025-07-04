%% ASSIGNMENT 1 B - EXPERIMENTAL MODAL ANALYSIS OF A LIGHT RAIL WHEEL

clear all
close all
clc
load('Data.mat');

f = freq;
FRF = frf;
ch = cohe;

for i=1:12
    magnitude(:,i) = abs(FRF(:,i));
    phase(:,i) = angle(FRF(:,i))*(180/pi);
end

figure()
subplot(3,1,1)
for i=1:12
    semilogy(f,magnitude(:,i))
    hold on
    grid on
end
legend('1','2','3','4','5','6','7','8','9','10','11','12','Location','northeastoutside')
ylim([1e-4 2])
subplot(3,1,2)
for i=1:12
    plot(f,phase(:,i))
    hold on
    grid on
end
subplot(3,1,3)
for i=1:12
    plot(f,ch(:,i))
    hold on
    grid on
end
mag_c = cell(1,12);
omega_c = cell(1,12);

for i = 1:12
    [mag_c{i}, omega_c{i}] = findpeaks(magnitude(:,i), f, 'MinPeakProminence', 0.1);
end

mag = nan(2, 12);    
omega = nan(2, 12);

for i = 1:12
    omega(1:2, i) = omega_c{i}(1:2);
    mag(1:2, i) = mag_c{i}(1:2);
end

omega_1 = omega(1,1)*2*pi;
omega_2 = omega(2,1)*2*pi;

psi = 0.008;

%% Numerical computation

% First mode

f_min = 640;
f_max = 700;
freq1 = linspace(f_min,f_max,500);
for i = 1:12
    FRF_mod1(:,i) = interp1(f, FRF(:,i), freq1, 'spline');
end

Var = zeros(12,4);
FRF_num1 =zeros(12,length(freq1));
magnitude_num1 = zeros(12,length(freq1));
phase_num1 = zeros(12,length(freq1));

function err = cost_function1(params,freq,G_exp)
    om1 = params(1);
    psi1= params(2);
    A1= params(3);
    Rh1=params(4);
    
    % FRF del singolo modo
    G_model = A1./ (-(freq*2*pi).^2 + 1i*2*psi1*om1*freq*2*pi + (om1)^2)+Rh1;
    diff = G_exp - G_model;
    % Errore: parte reale e immaginaria
    %err = [real(G_model - G_exp); imag(G_model - G_exp)];
    err = sum(real(diff).^2 + imag(diff).^2);
end

% Per la stima di Ajk partiamo dalla formula mag(1,i)*2*psi*(omega_1^2)

for i = 1:12
    params1 = [omega_1 psi mag(1,i)*2.8129e+05 0];

    err = @(params) cost_function1(params1,freq1,FRF_mod1(:,i));
    % err_v(i) = err;
    lb = [zeros(1,4)];
    ub = [Inf(1,4)];

    % Ottimizzazione
    opts = optimoptions('lsqnonlin','Display','iter','MaxFunctionEvaluations',5000);

    x_opt_1 = lsqnonlin(err, params1, lb, ub, opts);
    Var(i,:) = x_opt_1;
    %FRF_num1_s = x_opt_1(3)./ (-(freq1*2*pi).^2 + 1i*2*x_opt_1(2)*x_opt_1(1)*freq1*2*pi + x_opt_1(1)^2)+x_opt_1(4);
    den = -(freq1*2*pi).^2 + 1i * 2 * Var(i,2) * Var(i,1) * freq1*2*pi + Var(i,1)^2;
    FRF_num1_s = Var(i,3) ./ den + Var(i,4);
    FRF_num1(i,:) = FRF_num1_s;

    magnitude_num1(i,:) = abs(FRF_num1_s);
    phase_num1(i,:) = angle(FRF_num1_s)*(180/pi);
end

figure
subplot(2,1,1)
semilogy(f,magnitude(:,1),'-b')
hold on 
semilogy(freq1,magnitude_num1(1,:),'--','Color','#0072BD','LineWidth',2)
hold on
semilogy(f,magnitude(:,9),'-r')
hold on
semilogy(freq1,magnitude_num1(9,:),'--m','LineWidth',2')
title('FRF')
legend('FRF exp sensor 1','FRF num sensor 1','FRF exp sensor 9','FRF num sensor 9')
xlabel('Frequency [Hz]')
ylabel('Magnitude [m/N]')
xlim([f_min-10 f_max+10])
grid on

subplot(2,1,2)
plot(f,phase(:,1),'-b')
hold on
plot(freq1,(phase_num1(1,:)),'--','Color','#0072BD','LineWidth',2)
hold on
plot(f,phase(:,9),'-r')
hold on
plot(freq1,(180+phase_num1(9,:)),'--','Color','m','LineWidth',2)
xlabel ('Frequency [Hz]')
ylabel('Phase [°]')
xlim([f_min-10 f_max+10])
grid on

% Second mode

f_min = 1600;
f_max = 1660;
freq2 = linspace(f_min,f_max,500);
for i = 1:12
    FRF_mod2(:,i) = interp1(f, FRF(:,i), freq2, 'spline');
end

Var2 = zeros(12,5);
FRF_num2 = zeros(12,length(freq2));
magnitude_num2 = zeros(12,length(freq2));
phase_num2 = zeros(12,length(freq2));

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
    %err = [real(diff); imag(diff)];
    err = sum(sum((real(diff).^2 + imag(diff).^2)));
end

for i = 1:12
    params2 = [omega_2 0.0075 mag(2,i)*omega_2^2*2*0.0075 0 0];

    err = @(params) cost_function_FRF_seismic(params2,freq2,FRF_mod2(:,i));
    % err_v(i) = err;
    lb = [0 0 0 -Inf -Inf];
    ub = [Inf Inf Inf Inf Inf];
    % Ottimizzazione
    opts = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt', 'StepTolerance', 10e-20);

    x_opt_2 = lsqnonlin( err, params2, [], []);
    Var2(i,:) = x_opt_2;
    den = -(freq2*2*pi).^2 + 1i * 2 * Var2(i,2) * Var2(i,1) * freq2*2*pi + Var2(i,1)^2;
    FRF_num2_s = Var2(i,3) ./ den + Var2(i,4) + Var2(i,5)./((2*pi*freq2).^2);
    FRF_num2(i,:) = FRF_num2_s;

    magnitude_num2(i,:) = abs(FRF_num2_s);
    phase_num2(i,:) = angle(FRF_num2_s)*(180/pi);
end

figure
subplot(2,1,1)
semilogy(f,magnitude(:,1),'-b')
hold on 
semilogy(freq2,magnitude_num2(1,:),'--','Color','#0072BD','LineWidth',2)
hold on
semilogy(f,magnitude(:,9),'-r')
hold on
semilogy(freq2,magnitude_num2(9,:),'--m','LineWidth',2')
title('FRF')
legend('FRF exp sensor 1','FRF num sensor 1','FRF exp sensor 9','FRF num sensor 9')
xlabel('Frequency [Hz]')
ylabel('Magnitude [m/N]')
xlim([f_min-10 f_max+10])
grid on

subplot(2,1,2)
plot(f,phase(:,1),'-b')
hold on
plot(freq2,(180+phase_num2(1,:)),'--','Color','#0072BD','LineWidth',2)
hold on
plot(f,phase(:,9),'-r')
hold on
plot(freq2,(180+phase_num2(9,:)),'--m','LineWidth',2)
xlabel ('Frequency [Hz]')
ylabel('Phase [°]')
xlim([f_min-10 f_max+10])
grid on

%% POLAR PLOTS

% Numero di punti di misura (es. 12 sensori)
N = 12;

% Angoli (in radianti) equispaziati
theta = linspace(0, pi, N+1);  % aggiungiamo un punto per chiudere la curva
theta(end) = [];  % togliamo l'ultimo punto (duplicato a 0/360°)

% Forma modale identificata (esempio fittizio)
% Qui metti i valori reali ottenuti dal tuo calcolo modale
idx = 2003;
mode_shape = [real(FRF(idx,1)) real(FRF(idx,2)) real(FRF(idx,3)) real(FRF(idx,4)) real(FRF(idx,5)) real(FRF(idx,6)) real(FRF(idx,7)) real(FRF(idx,8)) real(FRF(idx,9)) real(FRF(idx,10)) real(FRF(idx,11)) real(FRF(idx,12))];
mode_shape_num = [real(FRF_num1(1,228)) real(FRF_num1(2,228)) real(FRF_num1(3,228)) real(FRF_num1(4,228)) real(FRF_num1(5,228)) real(FRF_num1(6,228)) real(FRF_num1(7,228)) real(FRF_num1(8,228)) real(FRF_num1(9,228)) real(FRF_num1(10,228)) real(FRF_num1(11,228)) real(FRF_num1(12,228))];

% Forma simmetrica (opzionale)
mode_sym = -mode_shape;

theta_dense = linspace(0,(pi-15*pi/180),300);
mode_shape_interp = spline(theta,mode_shape,theta_dense);
mode_sym_interp = -mode_shape_interp;


% Raggio del cerchio indeformato
r0 = 1;
circle = r0 * ones(size(theta_dense));

% --- Plot ---
figure; 
pax = polaraxes;
pax.ThetaZeroLocation = 'bottom';
hold on

% Cerchio indeformato
polarplot(theta_dense, circle, 'k--', 'LineWidth', 1.2)
polarplot(-theta_dense, circle, 'k--', 'LineWidth', 1.2,'HandleVisibility','off')

% Forma modale (deformata aggiunta al cerchio)
polarplot(theta_dense, r0 + mode_shape_interp, 'r', 'LineWidth', 2)

% Forma simmetrica (opzionale)
polarplot(-theta_dense, r0 - mode_sym_interp, 'b', 'LineWidth', 2)

polarplot(theta, r0 + mode_shape, 'or')

legend('Undeformed','Axial mode shape (identified)', 'Axial mode shape (symmetry)', 'Location','best')
title('Modal Shape in Polar Coordinates f = 667.324 xi = 0.008')


%% Polar plot del secondo modo

% Numero di punti di misura (es. 12 sensori)
N = 12;

% Angoli (in radianti) equispaziati
theta = linspace(0, pi, N+1);  % aggiungiamo un punto per chiudere la curva
theta(end) = [];  % togliamo l'ultimo punto (duplicato a 0/360°)

% Forma modale identificata (esempio fittizio)
% Qui metti i valori reali ottenuti dal tuo calcolo modale
idx2 = 4877;
mode_shape = [real(FRF(idx2,1)) real(FRF(idx2,2)) real(FRF(idx2,3)) real(FRF(idx2,4)) real(FRF(idx2,5)) real(FRF(idx2,6)) real(FRF(idx2,7)) real(FRF(idx2,8)) real(FRF(idx2,9)) real(FRF(idx2,10)) real(FRF(idx2,11)) real(FRF(idx2,12))];

% Forma simmetrica (opzionale)
mode_sym = -mode_shape;

theta_dense = linspace(0,(pi-15*pi/180),300);
mode_shape_interp = spline(theta,mode_shape,theta_dense);
mode_sym_interp = -mode_shape_interp;

mode_shape_num = [real(FRF_num2(1,211)) real(FRF_num2(2,211)) real(FRF_num2(3,211)) real(FRF_num2(4,211)) real(FRF_num2(5,211)) real(FRF_num2(6,211)) real(FRF_num2(7,211)) real(FRF_num2(8,211)) real(FRF_num2(9,211)) real(FRF_num2(10,211)) real(FRF_num2(11,211)) real(FRF_num2(12,211))];

% Raggio del cerchio indeformato
r0 = 1;
circle = r0 * ones(size(theta_dense));

% --- Plot ---
figure; 
pax = polaraxes;
pax.ThetaZeroLocation = 'bottom';
hold on

% Cerchio indeformato
polarplot(theta_dense, circle, 'k--', 'LineWidth', 1.2)
polarplot(-theta_dense, circle, 'k--', 'LineWidth', 1.2,'HandleVisibility','off')

% Forma modale (deformata aggiunta al cerchio)
polarplot(theta_dense, r0 + mode_shape_interp, 'r', 'LineWidth', 2)

% Forma simmetrica (opzionale)
polarplot(-theta_dense, r0 - mode_sym_interp, 'b', 'LineWidth', 2)

polarplot(theta, r0 + mode_shape, 'or')

legend('Undeformed','Axial mode shape (identified)', 'Axial mode shape (symmetry)', 'Location','best')
title('Modal Shape in Polar Coordinates f = 1625.3 xi = 0.0075')