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

omega_1 = omega(1,1);
omega_2 = omega(2,1);

psi = 0.01;

%% Numerical computation

% First mode

f_min = 650;
f_max = 680;
freq1 = linspace(f_min,f_max,500);
for i = 1:12
    FRF_mod1(:,i) = interp1(f, FRF(:,i), freq1, 'spline');
end

Var = zeros(12,4);
FRF_num1 =zeros(12,length(freq1));
magnitude_num1 = zeros(12,length(freq1));
phase_num1 = zeros(12,length(freq1));

for i = 1:12
    params1 = [omega_1 psi mag(1,i) 0];
    err = @(params) cost_function1(params1,freq1,FRF_mod1(:,i));
    % err_v(i) = err;
    lb = [zeros(1,4)];
    ub = [Inf(1,4)];

    % Ottimizzazione
    opts = optimoptions('lsqnonlin','Display','iter','MaxFunctionEvaluations',5000);

    x_opt_1 = lsqnonlin(err, params1, lb, ub, opts);
    Var(i,:) = x_opt_1;
    FRF_num1_s = Var(i,3)./ (-(2*pi*freq1).^2 + 1i*2*Var(i,2)*Var(i,1)*2*pi*freq1 + Var(i,1)^2)+Var(i,4);
    FRF_num1(i,:) = FRF_num1_s;

    magnitude_num1(i,:) = abs(FRF_num1_s);
    phase_num1(i,:) = angle(FRF_num1_s)*(180/pi);
end

figure
subplot(2,1,1)
semilogy(f,magnitude(:,2),'LineWidth',3)
hold on 
semilogy(freq1,magnitude_num1(2,:),'or')
title('FRF')
xlabel('Frequency [Hz]')
ylabel('Magnitude [m/N]')
xlim([f_min f_max])
grid on

subplot(2,1,2)
plot(f,phase(:,2),'LineWidth',3)
hold on 
plot(freq1,phase_num1(2,:),'or')
xlabel ('Frequency [Hz]')
ylabel('Phase [°]')
xlim([f_min f_max])
grid on