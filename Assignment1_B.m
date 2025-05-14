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
% mag=zeros(length(f),12);
% omega=zeros(length(f),12);
% for i=1:12
%     [mag(:,i), omega(:,i)]=findpeaks(magnitude(:,i),f,'MinPeakProminence', 0.05);
% end
% i_nat = [];
% 
% for i=2:length(f)-1
%     if magnitude(i,1) < magnitude(i-1,1) && magnitude(i,1) < magnitude(i+1,1)
%         i_nat(end+1) = i;
%     end
% end
% 
% omega = f (i_nat);
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


figure()
semilogy(f,magnitude(:,2))
grid on


