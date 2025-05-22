function [omega0, csi_exp, Gjk, FRF_exp] = modal_parameters_simplified_4(FRF,range, omega)
%The function takes as an input:
% - The experimental FRF vector or matrix;
% - The vector of the considered omegas;
% - The number of desired eigwnfrequencies to be computed
% The output is the experimental values of:
% - resonance frequencies;
% - damping;
% - mode shape at the selected resonance.

omega0 = zeros(size(FRF, 2), 1);
csi_exp = zeros(size(FRF, 2), 1);
Gjk = zeros(size(FRF, 2), 1);
f_lower = abs(omega-ones(length(omega), 1).*range(1));
f_upper = abs(omega-ones(length(omega), 1).*range(end));
a = find(f_lower == min(f_lower));
b = find(f_upper == min(f_upper));

for ii = 1:size(FRF, 2)
    FRF_temp = FRF(:,ii);
    peak_FRF = max(abs(FRF_temp(a:b)));
    index_w0 = find(abs(FRF_temp) == peak_FRF);
    omega0(ii) = omega(index_w0);
    FRF_csi = abs(abs(FRF_temp)-ones.*(peak_FRF/sqrt(2)));
    csi_exp(ii) = (omega(FRF_csi == min(FRF_csi(index_w0:b)))^2 ...
        - omega(FRF_csi == min(FRF_csi(a:index_w0)))^2) ./ (4*(omega0(ii).^2));
    Gjk(ii) = -imag(FRF_temp(index_w0)*2*csi_exp(ii)*omega0(ii)^2);
end
FRF_exp = FRF(a:b, :);
FRF_exp = FRF_exp';
end

