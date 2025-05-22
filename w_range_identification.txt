function [w] = w_range_identification(FRF, ww, W)

FRF_temp = abs(FRF(:, 1));
lim_inf   = 1;
lim_sup   = 0;
l = size(FRF, 1);
w = zeros(ww, l);

for ii = 1:ww
    peak_FRF = max(FRF_temp(lim_inf:end));
    index_w0 = find(FRF_temp == peak_FRF);

    der_FRF = diff(FRF);
    lim_vect = find(der_FRF(index_w0+1:end)>0);
    lim_sup = lim_sup + lim_vect(1);
    temp    = W(lim_inf:lim_sup);
    w(ii, 1:length(temp)) = temp;
    lim_inf = lim_sup+1;
end