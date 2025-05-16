function err = cost_function1(params,freq,G_exp)
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