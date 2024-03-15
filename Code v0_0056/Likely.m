%Likelihood function following Uturbey
function [L e] = Likely(y,ar_params, ma_params, sigma_e,p_max)
    ar_order = size(ar_params,1);
    ma_order = size(ma_params,1);
    y_lagged = zeros(size(y,1),ar_order);
    for k = 1:ar_order
        y_lagged(k+1:end, k) = y(1:end-k);
    end
        
    e = y-y_lagged*ar_params;
    %compute log-likelihood
    %DN -> AMG: Testing log-Likelihood with fixed p at p_max
    L = (-(size(y,1)-p_max)/2)*(log(2) + log(pi) + 2*log(sigma_e))-transpose(e)*e/(2*(sigma_e^2));
    %Log-Likelihood with variable p
    %L = (-(size(y,1)-ar_order)/2)*(log(2) + log(pi) + 2*log(sigma_e))-transpose(e)*e/(2*sigma_e^2);
    %This is normal Likelihood:
    %L = (2*pi*sigma_e^2)^(-(size(y,1)-ar_order)/2)*exp(transpose(e)*e/(2*sigma_e^2)); 

    return;
end