% function [L] = likelyARMA(y,state,settings)
% %     following Hamilton 1994
% %   Set initial values to unconditional expecataion
% 
%     %Compute Innovations
%     maxorder = max(state.ps, state.qs);
%     N = size(y,1);
%     e = zeros(N+maxorder,1);
%     y = [zeros(maxorder,1);y];
%     y_pred = zeros(N+maxorder,1);
%     
%     if state.qs > 0 && state.ps > 0
%         for cntr = maxorder + 1 : N + maxorder
%             y_pred(cntr) = state.arParameters(1:state.ps)'*y(cntr-state.ps:cntr-1) + (state.maParameters(1:state.qs))'*e(cntr-state.qs:cntr-1);
%             e(cntr) = y(cntr) - y_pred(cntr);
%         end;
%     else
%         if state.qs == 0 && state.ps > 0
%             for cntr = maxorder + 1 : N + maxorder
%                 y_pred(cntr) = state.arParameters(1:state.ps)'*y(cntr-state.ps:cntr-1);
%                 e(cntr) = y(cntr) - y_pred(cntr);             
%             end;
%         else
%             if state.ps == 0 && state.qs > 0
%                 for cntr = maxorder + 1 : N + maxorder
%                     y_pred(cntr) = (state.maParameters(1:state.qs))'*e(cntr-state.qs:cntr-1);
%                     e(cntr) = y(cntr) - y_pred(cntr);
%                 end;
%             else
%                 e = y;
%             end;
%         end;
%     end;
% 
% %     disp(['error' num2str(sum(y - y_pred - e))]);    
%     e = e(end-settings.pMax-N+1:end);
%     %compute log-likelihood following hamilton 1994
%     disp(num2str(length(e) - N - settings.pMax));
%     L = (-(N-settings.pMax)/2)*(log(2) + log(pi) + 2*log(state.sigmaEs))-transpose(e)*e\(2*(state.sigmaEs^2));
% %     if abs(L) > 1000
% %         disp('sth is wrong');
% %     end;
%     %Log-Likelihood with variable p
%     %L = (-(size(y,1)-ar_order)/2)*(log(2) + log(pi) + 2*log(sigma_e))-transpose(e)*e/(2*sigma_e^2);
%     %This is normal Likelihood:
%     %L = (2*pi*sigma_e^2)^(-(size(y,1)-ar_order)/2)*exp(transpose(e)*e/(2*sigma_e^2)); 
% 
% end


%Likelihood function following Brooks, Ehlers 2004
%NOT FUNCTIONAL OLD VERSION
function [L] = likelyARMA(y,state,settings)
%     ar_order = size(ar_params,1);
%     ma_order = size(ma_params,1);
%     y_lagged = zeros(size(y,1),ar_order);
%     for k = 1:ar_order
%         y_lagged(k+1:end, k) = y(1:end-k);
%     end
%         
%     e = y-y_lagged*ar_params;

    %Compute Innovations
    N = size(y,1);
    e = zeros(N,1);
    y_pred = zeros(N-settings.pMax,1);
    
%     try
    if state.qs > 0 && state.ps > 0
        for cntr = settings.pMax + 1 : N
            y_pred(cntr) = state.arParameters(state.ps:-1:1)'*y(cntr-state.ps:cntr-1) + state.maParameters(state.qs:-1:1)'*e(cntr-state.qs:cntr-1);
            e(cntr) = y(cntr) - y_pred(cntr);
        end;
    else
        if state.qs == 0 && state.ps > 0
            for cntr = settings.pMax + 1 : N
                y_pred(cntr) = state.arParameters(state.ps:-1:1)'*y(cntr-state.ps:cntr-1);
                e(cntr) = y(cntr) - y_pred(cntr);             
            end;
        else
            if state.ps == 0 && state.qs > 0
                for cntr = settings.pMax + 1 : N
                    y_pred(cntr) = state.maParameters(state.qs:-1:1)'*e(cntr-state.qs:cntr-1);
                    e(cntr) = y(cntr) - y_pred(cntr);
                end;
            else
                e = y;
            end;
        end;
    end;
%     catch error
%         disp('schj8ubi');
%     end;
    
%     disp(['error' num2str(sum(y - y_pred - e))]);    
    e = e(settings.pMax + 1:end);
    %compute log-likelihood (Brooks Ehlers 2004)
%     disp(num2str(length(e) - (N - settings.pMax)));
    L = (-( N-settings.pMax)/2)*(log(2) + log(pi) + 2*log(state.sigmaEs))-transpose(e)*e/(2*(state.sigmaEs^2));
    if abs(L) > 1000
%         disp('sth is wrong');
    end;
    %Log-Likelihood with variable p
    %L = (-(size(y,1)-ar_order)/2)*(log(2) + log(pi) + 2*log(sigma_e))-transpose(e)*e/(2*sigma_e^2);
    %This is normal Likelihood:
    %L = (2*pi*sigma_e^2)^(-(size(y,1)-ar_order)/2)*exp(transpose(e)*e/(2*sigma_e^2)); 

end