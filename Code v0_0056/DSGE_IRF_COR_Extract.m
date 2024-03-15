function [IRF COR] = DSGE_IRF_COR_Extract(modelinfo,state,horizon,parameters,do,noninvertible)
parameters.Z.Sigma = state.sigmaEs^2;
%Extract orders
parameters.Z.p = state.ps;
if state.ps > 0
    parameters.Z.P = transpose(state.arParameters(1:parameters.Z.p,:));
else
    parameters.Z.P = [];
end;

parameters.Z.q = state.qs;
if state.qs > 0
    parameters.Z.Q = [transpose(state.maParameters(1:parameters.Z.q,:))];
    if noninvertible
        [parameters]=sample_noninvertible(parameters);
    end;
else
    parameters.Z.Q = [];
end;


%Solve DSGE Model
            
            [modelinfo_solved,parameters_solved]=solve_recursive_dsge(modelinfo,parameters);


% GET IRFs
if do.IRF
IRF = DSGE_IRF(modelinfo_solved, parameters_solved, horizon.IRF);
else
    IRF =[];
end
if do.COR
COR = DSGE_COR(modelinfo_solved, parameters_solved, horizon.COR, do.HPLAMBDA);
else
    COR =[];
end
