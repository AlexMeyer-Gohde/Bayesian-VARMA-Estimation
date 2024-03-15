function [logPosterior, estimatedVariance, dropDraw, modelinfo, parameters] = evaluatePosterior(y, state, priorsARMA, settings, modelinfo, parameters)
    %Returns the log Posterior
    %Currently only for one Process
    [logPosterior, dropDraw] = evaluatePriorARMA(state, priorsARMA, settings);
    if ~dropDraw
        switch settings.likelihoodFunction
            case 0
                for cntrProcess = 1:settings.processCount
                    parameters.Z.p = state.ps;
                    if state.ps > 0 
                        parameters.Z.P = transpose(state.arParameters);
                    else
                        parameters.Z.P = [];
                    end;

                    parameters.Z.q = state.qs;
                    if state.qs > 0            
                        parameters.Z.Q = [transpose(state.maParameters)];
                    else
                        parameters.Z.Q = [];
                    end;

                    parameters.Z.Sigma = state.sigmaEs^2;

                    [modelinfo_solved,parameters_solved]=solve_recursive_dsge(modelinfo,parameters);

                    [logLikelihood, estimatedVariance] = evaluateLikelihood(modelinfo_solved,parameters_solved, settings.HPLambda);
                     logPosterior = logPosterior + logLikelihood;
                end;
            case 5
                for cntrProcess = 1:settings.processCount
                    parameters.Z.p = state.ps;
                    if state.ps > 0 
                        parameters.Z.P = transpose(state.arParameters);
                    else
                        parameters.Z.P = [];
                    end;

                    parameters.Z.q = state.qs;
                    if state.qs > 0            
                        parameters.Z.Q = [transpose(state.maParameters)];
                    else
                        parameters.Z.Q = [];
                    end;

                    parameters.Z.Sigma = state.sigmaEs^2;

                    [modelinfo_solved,parameters_solved]=solve_recursive_dsge(modelinfo,parameters);
                    modelinfo_solved.data = modelinfo_solved.data';

                    [logLikelihood] = kalman_log_likelihood_alt_mv_onesidedhp(modelinfo_solved,parameters_solved, settings.HPLambda);
                     logPosterior = logPosterior + logLikelihood;
                    estimatedVariance = (-1/0);       
                end;
            case 6
                for cntrProcess = 1:settings.processCount
                    parameters.Z.p = state.ps;
                    if state.ps > 0 
                        parameters.Z.P = transpose(state.arParameters);
                    else
                        parameters.Z.P = [];
                    end;

                    parameters.Z.q = state.qs;
                    if state.qs > 0            
                        parameters.Z.Q = [transpose(state.maParameters)];
                    else
                        parameters.Z.Q = [];
                    end;

                    parameters.Z.Sigma = state.sigmaEs^2;
                    
                    [modelinfo_solved,parameters_solved]=solve_recursive_dsge(modelinfo,parameters);
                    modelinfo_solved.data = modelinfo_solved.data';

                    [logLikelihood] = kalman_log_likelihood_alt_mv(modelinfo_solved,parameters_solved);
                     logPosterior = logPosterior + logLikelihood;
                    estimatedVariance = (-1/0);
                end;
            otherwise
                disp('Your choice of likelihood function is invalid');
                error('Your choice of likelihood function is invalid');
        end;
    else
%         disp('Draw dropped due to zero posterior');
        logPosterior = (-1/0);
        estimatedVariance = (-1/0);
    end;
end