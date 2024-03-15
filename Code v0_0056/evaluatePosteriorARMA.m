function [logPosterior, estimatedVariance, dropDraw] = evaluatePosteriorARMA(y, state, priorsARMA, settings)
    %Returns the log Posterior
    %Currently only for one Process
    [logPosterior, dropDraw] = evaluatePriorARMA(state, priorsARMA, settings);
    if ~dropDraw
        for cntrProcess = 1:settings.processCount
            switch settings.likelihoodFunction
                case 0
                    modelinfo = struct;
                    modelinfo.data = y;
                    modelinfo.num_observations = max(size(y));
                    modelinfo.num_observables = 1;
                    modelinfo.num_exog = 1;

                    parameters = struct;
                    if state.ps > 0
                        parameters.W.Phi = transpose(state.arParameters);
                    else
                        parameters.W.Phi = [];
                    end;
                    if state.qs > 0
                        parameters.W.Theta = [1 transpose(state.maParameters)];
                    else
                        parameters.W.Theta = 1;
                    end;

                    parameters.W.p = state.ps;
                    parameters.W.q = state.qs;
                    parameters.Sigma = state.sigmaEs^2;

                    if max(size(parameters.W.Theta)) - 1 ~= parameters.W.q
                        disp('Something is not right');
                    end;

                    [logLikelihood, estimatedVariance] = evaluateLikelihoodARMAOnlyOLD(modelinfo,parameters);
                    logPosterior = logPosterior + logLikelihood;
                case 1                    
                    modelinfo = struct;
                    modelinfo.data = y;
                    modelinfo.num_observations = max(size(y));
                    modelinfo.num_observables = 1;
                    modelinfo.num_exog = 1;

                    parameters = struct;
                    if state.ps > 0
                        parameters.W.Phi = transpose(state.arParameters);
                    else
                        parameters.W.Phi = [];
                    end;
                    if state.qs > 0
                        parameters.W.Theta = [1 transpose(state.maParameters)];
                    else
                        parameters.W.Theta = 1;
                    end;

                    parameters.W.p = state.ps;
                    parameters.W.q = state.qs;
                    parameters.Sigma = state.sigmaEs^2;

                    if max(size(parameters.W.Theta)) - 1 ~= parameters.W.q
                        disp('Something is not right');
                    end;
                    [logLikelihood] = kalman_log_likelihood(parameters, modelinfo);
                    
%                     [logLikelihood1, estimatedVariance] = evaluateLikelihoodARMAOnlyOLD(modelinfo,parameters);
%                     logLikelihood2 = likelyARMA(y,state,settings);
%                     
%                     disp(['L Kalman: ' num2str(logLikelihood) ' L Spec: ' num2str(logLikelihood1) ' L Whack: ' num2str(logLikelihood2)]);
%                     disp(['Proportionality Kalman/Spec: ' num2str(logLikelihood - logLikelihood1) 'Proportionality Whack/Spec: ' num2str(logLikelihood2 - logLikelihood1)  'Proportionality Whack/Kalman: ' num2str(logLikelihood2 - logLikelihood)]);
                    logPosterior = logPosterior + logLikelihood;
                     estimatedVariance = -1/0;
                    
                case 2
                    logLikelihood = likelyARMA(y,state,settings);
                    logPosterior = logPosterior + logLikelihood;
                    estimatedVariance = 0;
                case 3
                    modelinfo = struct;
                    modelinfo.data = y;
                    modelinfo.num_observations = max(size(y));
                    modelinfo.num_observables = 1;
                    modelinfo.num_exog = 1;

                    parameters = struct;
                    if state.ps > 0
                        parameters.W.Phi = transpose(state.arParameters);
                    else
                        parameters.W.Phi = [];
                    end;
                    if state.qs > 0
                        parameters.W.Theta = [1 transpose(state.maParameters)];
                    else
                        parameters.W.Theta = 1;
                    end;

                    parameters.W.p = state.ps;
                    parameters.W.q = state.qs;
                    parameters.Sigma = state.sigmaEs^2;

                    if max(size(parameters.W.Theta)) - 1 ~= parameters.W.q
                        disp('Something is not right');
                    end;
%                     [logLikelihood1] = kalman_log_likelihood(parameters, modelinfo);
                    tic;
                    [logLikelihood1, estimatedVariance1] = evaluateLikelihoodARMAOnlyOLD(modelinfo,parameters);
                    toc;
                    logPosterior1 = logPosterior + logLikelihood1;        
                    tic;
                    [logLikelihood2] = kalman_log_likelihood_alt(parameters, modelinfo);
                    toc;
                    
                    tic;
                    [logLikelihood] = kalman_log_likelihood(parameters, modelinfo);
                    toc;
                    
                    disp(['LogLikelihood Spectral          : ' num2str(logLikelihood1) ]);                    
                    disp(['LogLikelihood Kalman            : ' num2str(logLikelihood) ]);                    
                    disp(['LogLikelihood Kalman Alternative: ' num2str(logLikelihood2) ]);                    

                    disp(['Diff LogLikelihood Spectral - Kalman               : ' num2str(logLikelihood1 - logLikelihood) ]);                    
                    disp(['Diff LogLikelihood Spectral - Kalman Alternative   : ' num2str(logLikelihood1 - logLikelihood2) ]);                    
                    disp(['Diff LogLikelihood Kalman - Kalman Alternative     : ' num2str(logLikelihood  - logLikelihood2) ]);                    
                    
                    if logLikelihood  - logLikelihood2 == 0
                        disp('That worked somehow....');
                    end;
                    estimatedVariance = 0;
%                     logLikelihood = likelyARMA(y,state,settings);
%                     estimatedVariance = 0;
                    logPosterior2 = logPosterior + logLikelihood;

                    
                    logPosterior = logPosterior2;
                    if (abs(logLikelihood - logLikelihood1) < 5)
                        disp('check this');
                    end;
                    
%                     if (isfinite(logLikelihood*logLikelihood1))
%                         disp(['Difference in Likelihood1: ' num2str(logLikelihood - logLikelihood1)]);
%                     else                        
%                         disp('boom');
%                     end;
%                     disp(['Difference in Likelihood2: ' num2str(L - logLikelihood)]);
%                     disp(['Difference in Likelihood3: ' num2str(L - logLikelihood1)]);
%                     disp(['Difference in Variance1: ' num2str(estimatedVariance1 - estimatedVariance)]);
                case 4
                    modelinfo = struct;
                    modelinfo.data = y;
                    modelinfo.num_observations = max(size(y));
                    modelinfo.num_observables = 1;
                    modelinfo.num_exog = 1;

                    parameters = struct;
                    if state.ps > 0
                        parameters.W.Phi = transpose(state.arParameters);
                    else
                        parameters.W.Phi = [];
                    end;
                    if state.qs > 0
                        parameters.W.Theta = [1 transpose(state.maParameters)];
                    else
                        parameters.W.Theta = 1;
                    end;

                    parameters.W.p = state.ps;
                    parameters.W.q = state.qs;
                    parameters.Sigma = state.sigmaEs^2;

                    if max(size(parameters.W.Theta)) - 1 ~= parameters.W.q
                        disp('Something is not right');
                    end;
                    [logLikelihood] = kalman_log_likelihood_alt(parameters, modelinfo);
                    
%                     [logLikelihood1, estimatedVariance] = evaluateLikelihoodARMAOnlyOLD(modelinfo,parameters);
%                     logLikelihood2 = likelyARMA(y,state,settings);
%                     
%                     disp(['L Kalman: ' num2str(logLikelihood) ' L Spec: ' num2str(logLikelihood1) ' L Whack: ' num2str(logLikelihood2)]);
%                     disp(['Proportionality Kalman/Spec: ' num2str(logLikelihood - logLikelihood1) 'Proportionality Whack/Spec: ' num2str(logLikelihood2 - logLikelihood1)  'Proportionality Whack/Kalman: ' num2str(logLikelihood2 - logLikelihood)]);
                    logPosterior = logPosterior + logLikelihood;
                    estimatedVariance = -1/0;
                case 5
                    modelinfo = struct;
                    modelinfo.data = y;
                    modelinfo.num_observations = max(size(y));
                    modelinfo.num_observables = 1;
                    modelinfo.num_exog = 1;
                    modelinfo.num_endo = 1;

                    parameters = struct;
                    if state.ps > 0
                        parameters.W.Phi = transpose(state.arParameters);
                    else
                        parameters.W.Phi = [];
                    end;
                    if state.qs > 0
                        parameters.W.Theta = [1 transpose(state.maParameters)];
                    else
                        parameters.W.Theta = 1;
                    end;

                    parameters.W.p = state.ps;
                    parameters.W.q = state.qs;
                    parameters.Sigma = state.sigmaEs^2;

                    if max(size(parameters.W.Theta)) - 1 ~= parameters.W.q
                        disp('Something is not right');
                    end;
                    [logLikelihood] = kalman_log_likelihood_alt_mv_onesidedhp(parameters, modelinfo, settings.HPLambda);
                    
                    logPosterior = logPosterior + logLikelihood;
                    estimatedVariance = -1/0;
            end;
        end;
    else
%          disp('Draw dropped due to zero prior');
        logPosterior = (-1/0);
        estimatedVariance = (-1/0);
    end;
end