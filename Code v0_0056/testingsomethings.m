                    modelinfo = struct;
                    modelinfo.data = y;
                    modelinfo.num_observations = max(size(y));
                    modelinfo.num_observables = 1;
                    modelinfo.num_exog = 1;

                    p=3;
                    q=3;
                    
                    parameters = struct;
                    if p > 0
                        parameters.W.Phi = transpose([0.9; -0.23; 0.3]);
                    else
                        parameters.W.Phi = [];
                    end;
                    if q > 0
                        parameters.W.Theta = [1 transpose([-0.4; 0.6; -0.5])];
                    else
                        parameters.W.Theta = 1;
                    end;

                    parameters.W.p = p;
                    parameters.W.q = q;
                    parameters.Sigma = 1;

                    if max(size(parameters.W.Theta)) - 1 ~= parameters.W.q
                        disp('Something is not right');
                    end;

                    [logLikelihood, estimatedVariance] = evaluateLikelihoodARMAOnlyOLD(modelinfo,parameters);