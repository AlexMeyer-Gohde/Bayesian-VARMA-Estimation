function y = generateDataRBCUhligHansen(modelinfo, parameters)
        sampleSize = 250;

        z = zeros(sampleSize+1000,1);
        y = zeros(modelinfo.num_endo,sampleSize+1000);
        epsilon = normrnd(zeros(sampleSize+1000+ max(parameters.Z.p,parameters.Z.q)+1,1),sqrt(parameters.Z.Sigma));
        
        for genStep = (max(parameters.Z.p,parameters.Z.q)+1):sampleSize+1000+max(parameters.Z.p,parameters.Z.q)+1
            z(genStep) = parameters.Z.P*z(genStep - 1:-1:genStep - max(parameters.Z.p,1)) +...
                [1, parameters.Z.Q]*epsilon(genStep:-1:genStep-parameters.Z.q);
            y(:,genStep)=parameters.X.Lambda*y(:,genStep-1)+parameters.X.Phi*z(genStep :-1:genStep - max(parameters.Z.p,1)+1)...
                +parameters.X.Theta*epsilon(genStep:-1:genStep-parameters.Z.q+1);
        end;
        y=y(end-sampleSize+1:end);
end