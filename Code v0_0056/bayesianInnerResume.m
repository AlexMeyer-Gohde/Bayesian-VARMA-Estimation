function [states, accepted, dropped, modelinfo, parameters] = bayesianInnerResume(y, priorsARMA, proposalsARMA, settings, modelinfo, parameters, states, start)

%Start: The first new draw
progressbar

%Iterate until settings.draws is reached
for i = start:settings.draws
    [state, draw, modelinfo, parameters] = MH_Step_Inner(y, priorsARMA, proposalsARMA, states(i-1), settings, modelinfo, parameters);
    accepted = accepted + draw.accepted;
    dropped = dropped + draw.drawDropped;
    states(i) = state; 
    progressbar(i/settings.draws)
    if mod(i,500000) == 0
        save('temp.mat', '-v7.3');
    end;
%     if mod(i,1000) == 0
%         disp(['Iteration ' num2str(i)]);
%     end;
end;
end
