%Filename: Give Filename to load saved Workspace from
disp('Starting to load data. This could take a while...');
load('PUT FILENAME HERE');
disp('Data loaded. Now let us get to work');

settings.draws = 4000000;
start = 3000000;

[states, accepted, dropped, modelinfo, parameters] = bayesianInnerResume(y, priorsARMA, proposalsARMA, settings, modelinfo, parameters, start, filename);

try 
saveStuff(y, states, priorsARMA, proposalsARMA, settings, AR_order, MA_order);
catch(whateverhappened)
disp('HDD full?');
end;
plotResults(states,settings, accepted, y);
