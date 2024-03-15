availableSeries = [ 1:40 50:100];

for cntrData = 1:length(availableSeries)
    load(['Results Bayesian_Monte Carlo Synthetic Philippe (3,2)_Sample Size 100_Draws 1500000_Date 08-Jan-2014_Data Series ' num2str(availableSeries(cntrData)) '_Chain 1 Save Type 2.mat'],'y');
    
    csvwrite(['Data_Chain_' num2str(availableSeries(cntrData)) '.csv']);
end;