testing = zeros(10000000,1);
schubi = discreteLaplaceCDF(10,0.5);
for cntr=1:10000000 
     testing(cntr) = sampleDiscreteLaplace(5, schubi); 
end;