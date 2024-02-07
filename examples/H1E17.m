print "Example: H1E17";
ei := GetEchoInput();
SetEchoInput(true);

C := ZpHadamardCode(2, [2,1]);
v := Random(C);
HomogeneousWeight(v) eq LeeWeight(v); 
MinimumHomogeneousWeight(C) eq MinimumLeeWeight(C);
HomogeneousWeightDistribution(C) eq LeeWeightDistribution(C);      
DualHomogeneousWeightDistribution(C) eq DualLeeWeightDistribution(C); 

SetEchoInput(ei); 