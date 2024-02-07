print "Example: H1E16";
ei := GetEchoInput();
SetEchoInput(true);

p := 3; s := 3; t := 2; u := 1;

C := ZpHadamardCode(p, [2,0,1]);
N := OverDimension(CarletGrayMap(C)(Random(C)));
d := MinimumHomogeneousWeight(C);
d eq (p-1)*N/p;
d eq Min([HomogeneousWeight(c) : c in C | c ne C!0]);
HomogeneousWeightDistribution(C);

C := ZpSimplexAlphaCode(p, s, t);
N := OverDimension(CarletGrayMap(C)(Random(C)));
d := MinimumHomogeneousWeight(C);
d eq (p-1)*N/p;
HomogeneousWeightDistribution(C);

C := ZpSimplexBetaCode(p, s, t);
N := OverDimension(CarletGrayMap(C)(Random(C)));
d := MinimumHomogeneousWeight(C);
d eq (p-1)*N/p;
HomogeneousWeightDistribution(C);

C := ZpMacDonaldAlphaCode(p, s, t, u);
MinimumHomogeneousWeight(C);
HomogeneousWeightDistribution(C);

C := ZpMacDonaldBetaCode(p, s, t, u);
MinimumHomogeneousWeight(C);
HomogeneousWeightDistribution(C);

SetEchoInput(ei); 