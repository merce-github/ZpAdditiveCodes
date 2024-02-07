print "Example: H1E6";
ei := GetEchoInput();
SetEchoInput(true);

p := 2; s := 3; t := 2; 
Zps := Integers(p^s);
H, GH := ZpHadamardCode(p, [t+1] cat [0^^(s-1)]);

Salpha, GSalpha := ZpSimplexAlphaCode(p, s, t);
Salpha eq LinearCode(GSalpha);
GSalpha eq RemoveRow(GH, 1);   
_, A := ZpSimplexAlphaCode(p, s, t-1);
R1 := HorizontalJoin([A^^(p^s)]);
R2 := Matrix(Zps, [&cat[[i^^Ncols(A)] : i in [0..p^s-1]]]);
GSalpha eq VerticalJoin(R1, R2);
Length(Salpha) eq p^(s*t);
ZpType(Salpha) eq [t] cat [0^^(s-1)];  

Salpha_p := CarletGrayMapImage(Salpha);
N := OverDimension(Salpha_p[1]);
minDistance := MinimumHomogeneousWeight(Salpha);
N eq Length(Salpha)*p^(s-1);
#Salpha eq p^(s*t);
minDistance eq (p-1)*N/p;

Sbeta, GSbeta := ZpSimplexBetaCode(p, s, t);
Sbeta eq LinearCode(GSbeta);
_, A := ZpSimplexAlphaCode(p, s, t-1);
_, B := ZpSimplexBetaCode(p, s, t-1);
R1 := HorizontalJoin(A, HorizontalJoin([B^^(p^(s-1))]));
R2 := HorizontalJoin(Matrix(Zps, [[1^^Ncols(A)]]), 
                     Matrix(Zps, [[p*i : i in [0..p^(s-1)-1]]]));
GSbeta eq VerticalJoin(R1, R2); 
Length(Sbeta) eq p^((s-1)*(t-1))*(p^t-1)/(p-1);
ZpType(Sbeta) eq ZpType(Salpha);  

Sbeta_p := CarletGrayMapImage(Sbeta);
N := OverDimension(Sbeta_p[1]);
minDistance := MinimumHomogeneousWeight(Sbeta);
N eq Length(Sbeta)*p^(s-1);
#Sbeta eq p^(s*t);
minDistance eq (p-1)*N/p;

weightDistribution := {* Weight(c) : c in Sbeta_p *};
Multiplicity(weightDistribution, p^(s*t-t-1)*(p^t-1)) eq p^t*(p^((s-1)*t)-1);
Multiplicity(weightDistribution, p^(s*t-1)) eq p^t-1;
Multiplicity(weightDistribution, 0) eq 1;

SetEchoInput(ei); 