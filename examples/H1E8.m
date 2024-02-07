print "Example: H1E8";
ei := GetEchoInput();
SetEchoInput(true);

p := 2; s := 3; t := 2; u := 1;
Malpha, GMalpha := ZpMacDonaldAlphaCode(p, s, t, u);
Malpha eq LinearCode(GMalpha);
_, A := ZpSimplexAlphaCode(p, s, t);
GMalpha eq ColumnSubmatrixRange(A, (p^(s*u))+1, (p^s)^t);
ZpType(Malpha) eq [t] cat [0^^(s-1)];

N := OverDimension(CarletGrayMap(Malpha)(Random(Malpha)));
N eq Length(Malpha)*p^(s-1);
N eq p^(s*(t+1)-1)-p^(s*(u+1)-1);
#Malpha eq p^(s*t);

Mbeta, GMbeta := ZpMacDonaldBetaCode(p, s, t, u);
Mbeta eq LinearCode(GMbeta);
ZpType(Mbeta) eq ZpType(Malpha);

N := OverDimension(CarletGrayMap(Mbeta)(Random(Mbeta)));
N eq Length(Mbeta)*p^(s-1);
N eq p^(s*t-t)*(p^t-1)/(p-1)-p^(s*u-u)*(p^u-1)/(p-1);
#Mbeta eq p^(s*t);

SetEchoInput(ei); 