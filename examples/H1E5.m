print "Example: H1E5";
ei := GetEchoInput();
SetEchoInput(true);

t1 := 2; t2 := 1;
C := ZpHadamardCode(2, [t1, t2]);
C eq HadamardCodeZ4(t1, t2 + 2*t1 - 1);

Cbin := CarletGrayMapImage(C);
N := OverDimension(Cbin[1]);
minDistance := MinimumHomogeneousWeight(C);
(#Cbin eq 2*N) and (minDistance eq N/2);
HasLinearCarletGrayMapImage(C);

p := 3; t1 := 1; t2 := 1; t3 := 1;
C, Gc := ZpHadamardCode(p, [t1,t2,t3]);
Gc;
C eq LinearCode(Gc);
ZpType(C) eq [t1,t2,t3];

Cp := CarletGrayMapImage(C);
N := OverDimension(Cp[1]);
minDistance := MinimumHomogeneousWeight(C);
(#Cp eq p*N) and (minDistance eq (p-1)*N/p);
HasLinearCarletGrayMapImage(C);

SetEchoInput(ei); 