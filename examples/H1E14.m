print "Example: H1E14";
ei := GetEchoInput();
SetEchoInput(true);

type := [1,2,3,1]; n := 20;
C := RandomZpAdditiveCode(3, n, type);
PseudoDimension(C) eq &+type;
ZpPseudoDimension(C) eq Log(3, #C);
ZpInformationRate(C) eq ZpPseudoDimension(C)/(n*#type);
ZpType(C) eq type;
ZpTypeDual(C) eq [n-&+type] cat Reverse(type[2..#type]);
ZpTypeDual(C) eq ZpType(ZpDual(C));
 
type := [1,0,1,0,2];
C := ZpHadamardCode(2, [1,0,1,0,2]);
PseudoDimension(C) eq &+type;
ZpPseudoDimension(C) eq Log(2, #C);
ZpInformationRate(C) eq ZpPseudoDimension(C)/(Length(C)*#type);
ZpType(C) eq type;
ZpTypeDual(C) eq [Length(C)-&+type] cat Reverse(type[2..#type]);
ZpTypeDual(C) eq ZpType(ZpDual(C));

SetEchoInput(ei); 