print "Example: H1E12";
ei := GetEchoInput();
SetEchoInput(true);

C := ZpHadamardCode(2, [3,1]);

G_Zp, type := ZpMinRowsGeneratorMatrix(C);
G_Z4, t2, t1 := MinRowsGeneratorMatrix(C);
G_Zp eq G_Z4;
LinearCode(G_Zp) eq LinearCode(G_Z4);
type eq [t1, t2];

H_Zp := ZpMinRowsParityCheckMatrix(C);
H_Z4 := MinRowsParityCheckMatrix(C);
H_Zp eq H_Z4;
LinearCode(H_Zp) eq LinearCode(H_Z4);

ZpDual(C) eq DualZ4(C);

SetEchoInput(ei); 