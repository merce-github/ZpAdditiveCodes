print "Example: H1E11";
ei := GetEchoInput();
SetEchoInput(true);

C := RandomZpAdditiveCode(3, 1000, [2^^10]);

G := ZpMinRowsGeneratorMatrix(C);
H := ZpMinRowsParityCheckMatrix(C);
IsZero(G * Transpose(H));
LinearCode(G) eq C;
LinearCode(H) eq ZpDual(C);

Nrows(G) eq &+ZpType(C);
Nrows(H) eq &+ZpTypeDual(C);

time H := ZpMinRowsParityCheckMatrix(C);
time H2 := ZpMinRowsGeneratorMatrix(ZpDual(C));
time H3 := ParityCheckMatrix(C);
H eq H2;
H2 eq H3;
H eq H3;
LinearCode(H2) eq ZpDual(C);
LinearCode(H3) eq ZpDual(C);

SetEchoInput(ei); 