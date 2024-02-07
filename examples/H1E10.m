print "Example: H1E10";
ei := GetEchoInput();
SetEchoInput(true);

C := RandomZpAdditiveCode(3, 1000, [2^^10]);

time D := Dual(C);
time Dp := ZpDual(C);
D eq Dp;

G := GeneratorMatrix(C);
H := GeneratorMatrix(D);
IsZero(G * Transpose(H));

ZpTypeDual(C) eq ZpType(ZpDual(C));

SetEchoInput(ei); 