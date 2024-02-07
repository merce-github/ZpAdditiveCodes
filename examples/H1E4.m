print "Example: H1E4";
ei := GetEchoInput();
SetEchoInput(true);

C, G := RandomZpAdditiveCode(3, 16, [2,3,3]);
IsStandardFormMatrix(G);
LinearCode(G) eq C;
ZpType(C);

C := RandomZpAdditiveCode(2, 16, 3, 8);
ZpType(C);
C := RandomLinearCode(Integers(2^3), 16, 8);
ZpType(C);

SetEchoInput(ei); 