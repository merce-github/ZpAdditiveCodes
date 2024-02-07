print "Example: H1E13";
ei := GetEchoInput();
SetEchoInput(true);

C := ZpHadamardCode(3, [1,0,2]);
S, f, Gs, perm := ZpStandardForm(C);
S;
(Domain(f) eq C) and (Codomain(f) eq S); 
f(Random(C)) in S;
Gs eq GeneratorMatrix(C)^perm;
Gs eq ZpMinRowsGeneratorMatrix(C)^perm;
Gs eq GeneratorMatrix(S);
IsStandardFormMatrix(GeneratorMatrix(C));
IsStandardFormMatrix(Gs);

Ds, fDual, Hs, permDual := ZpStandardFormDual(C);
(Domain(fDual) eq ZpDual(C)) and (Codomain(fDual) eq Ds);
fDual(Random(ZpDual(C))) in Ds;
Hs eq ZpMinRowsParityCheckMatrix(C)^permDual;
perm eq permDual;
(Ds eq LinearCode(Hs)) and (Ds eq ZpDual(S));
IsStandardFormMatrix(Hs);
IsZero(Hs * Transpose(Gs));

SetEchoInput(ei); 