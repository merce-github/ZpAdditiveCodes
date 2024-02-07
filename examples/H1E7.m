print "Example: H1E7";
ei := GetEchoInput();
SetEchoInput(true);

ZpSimplexAlphaCode(2, 2, 3) eq SimplexAlphaCodeZ4(3);

ZpSimplexBetaCode(2, 2, 2) eq SimplexBetaCodeZ4(2);
S3, G3 := ZpSimplexBetaCode(2, 2, 3);
S3Z4 := SimplexBetaCodeZ4(3);
S3 eq S3Z4;
perm := Sym(28)!(17,21)(19,22)(23,27)(25,28);
G3b := G3^perm;
MultiplyColumn(~G3b, -1, 20);
MultiplyColumn(~G3b, -1, 26);
S3Z4 eq LinearCode(G3b);

SetEchoInput(ei); 