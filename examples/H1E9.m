print "Example: H1E9";
ei := GetEchoInput();
SetEchoInput(true);

C := RandomLinearCode(Integers(4), 10, 8);	
ZpResidueCode(C) eq BinaryResidueCode(C);
ZpTorsionCode(C, 1) eq ZpResidueCode(C);
ZpTorsionCode(C, 2) eq BinaryTorsionCode(C);

G := Matrix(Integers(8), [[1,0,7,6,5,4,3,2],
                            [0,1,2,3,4,5,6,7],
                            [0,0,2,2,4,4,6,6],
                            [0,0,0,4,0,4,0,4]]);
C := LinearCode(G);
Gdiv := Matrix(Integers(8), [[1,0,7,6,5,4,3,2],
                               [0,1,2,3,4,5,6,7],
                               [0,0,1,1,2,2,3,3],
                               [0,0,0,1,0,1,0,1]]);
G eq DiagonalMatrix(Integers(8), [1,1,2,4])*Gdiv;
Gbin := Matrix(GF(2), Gdiv);

CRes := ZpResidueCode(C);
CRes eq LinearCode(Matrix(GF(2), GeneratorMatrix(C)));

CTor1 := ZpTorsionCode(C, 1);
CTor2 := ZpTorsionCode(C, 2);
CTor3 := ZpTorsionCode(C, 3);
CTor1 eq CRes;
CTor1 eq LinearCode(RowSubmatrix(Gbin, 1, ZpType(C)[1]));
CTor2 eq LinearCode(RowSubmatrix(Gbin, 1, &+ZpType(C)[1..2]));
CTor3 eq LinearCode(RowSubmatrix(Gbin, 1, &+ZpType(C)[1..3]));
(CTor1 subset CTor2) and (CTor2 subset CTor3);

SetEchoInput(ei); 