print "Example: H1E1";
ei := GetEchoInput();
SetEchoInput(true);

Z8 := Integers(8);
C := LinearCode<Z8, 8 | [1,1,1,1,2,2,2,2],
                        [0,1,2,3,4,5,6,7],
                        [0,0,2,2,4,4,6,6] >;	
C;
V4 := VectorSpace(GF(2), 4);
mapGray := CarletGrayMap(2, 3);
mapGray(0);
mapGray(1);
mapGray(2);
mapGray(3);
mapGray(4); mapGray(4) eq mapGray(0) + V4![1,1,1,1];
mapGray(5); mapGray(5) eq mapGray(1) + V4![1,1,1,1];
mapGray(6); mapGray(6) eq mapGray(2) + V4![1,1,1,1];
mapGray(7); mapGray(7) eq mapGray(3) + V4![1,1,1,1];

mapGrayC := CarletGrayMap(C);
C.1;
mapGrayC(C.1);
mapGrayC(C.1) eq VectorSpace(GF(2), 8*4)!&cat[Eltseq(mapGray(i)) : i in Eltseq(C.1)];

HasLinearCarletGrayMapImage(C);
mapGrayU := CarletGrayMap(UniverseCode(Z8, Length(C)));
(mapGrayC(C.1) + mapGrayC(C.2)) @@ mapGrayU in C;

Z27 := Integers(27);
D := LinearCode<Z27, 10 | [1, 0, 0, 0, 0, 19, 18, 14, 10, 24],
                            [0, 1, 0, 0, 1, 10, 23,  3,  8,  6],
                            [0, 0, 1, 0, 1, 14, 23, 11, 18,  4], 
                            [0, 0, 0, 1, 2, 25, 20, 17,  2,  1],
                            [0, 0, 0, 0, 3, 21,  0,  9, 21,  0]>;	
#D;
HasLinearCarletGrayMapImage(D);
time HasLinearCarletGrayMapImage(D : AlgMethod := "StarProductMemory");

E := LinearCode< Integers(27), 6 | [[1,0,3,7,18,1],
                                    [0,3,0,0,21,18],
                                    [0,0,9,0,0,0],
                                    [0,0,0,9,0,0]] >;
time IsLinear, GrayCode := HasLinearCarletGrayMapImage(E);
time IsLinear := HasLinearCarletGrayMapImage(E : AlgMethod := "BruteForce");
IsLinear; Length(GrayCode); Dimension(GrayCode); Alphabet(GrayCode);

SetEchoInput(ei); 