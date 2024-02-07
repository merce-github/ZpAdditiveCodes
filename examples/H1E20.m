print "Example: H1E20";
ei := GetEchoInput();
SetEchoInput(true);

C := LinearCode<Integers(27), 5 | [[1,2,5,8,6],
                                   [0,3,6,12,21],
                                   [0,0,9,9,18]]>;
V, Vp, f, fp := ZpInformationSpace(C : IsSystematicEncoding := false);
mapGrayC := CarletGrayMap(C);
mapGrayV := CarletGrayMap(3, ZpType(C));
  
i := V![21,12,18];
ip := Vp![2,2,0,1,2,2];
ip eq mapGrayV(i);
Encoding(C, ip) eq Encoding(C, i);
c, cp := Encoding(C, i);
(c eq f(i)) and (cp eq fp(ip));
cp eq mapGrayC(c);

ibar := [i[1], i[2] div 3, i[3] div 9];
c eq Vector(ibar) * ZpMinRowsGeneratorMatrix(C);

L := Eltseq(i) cat Eltseq(V![13,3,9]); 
Encoding(C, L);
Encoding(C,L)[1] eq Encoding(C, i);

M := Matrix(Integers(27), 3, L);
Encoding(C, M) eq Encoding(C, L);

SetEchoInput(ei); 
