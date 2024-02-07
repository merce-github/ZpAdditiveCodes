print "Example: H1E2";
ei := GetEchoInput();
SetEchoInput(true);

Z8 := Integers(8);
C := LinearCode<Z8, 8 | [1,1,1,1,2,2,2,2],
                        [0,1,2,3,4,5,6,7],
                        [0,0,2,2,4,4,6,6] >;	
        
H := Matrix(GF(2), [[0,0,0,0], 
                    [0,1,0,1],
                    [0,0,1,1],
                    [0,1,1,0]]);
mapGrayH := GrayMap(H);
mapGrayCarlet := CarletGrayMap(2, 3);
HCarlet := Matrix([mapGrayCarlet(i) : i in [0..3]]);
H eq HCarlet;
a := Random(Z8); 
mapGrayH(a) eq mapGrayCarlet(a);

mapGrayC := GrayMap(C, H);
C.1;
mapGrayC(C.1);
C.2;
mapGrayC(C.2);

HasLinearGrayMapImage(C, H);
Cbin := GrayMapImage(C, H);
#Cbin eq #LinearCode(Matrix(Cbin));

H := Matrix(GF(2), [[0,0,0,0], 
                    [0,0,1,1],
                    [0,1,0,1],
                    [0,1,1,0]]);
mapGrayH := GrayMap(H);
H eq HCarlet;
mapGrayH(2) eq mapGrayCarlet(2);
mapGrayH(2) eq mapGrayCarlet(1);
HasLinearGrayMapImage(C, H);
Cbin := GrayMapImage(C, H);
#Cbin eq #LinearCode(Matrix(Cbin));
mapGrayC := GrayMap(C, H);
mapGrayU := GrayMap(UniverseCode(Z8, Length(C)), H);
(mapGrayC(C.1) + mapGrayC(C.2)) @@ mapGrayU in C;

SetEchoInput(ei); 