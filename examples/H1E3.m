print "Example: H1E3";
ei := GetEchoInput();
SetEchoInput(true);

C := ZpHadamardCode(2, [1,0,0,0,2]);
HasLinearCarletGrayMapImage(C);
Cbin := CarletGrayMapImage(C);
#Cbin eq #LinearCode(Matrix(Cbin));

D := HadamardDatabase();
H := Matrix(D, 16, 2);
HasLinearGrayMapImage(C, H);
Cbin := GrayMapImage(C, H);
#Cbin eq #LinearCode(Matrix(Cbin));

SetEchoInput(ei); 