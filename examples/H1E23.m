print "Example: H1E23";
ei := GetEchoInput();
SetEchoInput(true);

C := ZpHadamardCode(3, [2,0]);
n := Length(C);
np := n*3;

U := RSpace(Integers(9), n);
Up := VectorSpace(GF(3), np);
grayMap := CarletGrayMap(UniverseCode(Integers(9), n));

I, Ip := ZpInformationSet(C);

p1 := Sym(np)!(1,22,16,10,4,25,19,13,7)(2,23,17,11,5, 26,20,14,8)(3,24,18,12,6,27,21,15,9);
p2 := Sym(np)!(1,16,4,19,7,22,10,25,13)(2,17,5,20,8,23,11,26,14)(3,18,6,21,9,24,12,27,15);
p3 := Sym(np)!(1,10,19)(2,11,20)(3,12,21)(4,13,22)(5,14,23)(6,15,24)(7,16,25)(8,17,26)(9,18,27);
Sp := [Sym(np)!1, p1, p2, p3];

IsZpPermutationDecodeSet(C, Ip, Sp, 3);

c := C![0,1,2,3,4,5,6,7,8];
e := U![0,1,0,0,0,0,0,0,0];
u := c+e;
cp := grayMap(c);
up := grayMap(u);
(u in U) and (up in Up);

isDecoded, uDecoded, upDecoded := ZpPermutationDecode(C, Ip, Sp, 3, u);
isDecoded;
uDecoded eq c;
upDecoded eq cp;

isDecoded, uDecoded, upDecoded := ZpPermutationDecode(C, Ip, Sp, 3, up);
isDecoded;
uDecoded eq c;
upDecoded eq cp;

SetEchoInput(ei); 