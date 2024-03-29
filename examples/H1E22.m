print "Example: H1E22";
ei := GetEchoInput();
SetEchoInput(true);

C := ZpHadamardCode(2, [2,0,0]);

p1 := Sym(8)!(1,3,5,7)(2,4,6,8);
p2 := Sym(8)!(1,5)(2,6)(3,7)(4,8);
p3 := Sym(8)!(1,7,5,3)(2,8,6,4);
S := [Sym(8)!1, p1, p2, p3];

p1 := Sym(32)!(1,9,17,25)(2,10,18,26)(3,11,19,27)(4,12,20,28)(5,13,21,29)(6,14,22,30)(7,15,23,31)(8,16,24,32);
p2 := Sym(32)!(1,17)(2,18)(3,19)(4,20)(5,21)(6,22)(7,23)(8,24)(9,25)(10,26)(11,27)(12,28)(13,29)(14,30)(15,31)(16,32);
p3 := Sym(32)!(1,25,17,9)(2,26,18,10)(3,27,19,11)(4,28,20,12)(5,29,21,13)(6,30,22,14)(7,31,23,15)(8,32,24,16);
Sp := [Sym(32)!1, p1, p2, p3];

I, Ip := ZpInformationSet(C);

SetVerbose("IsPDSetFlag", 2);

IsZpPermutationDecodeSet(C, I, S, 3);

IsZpPermutationDecodeSet(C, Ip, Sp, 3);

SetEchoInput(ei); 
