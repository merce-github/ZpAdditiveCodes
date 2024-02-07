/************************************************************/
/*                                                          */
/* Project name: Z/p^s-additive codes in MAGMA              */
/* Test file name: ZpCarletGrayMap_BB_test.m                */
/*                                                          */
/* Comments: Black-box tests for the functions              */
/*           CarletGrayMap, CarletGrayMapImage, and         */  
/*           HasLinearCarletGrayMapImage                    */
/*           included in the ZpAdditiveCodes_Core.m file    */
/*                                                          */                                                
/* Authors: Guillermo Mosse, Noam von Rotberg and           */
/*          Mercè Villanueva                                */
/*                                                          */
/* Revision version and last date: v1.0    2018/06/28       */
/*                                 v1.1    2018/07/23       */
/*                                 v1.2    2018/07/26       */
/*                                 v1.3    2021/08/21       */
/*                                 v1.4    2023/12/14       */
/*                                                          */
/************************************************************/

SetAssertions(true);
//SetAssertions(true);
//Alarm(30*60);

/****************************************************************/
print "test 1: Zero code over Z16 of length 10";

p := 2;
s := 4;
n := 10;
nZp := n*p^(s-1);
Zps := Integers(p^s);
V := VectorSpace(GF(p), p^(s-1));
Vn := VectorSpace(GF(p), nZp);
U1 := UniverseCode(Zps, 1);

C := ZeroCode(Zps, n);

// Test CarletGrayMap function
OutputMap := CarletGrayMap(p, s);
assert Domain(OutputMap) eq Zps;
assert Codomain(OutputMap) eq V;
assert OutputMap(Zps!0) eq V!0;
assert {OutputMap(a) : a in Zps} eq Set(ReedMullerCode(1, s-1));
allones := V![1^^p^(s-1)];
for a in [0..p^(s-1)-1] do
    assert OutputMap(Zps!a) eq OutputMap(Zps!(a+p^(s-1))) + allones;
    assert OutputMap(Zps!a) eq CarletGrayMap(U1)(U1![a]);
    assert OutputMap(Zps!a) @@ OutputMap eq Zps!a; 
end for;

OutputMap := CarletGrayMap(C);
assert Domain(OutputMap) eq C;
assert Codomain(OutputMap) eq Vn;
assert OutputMap(C!0) eq Vn!0;

// Test CarletGrayMapImage function
OutputMapImage := CarletGrayMapImage(C);
expectedOutputMapImage := Setseq(Set(ZeroCode(GF(p), nZp)));
assert #OutputMapImage eq #expectedOutputMapImage;
assert Set(OutputMapImage) eq Set(expectedOutputMapImage);

/****************************************************************/
print "test 2: Universe code over Z8 of length 3";

p := 2;
s := 3;
n := 3;
nZp := n*p^(s-1);
Zps := Integers(p^s);
V := VectorSpace(GF(p), p^(s-1));
Vn := VectorSpace(GF(p), nZp);
U1 := UniverseCode(Zps, 1);

C := UniverseCode(Zps, n);

// Test CarletGrayMap function
OutputMap := CarletGrayMap(p, s);
assert Domain(OutputMap) eq Zps;
assert Codomain(OutputMap) eq V;
assert OutputMap(Zps!0) eq V!0;
assert {OutputMap(a) : a in Zps} eq Set(ReedMullerCode(1, s-1));
allones := V![1^^p^(s-1)];
for a in [0..p^(s-1)-1] do
    assert OutputMap(Zps!a) eq OutputMap(Zps!(a+p^(s-1))) + allones;
    assert OutputMap(Zps!a) eq CarletGrayMap(U1)(U1![a]); 
    assert OutputMap(Zps!a) @@ OutputMap eq Zps!a;
end for;

OutputMap := CarletGrayMap(C);
assert Domain(OutputMap) eq C;
assert Codomain(OutputMap) eq Vn;
assert OutputMap(C!0) eq Vn!0;

// Test CarletGrayMapImage function
OutputMapImage := CarletGrayMapImage(C);
expectedOutputMapImage := [OutputMap(c) : c in C];
assert #OutputMapImage eq #C;
assert Set(OutputMapImage) eq Set(expectedOutputMapImage);

/****************************************************************/
print "test 3: Repetition code over Z8 of length 5";

p := 2;
s := 3;
n := 5;
nbin := n*p^(s-1);
Z2s := Integers(p^s);
V := VectorSpace(GF(p), p^(s-1));
Vn := VectorSpace(GF(p), nbin);
U1 := UniverseCode(Z2s, 1);

C := RepetitionCode(Z2s, n);

// Test CarletGrayMap function
OutputMap := CarletGrayMap(p, s);
assert Domain(OutputMap) eq Z2s;
assert Codomain(OutputMap) eq V;
assert OutputMap(Z2s!0) eq V!0;
assert {OutputMap(a) : a in Z2s} eq Set(ReedMullerCode(1, s-1));
allones := V![1^^p^(s-1)];
for a in [0..p^(s-1)-1] do
    assert OutputMap(Z2s!a) eq OutputMap(Z2s!(a+p^(s-1))) + allones;
    assert OutputMap(Z2s!a) eq CarletGrayMap(U1)(U1![a]);
    assert OutputMap(Zps!a) @@ OutputMap eq Zps!a;
end for;

OutputMap := CarletGrayMap(C);
assert Domain(OutputMap) eq C;
assert Codomain(OutputMap) eq Vn;
assert OutputMap(C!0) eq Vn!0;

// Test CarletGrayMapImage function
OutputMapImage := CarletGrayMapImage(C);
expectedOutputMapImage := [OutputMap(c) : c in C];
assert #OutputMapImage eq #expectedOutputMapImage;
assert Set(OutputMapImage) eq Set(expectedOutputMapImage);

// Test ZpHadamardCode function cannot be tested with this code

/****************************************************************/
print "test 4: Z4-linear Hadamard code of type [1,0] and length 2";

p := 2;
L := [1, 0];
s := #L;
t := &+[(s-i+1)*L[i]: i in [1..s]] - 1;
n := p^(t-s+1);
nbin := n*p^(s-1); // It is equal to p^t
Z2s := Integers(p^s);
V := VectorSpace(GF(p), p^(s-1));
Vn := VectorSpace(GF(p), n*p^(s-1));

C, A := HadamardCodeZ4(L[1], t);

// Test CarletGrayMap function
OutputMap := CarletGrayMap(C);
assert Domain(OutputMap) eq C;
assert Codomain(OutputMap) eq Vn;
assert OutputMap(C!0) eq Vn!0;
for i in [0..10] do
    codeword := Random(C);
    assert OutputMap(codeword) @@ OutputMap eq codeword;
end for;

// Test CarletGrayMapImage function
OutputMapImage := CarletGrayMapImage(C);
assert Ncols(Random(OutputMapImage)) eq nbin; 
assert #OutputMapImage eq 2*nbin;
assert Min([Weight(c) : c in OutputMapImage | c ne Vn!0]) eq (nbin div 2);

// Test HasLinearCarletGrayMapImage function
OutputHasLinImage, OutputLinCode, OutputLinMap := HasLinearCarletGrayMapImage(C);
assert OutputHasLinImage;
assert OutputLinCode eq LinearCode(Matrix(GF(p), [[1,0],[0,1]]));
assert Domain(OutputLinMap) eq C;
assert Codomain(OutputLinMap) eq OutputLinCode;
assert [OutputLinMap(c) : c in C] eq [OutputMap(c) : c in C];
OutputHasLinImageBF, OutputLinCodeBF, OutputLinMapBF := HasLinearCarletGrayMapImage(C : 
                                                           AlgMethod := "BruteForce");
assert OutputHasLinImageBF eq OutputHasLinImage;
assert OutputLinCodeBF eq OutputLinCode;
assert [OutputLinMapBF(c) : c in C] eq [OutputMap(c) : c in C];
OutputHasLinImageSPM, OutputLinCodeSPM, OutputLinMapSPM := HasLinearCarletGrayMapImage(C : 
                                                           AlgMethod := "StarProductMemory");
assert OutputHasLinImageSPM eq OutputHasLinImage;
assert OutputLinCodeSPM eq OutputLinCode;
assert [OutputLinMapSPM(c) : c in C] eq [OutputMap(c) : c in C];

//Test ZpHadamardCode function
OutputCode, OutputMatrix := ZpHadamardCode(p, L);
assert OutputCode eq C;
assert OutputMatrix eq A;

//Test ZpType function
OutputType := ZpType(C);
assert OutputType eq L;

/****************************************************************/
print "test 5: Z4-linear Hadamard code of type [3,0] and length 2^5";

p := 2;
L := [3, 0];
s := #L;
t := &+[(s-i+1)*L[i]: i in [1..s]] - 1;
n := p^(t-s+1);
nbin := n*p^(s-1); // It is equal to p^t
Z2s := Integers(p^s);
V := VectorSpace(GF(p), p^(s-1));
Vn := VectorSpace(GF(p), n*p^(s-1));

C, A := HadamardCodeZ4(L[1], t);

// Test CarletGrayMap function
OutputMap := CarletGrayMap(C);
assert Domain(OutputMap) eq C;
assert Codomain(OutputMap) eq Vn;
assert OutputMap(C!0) eq Vn!0;
for i in [0..10] do
    codeword := Random(C);
    assert OutputMap(codeword) @@ OutputMap eq codeword;
end for;

// Test CarletGrayMapImage function
OutputMapImage := CarletGrayMapImage(C);
assert Ncols(Random(OutputMapImage)) eq nbin; 
assert #OutputMapImage eq 2*nbin;
assert Min([Weight(c) : c in OutputMapImage | c ne Vn!0]) eq (nbin div 2);

// Test HasLinearCarletGrayMapImage function
OutputHasLinImage, OutputLinCode, OutputLinMap := HasLinearCarletGrayMapImage(C);
assert not OutputHasLinImage;
assert OutputLinCode eq 0;
assert OutputLinMap eq 0;
OutputHasLinImageBF, OutputLinCodeBF, OutputLinMapBF := HasLinearCarletGrayMapImage(C : 
                                                           AlgMethod := "BruteForce");
assert OutputHasLinImageBF eq OutputHasLinImage;
assert OutputLinCodeBF eq 0;
assert OutputLinMapBF eq 0;
OutputHasLinImageSPM, OutputLinCodeSPM, OutputLinMapSPM := HasLinearCarletGrayMapImage(C : 
                                                           AlgMethod := "StarProductMemory");
assert OutputHasLinImageSPM eq OutputHasLinImage;
assert OutputLinCodeSPM eq 0;
assert OutputLinMapSPM eq 0;

//Test ZpHadamardCode function
OutputCode, OutputMatrix := ZpHadamardCode(p, L);
assert OutputCode eq C;
assert OutputMatrix eq A;

//Test ZpType function
OutputType := ZpType(C);
assert OutputType eq L;

/****************************************************************/
print "test 6: Z8-linear Hadamard code of type [1,0,1]
      (Example 1.1 of On Z2s-Linear Hadamard Codes: kernel and partial classification)";

p := 2;
L := [1, 0, 1];
s := #L;
t := &+[(s-i+1)*L[i]: i in [1..s]] - 1;
n := p^(t-s+1);
nbin := n*p^(s-1); // It is equal to p^t
Z2s := Integers(p^s);
V := VectorSpace(GF(p), p^(s-1));
Vn := VectorSpace(GF(p), n*p^(s-1));

M := [[1, 1], 
      [0, 4]];
A := Matrix(Z2s, M);
C := LinearCode(A);

// Test CarletGrayMap function
OutputMap := CarletGrayMap(C);
assert Domain(OutputMap) eq C;
assert Codomain(OutputMap) eq Vn;
assert OutputMap(C!0) eq Vn!0;
for i in [0..10] do
    codeword := Random(C);
    assert OutputMap(codeword) @@ OutputMap eq codeword;
end for;

// Test CarletGrayMapImage function
OutputMapImage := CarletGrayMapImage(C);
assert Ncols(Random(OutputMapImage)) eq nbin; 
assert #OutputMapImage eq 2*nbin;
assert Min([Weight(c) : c in OutputMapImage | c ne Vn!0]) eq (nbin div 2);

// Test HasLinearCarletGrayMapImage function
OutputHasLinImage, OutputLinCode, OutputLinMap := HasLinearCarletGrayMapImage(C);
assert OutputHasLinImage;
assert OutputLinCode eq LinearCode(Matrix([OutputMap(C![1,1]), OutputMap(C![2,2]),
                                           OutputMap(C![4,4]), OutputMap(C![0,4])]));
assert Domain(OutputLinMap) eq C;
assert Codomain(OutputLinMap) eq OutputLinCode;
assert [OutputLinMap(c) : c in C] eq [OutputMap(c) : c in C];
OutputHasLinImageBF, OutputLinCodeBF, OutputLinMapBF := HasLinearCarletGrayMapImage(C : 
                                                           AlgMethod := "BruteForce");
assert OutputHasLinImageBF eq OutputHasLinImage;
assert OutputLinCodeBF eq OutputLinCode;
assert [OutputLinMapBF(c) : c in C] eq [OutputMap(c) : c in C];
OutputHasLinImageSPM, OutputLinCodeSPM, OutputLinMapSPM := HasLinearCarletGrayMapImage(C : 
                                                           AlgMethod := "StarProductMemory");
assert OutputHasLinImageSPM eq OutputHasLinImage;
assert OutputLinCodeSPM eq OutputLinCode;
assert [OutputLinMapSPM(c) : c in C] eq [OutputMap(c) : c in C];

//Test ZpHadamardCode function
OutputCode, OutputMatrix := ZpHadamardCode(p, L);
assert OutputCode eq C;
assert OutputMatrix eq A;

//Test ZpType function
OutputType := ZpType(C);
assert OutputType eq L;

/****************************************************************/
print "test 7: Z8-linear Hadamard code of type [1,1,0] 
       (Example 1.2 of On Z2s-Linear Hadamard Codes: kernel and partial classification)";

p := 2;
L := [1, 1, 0];
s := #L;
t := &+[(s-i+1)*L[i]: i in [1..s]] - 1;
n := p^(t-s+1);
nbin := n*p^(s-1); // It is equal to p^t
Z2s := Integers(p^s);
V := VectorSpace(GF(p), p^(s-1));
Vn := VectorSpace(GF(p), n*p^(s-1));

M := [[1, 1, 1, 1], 
      [0, 2, 4, 6]];
A := Matrix(Z2s, M);
C := LinearCode(A);

// Test CarletGrayMap function
OutputMap := CarletGrayMap(C);
assert Domain(OutputMap) eq C;
assert Codomain(OutputMap) eq Vn;
assert OutputMap(C!0) eq Vn!0;
for i in [0..10] do
    codeword := Random(C);
    assert OutputMap(codeword) @@ OutputMap eq codeword;
end for;

// Test CarletGrayMapImage function
OutputMapImage := CarletGrayMapImage(C);
assert Ncols(Random(OutputMapImage)) eq nbin; 
assert #OutputMapImage eq 2*nbin;
assert Min([Weight(c) : c in OutputMapImage | c ne Vn!0]) eq (nbin div 2);

// Test HasLinearCarletGrayMapImage function
OutputHasLinImage, OutputLinCode, OutputLinMap := HasLinearCarletGrayMapImage(C);
assert OutputHasLinImage;
assert OutputLinCode eq LinearCode(Matrix([OutputMap(C![1,1,1,1]), 
                                           OutputMap(C![2,2,2,2]),
                                           OutputMap(C![4,4,4,4]), 
                                           OutputMap(C![0,2,4,6]),
                                           OutputMap(C![0,4,0,4])]));
assert Domain(OutputLinMap) eq C;
assert Codomain(OutputLinMap) eq OutputLinCode;
assert [OutputLinMap(c) : c in C] eq [OutputMap(c) : c in C];
OutputHasLinImageBF, OutputLinCodeBF, OutputLinMapBF := HasLinearCarletGrayMapImage(C : 
                                                           AlgMethod := "BruteForce");
assert OutputHasLinImageBF eq OutputHasLinImage;
assert OutputLinCodeBF eq OutputLinCode;
assert [OutputLinMapBF(c) : c in C] eq [OutputMap(c) : c in C];
OutputHasLinImageSPM, OutputLinCodeSPM, OutputLinMapSPM := HasLinearCarletGrayMapImage(C : 
                                                           AlgMethod := "StarProductMemory");
assert OutputHasLinImageSPM eq OutputHasLinImage;
assert OutputLinCodeSPM eq OutputLinCode;
assert [OutputLinMapSPM(c) : c in C] eq [OutputMap(c) : c in C];

//Test ZpHadamardCode function
OutputCode, OutputMatrix := ZpHadamardCode(p, L);
assert OutputCode eq C;
assert OutputMatrix eq A;

//Test ZpType function
OutputType := ZpType(C);
assert OutputType eq L;

/****************************************************************/
print "test 8: Z8-linear Hadamard code of type [2,0,0] 
       (Example 1.3 of On Z2s-Linear Hadamard Codes: kernel and partial classification)";

p := 2;
L := [2, 0, 0];
s := #L;
t := &+[(s-i+1)*L[i]: i in [1..s]] - 1;
n := p^(t-s+1);
nbin := n*p^(s-1); // It is equal to p^t
Z2s := Integers(p^s);
V := VectorSpace(GF(p), p^(s-1));
Vn := VectorSpace(GF(p), n*p^(s-1));

M := [[1,1,1,1, 1,1,1,1], 
      [0,1,2,3, 4,5,6,7]];
A := Matrix(Z2s, M);
C := LinearCode(A);

// Test CarletGrayMap function
OutputMap := CarletGrayMap(C);
assert Domain(OutputMap) eq C;
assert Codomain(OutputMap) eq Vn;
assert OutputMap(C!0) eq Vn!0;
for i in [0..10] do
    codeword := Random(C);
    assert OutputMap(codeword) @@ OutputMap eq codeword;
end for;

// Test CarletGrayMapImage function
OutputMapImage := CarletGrayMapImage(C);
assert Ncols(Random(OutputMapImage)) eq nbin; 
assert #OutputMapImage eq 2*nbin;
assert Min([Weight(c) : c in OutputMapImage | c ne Vn!0]) eq (nbin div 2);

// Test HasLinearCarletGrayMapImage function
OutputHasLinImage, OutputLinCode, OutputLinMap := HasLinearCarletGrayMapImage(C);
assert not OutputHasLinImage;
assert OutputLinCode eq 0;
assert OutputLinMap eq 0;
OutputHasLinImageBF, OutputLinCodeBF, OutputLinMapBF := HasLinearCarletGrayMapImage(C : 
                                                           AlgMethod := "BruteForce");
assert OutputHasLinImageBF eq OutputHasLinImage;
assert OutputLinCodeBF eq 0;
assert OutputLinMapBF eq 0;
OutputHasLinImageSPM, OutputLinCodeSPM, OutputLinMapSPM := HasLinearCarletGrayMapImage(C : 
                                                           AlgMethod := "StarProductMemory");
assert OutputHasLinImageSPM eq OutputHasLinImage;
assert OutputLinCodeSPM eq 0;
assert OutputLinMapSPM eq 0;

//Test ZpHadamardCode function
OutputCode, OutputMatrix := ZpHadamardCode(p, L);
assert OutputCode eq C;
assert OutputMatrix eq A;

//Test ZpType function
OutputType := ZpType(C);
assert OutputType eq L;

/****************************************************************/
print "test 9: Z8-linear Hadamard code of type [1,1,1] 
       (Example 1.4 of On Z2s-Linear Hadamard Codes: kernel and partial classification)";

p := 2;
L := [1, 1, 1];
s := #L;
t := &+[(s-i+1)*L[i]: i in [1..s]] - 1;
n := p^(t-s+1);
nbin := n*p^(s-1); // It is equal to p^t
Z2s := Integers(p^s);
V := VectorSpace(GF(p), p^(s-1));
Vn := VectorSpace(GF(p), n*p^(s-1));

M := [[1,1,1,1, 1,1,1,1], 
      [0,2,4,6, 0,2,4,6], 
      [0,0,0,0, 4,4,4,4]];
A := Matrix(Z2s, M);
C := LinearCode(A);

// Test CarletGrayMap function
OutputMap := CarletGrayMap(C);
assert Domain(OutputMap) eq C;
assert Codomain(OutputMap) eq Vn;
assert OutputMap(C!0) eq Vn!0;
for i in [0..10] do
    codeword := Random(C);
    assert OutputMap(codeword) @@ OutputMap eq codeword;
end for;

// Test CarletGrayMapImage function
OutputMapImage := CarletGrayMapImage(C);
assert Ncols(Random(OutputMapImage)) eq nbin; 
assert #OutputMapImage eq 2*nbin;
assert Min([Weight(c) : c in OutputMapImage | c ne Vn!0]) eq (nbin div 2);

// Test HasLinearCarletGrayMapImage function
OutputHasLinImage, OutputLinCode, OutputLinMap := HasLinearCarletGrayMapImage(C);
assert OutputHasLinImage;
assert OutputLinCode eq LinearCode(Matrix([OutputMap(C![1,1,1,1,1,1,1,1]), 
                                           OutputMap(C![2,2,2,2,2,2,2,2]),
                                           OutputMap(C![4,4,4,4,4,4,4,4]), 
                                           OutputMap(C![0,2,4,6,0,2,4,6]),
                                           OutputMap(C![0,4,0,4,0,4,0,4]),
                                           OutputMap(C![0,0,0,0,4,4,4,4])]));
assert Domain(OutputLinMap) eq C;
assert Codomain(OutputLinMap) eq OutputLinCode;
assert [OutputLinMap(c) : c in C] eq [OutputMap(c) : c in C];
OutputHasLinImageBF, OutputLinCodeBF, OutputLinMapBF := HasLinearCarletGrayMapImage(C : 
                                                           AlgMethod := "BruteForce");
assert OutputHasLinImageBF eq OutputHasLinImage;
assert OutputLinCodeBF eq OutputLinCode;
assert [OutputLinMapBF(c) : c in C] eq [OutputMap(c) : c in C];
OutputHasLinImageSPM, OutputLinCodeSPM, OutputLinMapSPM := HasLinearCarletGrayMapImage(C : 
                                                           AlgMethod := "StarProductMemory");
assert OutputHasLinImageSPM eq OutputHasLinImage;
assert OutputLinCodeSPM eq OutputLinCode;
assert [OutputLinMapSPM(c) : c in C] eq [OutputMap(c) : c in C];

//Test ZpHadamardCode function
OutputCode, OutputMatrix := ZpHadamardCode(p, L);
assert OutputCode eq C;
assert OutputMatrix eq A;

//Test ZpType function
OutputType := ZpType(C);
assert OutputType eq L;

/****************************************************************/
print "test 10: Z8-linear Hadamard code of type [2,0,0] 
       (Example 1.3 of On Z2s-Linear Hadamard Codes: kernel and partial classification)";

p := 2;
L := [2, 0, 0];
s := #L;
t := &+[(s-i+1)*L[i]: i in [1..s]] - 1;
n := p^(t-s+1);
nbin := n*p^(s-1); // It is equal to p^t
Z2s := Integers(p^s);
V := VectorSpace(GF(p), p^(s-1));
Vn := VectorSpace(GF(p), n*p^(s-1));

M := [[1,1,1,1, 1,1,1,1], 
      [0,1,2,3, 4,5,6,7]];
A := Matrix(Z2s, M);
C := LinearCode(A);

// Test CarletGrayMap function
OutputMap := CarletGrayMap(C);
assert Domain(OutputMap) eq C;
assert Codomain(OutputMap) eq Vn;
assert OutputMap(C!0) eq Vn!0;
for i in [0..10] do
    codeword := Random(C);
    assert OutputMap(codeword) @@ OutputMap eq codeword;
end for;

// Test CarletGrayMapImage function
OutputMapImage := CarletGrayMapImage(C);
assert Ncols(Random(OutputMapImage)) eq nbin; 
assert #OutputMapImage eq 2*nbin;
assert Min([Weight(c) : c in OutputMapImage | c ne Vn!0]) eq (nbin div 2);

// Test HasLinearCarletGrayMapImage function
OutputHasLinImage, OutputLinCode, OutputLinMap := HasLinearCarletGrayMapImage(C);
assert not OutputHasLinImage;
assert OutputLinCode eq 0;
assert OutputLinMap eq 0;
OutputHasLinImageBF, OutputLinCodeBF, OutputLinMapBF := HasLinearCarletGrayMapImage(C : 
                                                           AlgMethod := "BruteForce");
assert OutputHasLinImageBF eq OutputHasLinImage;
assert OutputLinCodeBF eq 0;
assert OutputLinMapBF eq 0;
OutputHasLinImageSPM, OutputLinCodeSPM, OutputLinMapSPM := HasLinearCarletGrayMapImage(C : 
                                                           AlgMethod := "StarProductMemory");
assert OutputHasLinImageSPM eq OutputHasLinImage;
assert OutputLinCodeSPM eq 0;
assert OutputLinMapSPM eq 0;

//Test ZpHadamardCode function
OutputCode, OutputMatrix := ZpHadamardCode(p, L);
assert OutputCode eq C;
assert OutputMatrix eq A;

//Test ZpType function
OutputType := ZpType(C);
assert OutputType eq L;

/****************************************************************/
print "test 11: Z8-linear Hadamard code of type [2,0,1] 
       (Example 1.5 of On Z2s-Linear Hadamard Codes: kernel and partial classification)";

p := 2;
L := [2, 0, 1];
s := #L;
t := &+[(s-i+1)*L[i]: i in [1..s]] - 1;
n := p^(t-s+1);
nbin := n*p^(s-1); // It is equal to p^t
Z2s := Integers(p^s);
V := VectorSpace(GF(p), p^(s-1));
Vn := VectorSpace(GF(p), n*p^(s-1));

M := [[1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1], 
      [0,1,2,3,4,5,6,7, 0,1,2,3,4,5,6,7], 
      [0,0,0,0,0,0,0,0, 4,4,4,4,4,4,4,4]];
A := Matrix(Z2s, M);
C := LinearCode(A);

// Test CarletGrayMap function
OutputMap := CarletGrayMap(C);
assert Domain(OutputMap) eq C;
assert Codomain(OutputMap) eq Vn;
assert OutputMap(C!0) eq Vn!0;
for i in [0..10] do
    codeword := Random(C);
    assert OutputMap(codeword) @@ OutputMap eq codeword;
end for;

// Test CarletGrayMapImage function
OutputMapImage := CarletGrayMapImage(C);
assert Ncols(Random(OutputMapImage)) eq nbin; 
assert #OutputMapImage eq 2*nbin;
assert Min([Weight(c) : c in OutputMapImage | c ne Vn!0]) eq (nbin div 2);

// Test HasLinearCarletGrayMapImage function
OutputHasLinImage, OutputLinCode, OutputLinMap := HasLinearCarletGrayMapImage(C);
assert not OutputHasLinImage;
assert OutputLinCode eq 0;
assert OutputLinMap eq 0;
OutputHasLinImageBF, OutputLinCodeBF, OutputLinMapBF := HasLinearCarletGrayMapImage(C : 
                                                           AlgMethod := "BruteForce");
assert OutputHasLinImageBF eq OutputHasLinImage;
assert OutputLinCodeBF eq 0;
assert OutputLinMapBF eq 0;
OutputHasLinImageSPM, OutputLinCodeSPM, OutputLinMapSPM := HasLinearCarletGrayMapImage(C : 
                                                           AlgMethod := "StarProductMemory");
assert OutputHasLinImageSPM eq OutputHasLinImage;
assert OutputLinCodeSPM eq 0;
assert OutputLinMapSPM eq 0;

//Test ZpHadamardCode function
OutputCode, OutputMatrix := ZpHadamardCode(p, L);
assert OutputCode eq C;
assert OutputMatrix eq A;

//Test ZpType function
OutputType := ZpType(C);
assert OutputType eq L;

/****************************************************************/
print "test 12: Z8-linear Hadamard code of type [2,1,0] 
       (Example 1.6 of On Z2s-Linear Hadamard Codes: kernel and partial classification)";

p := 2;
L := [2, 1, 0];
s := #L;
t := &+[(s-i+1)*L[i]: i in [1..s]] - 1;
n := p^(t-s+1);
nbin := n*p^(s-1); // It is equal to p^t
Z2s := Integers(p^s);
V := VectorSpace(GF(p), p^(s-1));
Vn := VectorSpace(GF(p), n*p^(s-1));

M := [[1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1],
      [0,1,2,3,4,5,6,7, 0,1,2,3,4,5,6,7, 0,1,2,3,4,5,6,7, 0,1,2,3,4,5,6,7],
      [0,0,0,0,0,0,0,0, 2,2,2,2,2,2,2,2, 4,4,4,4,4,4,4,4, 6,6,6,6,6,6,6,6]];
A := Matrix(Z2s, M);
C := LinearCode(A);

// Test CarletGrayMap function
OutputMap := CarletGrayMap(C);
assert Domain(OutputMap) eq C;
assert Codomain(OutputMap) eq Vn;
assert OutputMap(C!0) eq Vn!0;
for i in [0..10] do
    codeword := Random(C);
    assert OutputMap(codeword) @@ OutputMap eq codeword;
end for;

// Test CarletGrayMapImage function
OutputMapImage := CarletGrayMapImage(C);
assert Ncols(Random(OutputMapImage)) eq nbin; 
assert #OutputMapImage eq 2*nbin;
assert Min([Weight(c) : c in OutputMapImage | c ne Vn!0]) eq (nbin div 2);

// Test HasLinearCarletGrayMapImage function
OutputHasLinImage, OutputLinCode, OutputLinMap := HasLinearCarletGrayMapImage(C);
assert not OutputHasLinImage;
assert OutputLinCode eq 0;
assert OutputLinMap eq 0;
OutputHasLinImageBF, OutputLinCodeBF, OutputLinMapBF := HasLinearCarletGrayMapImage(C : 
                                                           AlgMethod := "BruteForce");
assert OutputHasLinImageBF eq OutputHasLinImage;
assert OutputLinCodeBF eq 0;
assert OutputLinMapBF eq 0;
OutputHasLinImageSPM, OutputLinCodeSPM, OutputLinMapSPM := HasLinearCarletGrayMapImage(C : 
                                                           AlgMethod := "StarProductMemory");
assert OutputHasLinImageSPM eq OutputHasLinImage;
assert OutputLinCodeSPM eq 0;
assert OutputLinMapSPM eq 0;

//Test ZpHadamardCode function
OutputCode, OutputMatrix := ZpHadamardCode(p, L);
assert OutputCode eq C;
assert OutputMatrix eq A;

//Test ZpType function
OutputType := ZpType(C);
assert OutputType eq L;

/****************************************************************/
print "test 13: Z4-linear Hadamard code of type [3,2]";

p := 2;
L := [3,2];
s := #L;
t := &+[(s-i+1)*L[i]: i in [1..s]] - 1;
n := p^(t-s+1);
nbin := n*p^(s-1); // It is equal to p^t
Z2s := Integers(p^s);
V := VectorSpace(GF(p), p^(s-1));
Vn := VectorSpace(GF(p), n*p^(s-1));

C, A := ZpHadamardCode(2, L);

// Test CarletGrayMap function
OutputMap := CarletGrayMap(C);
assert Domain(OutputMap) eq C;
assert Codomain(OutputMap) eq Vn;
assert OutputMap(C!0) eq Vn!0;
for i in [0..10] do
    codeword := Random(C);
    assert OutputMap(codeword) @@ OutputMap eq codeword;
end for;

// Test CarletGrayMapImage function
OutputMapImage := CarletGrayMapImage(C);
assert Ncols(Random(OutputMapImage)) eq nbin; 
assert #OutputMapImage eq 2*nbin;
assert Min([Weight(c) : c in OutputMapImage | c ne Vn!0]) eq (nbin div 2);

// Test HasLinearCarletGrayMapImage function
OutputHasLinImage, OutputLinCode, OutputLinMap := HasLinearCarletGrayMapImage(C);
assert not OutputHasLinImage;
assert OutputLinCode eq 0;
assert OutputLinMap eq 0;
OutputHasLinImageBF, OutputLinCodeBF, OutputLinMapBF := HasLinearCarletGrayMapImage(C : 
                                                           AlgMethod := "BruteForce");
assert OutputHasLinImageBF eq OutputHasLinImage;
assert OutputLinCodeBF eq 0;
assert OutputLinMapBF eq 0;
OutputHasLinImageSPM, OutputLinCodeSPM, OutputLinMapSPM := HasLinearCarletGrayMapImage(C : 
                                                           AlgMethod := "StarProductMemory");
assert OutputHasLinImageSPM eq OutputHasLinImage;
assert OutputLinCodeSPM eq 0;
assert OutputLinMapSPM eq 0;

//Test ZpType function
OutputType := ZpType(C);
assert OutputType eq L;

/****************************************************************/
print "test 14: Z8-linear Hadamard code of type [2,1,3]";

p := 2;
L := [2,1,3];
s := #L;
t := &+[(s-i+1)*L[i]: i in [1..s]] - 1;
n := p^(t-s+1);
nbin := n*p^(s-1); // It is equal to p^t
Z2s := Integers(p^s);
V := VectorSpace(GF(p), p^(s-1));
Vn := VectorSpace(GF(p), n*p^(s-1));

C, A := ZpHadamardCode(2, L);

// Test CarletGrayMap function
OutputMap := CarletGrayMap(C);
assert Domain(OutputMap) eq C;
assert Codomain(OutputMap) eq Vn;
assert OutputMap(C!0) eq Vn!0;
for i in [0..10] do
    codeword := Random(C);
    assert OutputMap(codeword) @@ OutputMap eq codeword;
end for;

// Test CarletGrayMapImage function
OutputMapImage := CarletGrayMapImage(C);
assert Ncols(Random(OutputMapImage)) eq nbin; 
assert #OutputMapImage eq 2*nbin;
assert Min([Weight(c) : c in OutputMapImage | c ne Vn!0]) eq (nbin div 2);

// Test HasLinearCarletGrayMapImage function
OutputHasLinImage, OutputLinCode, OutputLinMap := HasLinearCarletGrayMapImage(C);
assert not OutputHasLinImage;
assert OutputLinCode eq 0;
assert OutputLinMap eq 0;
OutputHasLinImageBF, OutputLinCodeBF, OutputLinMapBF := HasLinearCarletGrayMapImage(C : 
                                                           AlgMethod := "BruteForce");
assert OutputHasLinImageBF eq OutputHasLinImage;
assert OutputLinCodeBF eq 0;
assert OutputLinMapBF eq 0;
OutputHasLinImageSPM, OutputLinCodeSPM, OutputLinMapSPM := HasLinearCarletGrayMapImage(C : 
                                                           AlgMethod := "StarProductMemory");
assert OutputHasLinImageSPM eq OutputHasLinImage;
assert OutputLinCodeSPM eq 0;
assert OutputLinMapSPM eq 0;

//Test ZpType function
OutputType := ZpType(C);
assert OutputType eq L;

/****************************************************************/
print "test 15: Z16-linear Hadamard code of type [1,0,1,2]";

p := 2;
L := [1,0,1,2];
s := #L;
t := &+[(s-i+1)*L[i]: i in [1..s]] - 1;
n := p^(t-s+1);
nbin := n*p^(s-1); // It is equal to p^t
Z2s := Integers(p^s);
V := VectorSpace(GF(p), p^(s-1));
Vn := VectorSpace(GF(p), n*p^(s-1));

C, A := ZpHadamardCode(2, L);

// Test CarletGrayMap function
OutputMap := CarletGrayMap(C);
assert Domain(OutputMap) eq C;
assert Codomain(OutputMap) eq Vn;
assert OutputMap(C!0) eq Vn!0;
for i in [0..10] do
    codeword := Random(C);
    assert OutputMap(codeword) @@ OutputMap eq codeword;
end for;

// Test CarletGrayMapImage function
OutputMapImage := CarletGrayMapImage(C);
assert Ncols(Random(OutputMapImage)) eq nbin; 
assert #OutputMapImage eq 2*nbin;
assert Min([Weight(c) : c in OutputMapImage | c ne Vn!0]) eq (nbin div 2);

// Test HasLinearCarletGrayMapImage function
OutputHasLinImage, OutputLinCode, OutputLinMap := HasLinearCarletGrayMapImage(C);
assert OutputHasLinImage;
assert OutputLinCode eq LinearCode(Matrix([OutputMap(A[1]), OutputMap(2*A[1]), OutputMap(4*A[1]),
                                           OutputMap(A[2]), OutputMap(2*A[2]),
                                           OutputMap(RSpace(Z2s,n)![15,3,7,11,7,11,15,3,7,11,15,3,15,3,7,11]),
                                           OutputMap(A[3]), 
                                           OutputMap(A[4])]));
assert Domain(OutputLinMap) eq C;
assert Codomain(OutputLinMap) eq OutputLinCode;
assert [OutputLinMap(c) : c in C] eq [OutputMap(c) : c in C];
OutputHasLinImageBF, OutputLinCodeBF, OutputLinMapBF := HasLinearCarletGrayMapImage(C : 
                                                           AlgMethod := "BruteForce");
assert OutputHasLinImageBF eq OutputHasLinImage;
assert OutputLinCodeBF eq OutputLinCode;
assert [OutputLinMapBF(c) : c in C] eq [OutputMap(c) : c in C];
OutputHasLinImageSPM, OutputLinCodeSPM, OutputLinMapSPM := HasLinearCarletGrayMapImage(C : 
                                                           AlgMethod := "StarProductMemory");
assert OutputHasLinImageSPM eq OutputHasLinImage;
assert OutputLinCodeSPM eq OutputLinCode;
assert [OutputLinMapSPM(c) : c in C] eq [OutputMap(c) : c in C];

//Test ZpType function
OutputType := ZpType(C);
assert OutputType eq L;


