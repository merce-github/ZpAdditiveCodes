/************************************************************/
/*                                                          */
/* Project name: Z/p^s-additive codes in MAGMA              */
/* Test file name: ZpHadamardGrayMap_BB_test.m              */
/*                                                          */
/* Comments: Black-box tests for the functions GrayMap(H),  */
/*           GrayMap(C,H) and GrayMapImage(C,H) included    */
/*           in the ZpAdditiveCodes_Core.m file             */
/*                                                          */
/* Authors: Javier Esmoris and M. Villanueva                */
/*                                                          */
/* Revision version and last date: v1.0    2023/07/25       */
/*                                                          */
/************************************************************/

SetAssertions(true);
//SetAssertions(true);
//Alarm(30*60);

/****************************************************************/
print "test 1: Zero code over Z32 of length 10";

p := 2;
s := 5;
n := 10;

nZp := n*p^(s-1);
Zps := Integers(p^s);
V := VectorSpace(GF(p), p^(s-1));
Vn := VectorSpace(GF(p), nZp);

C := ZeroCode(Zps, n);

CarletMap := CarletGrayMap(p, s);
H := Matrix(GF(p), [Eltseq(CarletMap(u)) : u in [Zps!i : i in [0..p^(s-1)-1]]]);

// Test GrayMap(H) function
OutputMap := GrayMap(H);
assert Domain(OutputMap) eq Zps;
assert Codomain(OutputMap) eq V;
assert OutputMap(Zps!0) eq V!0;
assert [OutputMap(a) : a in Zps] eq [CarletMap(a) : a in Zps];
for a in [0..p^(s-1)-1] do
    for lambda in [0..p-1] do
        lambdaVector := V![(p-lambda)^^p^(s-1)];
        assert OutputMap(Zps!a) eq OutputMap(Zps!(a+lambda*p^(s-1))) + lambdaVector;
    end for;
end for;

// Test GrayMap(C, H) function
OutputMapCode := GrayMap(C, H);
assert Domain(OutputMapCode) eq C;
assert Codomain(OutputMapCode) eq Vn;
assert OutputMapCode(C!0) eq Vn!0;
randomWord := Random(C);
assert OutputMapCode(randomWord) eq Vn!&cat[Eltseq(OutputMap(u)) : 
                                            u in Eltseq(randomWord)];

// Test GrayMapImage(C, H) function
OutputMapCodeImage := GrayMapImage(C, H);
expectedOutputMapCodeImage := CarletGrayMapImage(C);
assert #OutputMapCodeImage eq #expectedOutputMapCodeImage;
assert Set(OutputMapCodeImage) eq Set(expectedOutputMapCodeImage);

/****************************************************************/
print "test 2: Universe code over Z25 of length 3";

p := 5;
s := 2;
n := 3;

nZp := n*p^(s-1);
Zps := Integers(p^s);
V := VectorSpace(GF(p), p^(s-1));
Vn := VectorSpace(GF(p), nZp);

C := UniverseCode(Zps, n);

CarletMap := CarletGrayMap(p, s);
H := Matrix(GF(p), [Eltseq(CarletMap(u)) : u in [Zps!i : i in [0..p^(s-1)-1]]]);

// Test GrayMap(H) function
OutputMap := GrayMap(H);
assert Domain(OutputMap) eq Zps;
assert Codomain(OutputMap) eq V;
assert OutputMap(Zps!0) eq V!0;
assert {OutputMap(a) : a in Zps} eq {CarletMap(a) : a in Zps};
for a in [0..p^(s-1)-1] do
    for lambda in [0..p-1] do
        lambdaVector := V![(p-lambda)^^p^(s-1)];
        assert OutputMap(Zps!a) eq OutputMap(Zps!(a+lambda*p^(s-1))) + lambdaVector;
    end for;
end for;

// Test GrayMap(C, H) function
OutputMapCode := GrayMap(C, H);
assert Domain(OutputMapCode) eq C;
assert Codomain(OutputMapCode) eq Vn;
assert OutputMapCode(C!0) eq Vn!0;
randomWord := Random(C);
assert OutputMapCode(randomWord) eq Vn!&cat[Eltseq(OutputMap(u)) : 
                                            u in Eltseq(randomWord)];

// Test GrayMapImage(C, H) function
OutputMapCodeImage := GrayMapImage(C, H);
expectedOutputMapCodeImage := CarletGrayMapImage(C);
assert #OutputMapCodeImage eq #expectedOutputMapCodeImage;
assert Set(OutputMapCodeImage) eq Set(expectedOutputMapCodeImage);

/****************************************************************/
print "test 3: Repetition code over Z27 of length 5";

p := 3;
s := 3;
n := 5;

nZp := n*p^(s-1);
Zps := Integers(p^s);
V := VectorSpace(GF(p), p^(s-1));
Vn := VectorSpace(GF(p), nZp);

C := RepetitionCode(Zps, n);

CarletMap := CarletGrayMap(p, s);
H := Matrix(GF(p), p^(s-1) , p^(s-1), [ 
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 2, 0, 1, 2, 0, 1, 2],
    [0, 2, 1, 0, 2, 1, 0, 2, 1],
    [0, 0, 0, 1, 1, 1, 2, 2, 2],
    [0, 1, 2, 1, 2, 0, 2, 0, 1],
    [0, 2, 1, 1, 0, 2, 2, 1, 0],
    [0, 0, 0, 2, 2, 2, 1, 1, 1],
    [0, 1, 2, 2, 0, 1, 1, 2, 0],
    [0, 2, 1, 2, 1, 0, 1, 0, 2]
]);

// Test GrayMap(H) function
OutputMap := GrayMap(H);
assert Domain(OutputMap) eq Zps;
assert Codomain(OutputMap) eq V;
assert OutputMap(Zps!0) eq V!0;
assert {OutputMap(a) : a in Zps} eq {CarletMap(a) : a in Zps};
for a in [0..p^(s-1)-1] do
    for lambda in [0..p-1] do
        lambdaVector := V![(p-lambda)^^p^(s-1)];
        assert OutputMap(Zps!a) eq OutputMap(Zps!(a+lambda*p^(s-1))) + lambdaVector;
    end for;
end for;

// Test GrayMap(C, H) function
OutputMapCode := GrayMap(C, H);
assert Domain(OutputMapCode) eq C;
assert Codomain(OutputMapCode) eq Vn;
assert OutputMapCode(C!0) eq Vn!0;
randomWord := Random(C);
assert OutputMapCode(randomWord) eq Vn!&cat[Eltseq(OutputMap(u)) : 
                                            u in Eltseq(randomWord)];

// Test GrayMapImage(C, H) function
OutputMapCodeImage := GrayMapImage(C, H);
expectedOutputMapCodeImage := CarletGrayMapImage(C);
assert #OutputMapCodeImage eq #expectedOutputMapCodeImage;
assert Set(OutputMapCodeImage) eq Set(expectedOutputMapCodeImage);

/****************************************************************/
print "test 4: Linear code over Z49 of type (3; 2,0)";

p := 7;
s := 2;
n := 3;

nZp := n*p^(s-1);
Zps := Integers(p^s);
V := VectorSpace(GF(p), p^(s-1));
Vn := VectorSpace(GF(p), nZp);

G := [ [1,2,3],
       [2,1,4] ];
C := LinearCode(Matrix(Zps, G));

H := Matrix(GF(p), p^(s-1) , p^(s-1), [
    [ 0, 0, 0, 0, 0, 0, 0 ],
    [ 0, 4, 1, 5, 2, 6, 3 ],
    [ 0, 2, 4, 6, 1, 3, 5 ],
    [ 0, 5, 3, 1, 6, 4, 2 ],
    [ 0, 3, 6, 2, 5, 1, 4 ],
    [ 0, 6, 5, 4, 3, 2, 1 ],
    [ 0, 1, 2, 3, 4, 5, 6 ]
]);

// Test GrayMap(H) function
OutputMap := GrayMap(H);
assert Domain(OutputMap) eq Zps;
assert Codomain(OutputMap) eq V;
assert OutputMap(Zps!0) eq V!0;
assert [OutputMap(a) : a in [0..p^(s-1)-1]] eq Rows(H);
for a in [0..p^(s-1)-1] do
    for lambda in [0..p-1] do
        lambdaVector := V![(p-lambda)^^p^(s-1)];
        assert OutputMap(Zps!a) eq OutputMap(Zps!(a+lambda*p^(s-1))) + lambdaVector;
    end for;
end for;

// Test GrayMap(C, H) function
OutputMapCode := GrayMap(C, H);
assert Domain(OutputMapCode) eq C;
assert Codomain(OutputMapCode) eq Vn;
assert OutputMapCode(C!0) eq Vn!0;
randomWord := Random(C);
assert OutputMapCode(randomWord) eq Vn!&cat[Eltseq(OutputMap(u)) : 
                                            u in Eltseq(randomWord)];

// Test GrayMapImage(C, H) function
OutputMapCodeImage := GrayMapImage(C, H);
expectedOutputMapCodeImage := [OutputMapCode(c) : c in C];
assert #OutputMapCodeImage eq #expectedOutputMapCodeImage;
assert OutputMapCodeImage eq expectedOutputMapCodeImage;

/****************************************************************/
print "test 5: Linear code over Z32 of type (5; 2,0,0,0,0)";

p := 2;
s := 5;
n := 5;

nZp := n*p^(s-1);
Zps := Integers(p^s);
V := VectorSpace(GF(p), p^(s-1));
Vn := VectorSpace(GF(p), nZp);

G := [ [1,0,1,0,1],
       [0,1,2,4,8] ];
C := LinearCode(Matrix(Zps, G));

H := Matrix(GF(p), p^(s-1) , p^(s-1), [
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1],
    [0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1],
    [0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1],
    [0,0,1,1,1,0,0,1,0,1,1,0,1,1,0,0],
    [0,0,1,1,0,1,1,0,1,0,0,1,1,1,0,0],
    [0,0,1,1,1,1,0,0,1,1,0,0,0,0,1,1],
    [0,1,1,0,1,1,0,0,0,0,1,1,0,1,1,0],
    [0,1,1,0,0,0,1,1,1,1,0,0,0,1,1,0],
    [0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1],
    [0,1,0,1,0,1,1,0,0,1,1,0,1,0,1,0],
    [0,1,0,1,1,0,0,1,1,0,0,1,1,0,1,0],
    [0,1,0,1,1,0,1,0,1,0,1,0,0,1,0,1],
    [0,1,1,0,0,1,0,1,1,0,1,0,1,0,0,1],
    [0,1,1,0,1,0,1,0,0,1,0,1,1,0,0,1],
    [0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0]
]);

// Test GrayMap(H) function
OutputMap := GrayMap(H);
assert Domain(OutputMap) eq Zps;
assert Codomain(OutputMap) eq V;
assert OutputMap(Zps!0) eq V!0;
assert [OutputMap(a) : a in [0..p^(s-1)-1]] eq Rows(H);
for a in [0..p^(s-1)-1] do
    for lambda in [0..p-1] do
        lambdaVector := V![(p-lambda)^^p^(s-1)];
        assert OutputMap(Zps!a) eq OutputMap(Zps!(a+lambda*p^(s-1))) + lambdaVector;
    end for;
end for;

// Test GrayMap(C, H) function
OutputMapCode := GrayMap(C, H);
assert Domain(OutputMapCode) eq C;
assert Codomain(OutputMapCode) eq Vn;
assert OutputMapCode(C!0) eq Vn!0;
randomWord := Random(C);
assert OutputMapCode(randomWord) eq Vn!&cat[Eltseq(OutputMap(u)) : 
                                            u in Eltseq(randomWord)];

// Test GrayMapImage(C, H) function
OutputMapCodeImage := GrayMapImage(C, H);
expectedOutputMapCodeImage := [OutputMapCode(c) : c in C];
assert #OutputMapCodeImage eq #expectedOutputMapCodeImage;
assert OutputMapCodeImage eq expectedOutputMapCodeImage;

/****************************************************************/
print "test 6: Linear code over Z81 of type (4; 1,1,1,1)";

p := 3;
s := 4;
n := 4;

nZp := n*p^(s-1);
Zps := Integers(p^s);
V := VectorSpace(GF(p), p^(s-1));
Vn := VectorSpace(GF(p), nZp);

G := [ [1,1,1,1],
       [0,3,3,3],
       [0,0,9,9],
       [0,0,0,27] ];
C := LinearCode(Matrix(Zps, G));

H := Matrix(GF(p), p^(s-1) , p^(s-1), [
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,1,2,0,0,0,2,1,0,2,0,1,2,2,2,1,0,2,1,2,0,1,1,1,0,2,1],
    [0,2,1,0,0,0,2,0,1,1,0,2,1,1,1,0,1,2,2,1,0,2,2,2,1,2,0],
    [0,0,0,0,1,2,0,2,1,1,1,1,1,2,0,1,0,2,2,2,2,2,0,1,2,1,0],
    [0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2],
    [0,2,1,0,1,2,0,0,0,2,1,0,2,0,1,2,2,2,1,0,2,1,2,0,1,1,1],
    [0,0,0,0,2,1,1,2,0,2,2,2,2,1,0,0,1,2,1,1,1,1,0,2,2,0,1],
    [0,1,2,0,2,1,1,1,1,1,2,0,1,0,2,2,2,2,2,0,1,2,1,0,0,0,0],
    [0,2,1,0,2,1,0,2,1,0,2,1,0,2,1,0,2,1,0,2,1,0,2,1,0,2,1],
    [0,0,0,1,1,1,2,2,2,0,0,0,1,1,1,2,2,2,0,0,0,1,1,1,2,2,2],
    [0,1,2,1,1,1,1,0,2,2,0,1,0,0,0,0,2,1,1,2,0,2,2,2,2,1,0],
    [0,2,1,1,1,1,1,2,0,1,0,2,2,2,2,2,0,1,2,1,0,0,0,0,0,1,2],
    [0,0,0,1,2,0,2,1,0,1,1,1,2,0,1,0,2,1,2,2,2,0,1,2,1,0,2],
    [0,1,2,1,2,0,2,0,1,0,1,2,1,2,0,2,0,1,0,1,2,1,2,0,2,0,1],
    [0,2,1,1,2,0,2,2,2,2,1,0,0,1,2,1,1,1,1,0,2,2,0,1,0,0,0],
    [0,0,0,1,0,2,0,1,2,2,2,2,0,2,1,2,0,1,1,1,1,2,1,0,1,2,0],
    [0,1,2,1,0,2,0,0,0,1,2,0,2,1,0,1,1,1,2,0,1,0,2,1,2,2,2],
    [0,2,1,1,0,2,2,1,0,0,2,1,1,0,2,2,1,0,0,2,1,1,0,2,2,1,0],
    [0,0,0,2,2,2,1,1,1,0,0,0,2,2,2,1,1,1,0,0,0,2,2,2,1,1,1],
    [0,1,2,2,2,2,0,2,1,2,0,1,1,1,1,2,1,0,1,2,0,0,0,0,1,0,2],
    [0,2,1,2,2,2,0,1,2,1,0,2,0,0,0,1,2,0,2,1,0,1,1,1,2,0,1],
    [0,0,0,2,0,1,1,0,2,1,1,1,0,1,2,2,1,0,2,2,2,1,2,0,0,2,1],
    [0,1,2,2,0,1,1,2,0,0,1,2,2,0,1,1,2,0,0,1,2,2,0,1,1,2,0],
    [0,2,1,2,0,1,1,1,1,2,1,0,1,2,0,0,0,0,1,0,2,0,1,2,2,2,2],
    [0,0,0,2,1,0,2,0,1,2,2,2,1,0,2,1,2,0,1,1,1,0,2,1,0,1,2],
    [0,1,2,2,1,0,2,2,2,1,2,0,0,2,1,0,0,0,2,0,1,1,0,2,1,1,1],
    [0,2,1,2,1,0,1,0,2,0,2,1,2,1,0,1,0,2,0,2,1,2,1,0,1,0,2]
    ]);

// Test GrayMap(H) function
OutputMap := GrayMap(H);
assert Domain(OutputMap) eq Zps;
assert Codomain(OutputMap) eq V;
assert OutputMap(Zps!0) eq V!0;
assert [OutputMap(a) : a in [0..p^(s-1)-1]] eq Rows(H);
for a in [0..p^(s-1)-1] do
    for lambda in [0..p-1] do
        lambdaVector := V![(p-lambda)^^p^(s-1)];
        assert OutputMap(Zps!a) eq OutputMap(Zps!(a+lambda*p^(s-1))) + lambdaVector;
    end for;
end for;

// Test GrayMap(C, H) function
OutputMapCode := GrayMap(C, H);
assert Domain(OutputMapCode) eq C;
assert Codomain(OutputMapCode) eq Vn;
assert OutputMapCode(C!0) eq Vn!0;
randomWord := Random(C);
assert OutputMapCode(randomWord) eq Vn!&cat[Eltseq(OutputMap(u)) : 
                                            u in Eltseq(randomWord)];

// Test GrayMapImage(C, H) function
OutputMapCodeImage := GrayMapImage(C, H);
expectedOutputMapCodeImage := [OutputMapCode(c) : c in C];
assert #OutputMapCodeImage eq #expectedOutputMapCodeImage;
assert OutputMapCodeImage eq expectedOutputMapCodeImage;

/****************************************************************/
print "test 7: Linear code over Z16 of type (6; 0,3,0,0)";

p := 2;
s := 4;
n := 6;

nZp := n*p^(s-1);
Zps := Integers(p^s);
V := VectorSpace(GF(p), p^(s-1));
Vn := VectorSpace(GF(p), nZp);

G := [ [2,0,0,2,0,0],
       [0,2,0,0,2,0],
       [0,0,2,0,0,2] ];
C := LinearCode(Matrix(Zps, G));

H := Matrix(GF(p), p^(s-1) , p^(s-1), [
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,1,1,1,1],
    [0,0,1,1,0,0,1,1],
    [0,1,0,1,0,1,0,1],
    [0,1,0,1,1,0,1,0],
    [0,1,1,0,0,1,1,0],
    [0,0,1,1,1,1,0,0],
    [0,1,1,0,1,0,0,1]
    ]);

// Test GrayMap(H) function
OutputMap := GrayMap(H);
assert Domain(OutputMap) eq Zps;
assert Codomain(OutputMap) eq V;
assert OutputMap(Zps!0) eq V!0;
assert [OutputMap(a) : a in [0..p^(s-1)-1]] eq Rows(H);
for a in [0..p^(s-1)-1] do
    for lambda in [0..p-1] do
        lambdaVector := V![(p-lambda)^^p^(s-1)];
        assert OutputMap(Zps!a) eq OutputMap(Zps!(a+lambda*p^(s-1))) + lambdaVector;
    end for;
end for;

// Test GrayMap(C, H) function
OutputMapCode := GrayMap(C, H);
assert Domain(OutputMapCode) eq C;
assert Codomain(OutputMapCode) eq Vn;
assert OutputMapCode(C!0) eq Vn!0;
randomWord := Random(C);
assert OutputMapCode(randomWord) eq Vn!&cat[Eltseq(OutputMap(u)) : 
                                            u in Eltseq(randomWord)];

// Test GrayMapImage(C, H) function
OutputMapCodeImage := GrayMapImage(C, H);
expectedOutputMapCodeImage := [OutputMapCode(c) : c in C];
assert #OutputMapCodeImage eq #expectedOutputMapCodeImage;
assert OutputMapCodeImage eq expectedOutputMapCodeImage;

/****************************************************************/
print "test 8: Linear code over Z64 of type (4; 1,0,1,1,0,1)";

p := 2;
s := 6;
n := 4;

nZp := n*p^(s-1);
Zps := Integers(p^s);
V := VectorSpace(GF(p), p^(s-1));
Vn := VectorSpace(GF(p), nZp);

G := [ [1,0,0,0],
       [0,4,0,0],
       [0,0,8,0],
       [0,0,0,32] ];
C := LinearCode(Matrix(Zps, G));

H := Matrix(GF(p), p^(s-1) , p^(s-1), [
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
    [0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1],
    [0,0,0,0,0,0,1,1,0,0,1,1,1,1,1,1,0,0,1,1,1,1,1,1,0,0,0,0,0,0,1,1],
    [0,0,1,1,1,1,1,1,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1,1,1],
    [0,0,0,1,0,1,0,0,1,1,1,1,0,1,0,1,1,1,1,1,0,1,0,1,0,0,0,1,0,1,0,0],
    [0,0,1,0,1,0,0,0,1,1,1,1,1,0,1,0,1,1,1,1,1,0,1,0,0,0,1,0,1,0,0,0],
    [0,0,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1],
    [0,0,1,1,1,1,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1],
    [0,0,1,0,0,1,1,1,1,1,0,0,0,1,1,0,1,1,0,0,0,1,1,0,0,0,1,0,0,1,1,1],
    [0,0,0,1,1,0,1,1,1,1,0,0,1,0,0,1,1,1,0,0,1,0,0,1,0,0,0,1,1,0,1,1],
    [0,0,1,0,0,1,1,1,1,1,0,0,1,0,0,1,0,0,1,1,0,1,1,0,1,1,0,1,1,0,0,0],
    [0,0,0,1,1,0,1,1,1,1,0,0,0,1,1,0,0,0,1,1,1,0,0,1,1,1,1,0,0,1,0,0],
    [0,0,1,0,1,0,1,1,0,0,1,1,0,1,0,1,1,1,0,0,1,0,1,0,1,1,0,1,0,1,0,0],
    [0,0,0,1,0,1,1,1,0,0,1,1,1,0,1,0,1,1,0,0,0,1,0,1,1,1,1,0,1,0,0,0],
    [0,1,0,0,1,1,0,1,0,1,0,1,0,0,1,1,0,1,0,1,0,0,1,1,0,1,0,0,1,1,0,1],
    [0,1,1,1,0,0,0,1,0,1,0,1,1,1,0,0,0,1,0,1,1,1,0,0,0,1,1,1,0,0,0,1],
    [0,1,1,1,0,0,0,1,1,0,1,0,0,0,1,1,1,0,1,0,0,0,1,1,0,1,1,1,0,0,0,1],
    [0,1,1,1,0,0,1,0,1,0,0,1,1,1,0,0,0,1,1,0,0,0,1,1,1,0,0,0,1,1,0,1],
    [0,1,0,1,1,0,0,1,0,1,0,1,1,0,0,1,1,0,1,0,0,1,1,0,1,0,1,0,0,1,1,0],
    [0,1,1,0,0,1,0,1,0,1,0,1,0,1,1,0,1,0,1,0,1,0,0,1,1,0,0,1,1,0,1,0],
    [0,1,0,1,0,1,1,0,0,1,1,0,1,0,1,0,0,1,1,0,1,0,1,0,0,1,0,1,0,1,1,0],
    [0,1,1,0,1,0,1,0,0,1,1,0,0,1,0,1,0,1,1,0,0,1,0,1,0,1,1,0,1,0,1,0],
    [0,1,0,1,0,1,1,0,0,1,1,0,0,1,0,1,1,0,0,1,1,0,1,0,1,0,1,0,1,0,0,1],
    [0,1,1,0,1,0,1,0,0,1,1,0,1,0,1,0,1,0,0,1,0,1,0,1,1,0,0,1,0,1,0,1],
    [0,1,0,1,1,0,1,0,1,0,0,1,0,1,1,0,1,0,0,1,0,1,1,0,0,1,0,1,1,0,1,0],
    [0,1,1,0,0,1,1,0,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,0,1,1,0,0,1,1,0],
    [0,1,0,0,1,1,0,1,1,0,1,0,1,1,0,0,1,0,1,0,1,1,0,0,0,1,0,0,1,1,0,1],
    [0,1,0,0,1,1,1,0,1,0,0,1,0,0,1,1,0,1,1,0,1,1,0,0,1,0,1,1,0,0,0,1],
    [0,0,1,1,1,1,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,1,1,1,1,0,0],
    [0,1,1,1,0,0,0,1,1,0,1,0,0,0,1,1,0,1,0,1,1,1,0,0,1,0,0,0,1,1,1,0],
    [0,1,0,0,1,1,0,1,1,0,1,0,1,1,0,0,0,1,0,1,0,0,1,1,1,0,1,1,0,0,1,0]
    ]);

// Test GrayMap(H) function
OutputMap := GrayMap(H);
assert Domain(OutputMap) eq Zps;
assert Codomain(OutputMap) eq V;
assert OutputMap(Zps!0) eq V!0;
assert [OutputMap(a) : a in [0..p^(s-1)-1]] eq Rows(H);
for a in [0..p^(s-1)-1] do
    for lambda in [0..p-1] do
        lambdaVector := V![(p-lambda)^^p^(s-1)];
        assert OutputMap(Zps!a) eq OutputMap(Zps!(a+lambda*p^(s-1))) + lambdaVector;
    end for;
end for;

// Test GrayMap(C, H) function
OutputMapCode := GrayMap(C, H);
assert Domain(OutputMapCode) eq C;
assert Codomain(OutputMapCode) eq Vn;
assert OutputMapCode(C!0) eq Vn!0;
randomWord := Random(C);
assert OutputMapCode(randomWord) eq Vn!&cat[Eltseq(OutputMap(u)) : 
                                            u in Eltseq(randomWord)];

// Test GrayMapImage(C, H) function
OutputMapCodeImage := GrayMapImage(C, H);
expectedOutputMapCodeImage := [OutputMapCode(c) : c in C];
assert #OutputMapCodeImage eq #expectedOutputMapCodeImage;
assert OutputMapCodeImage eq expectedOutputMapCodeImage;
