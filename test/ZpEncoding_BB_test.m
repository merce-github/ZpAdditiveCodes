/**************************************************************/
/*                                                            */
/* Project name: Z/p^s-additive codes in MAGMA                */
/* Test file name: ZpInformation_BB_test.m                    */
/*                                                            */
/* Comments: Black-box tests for the functions                */
/*           ZpInformationSpace, ZpInformationSet, ZpType     */  
/*           IsZpInformationSet, Encoding and SystematicEncoding*/
/*           included in the ZpAdditiveCodes_Encoding.m file  */
/*                                                            */                                                
/* Authors: Adrián Torres and Mercè Villanueva                */
/*                                                            */
/* Revision version and last date: v1.0    2021/07/22         */
/*                                 v1.1    2021/08/27         */
/*                                 v1.2    2021/09/29         */
/*                                 v1.3    2023/10/17         */
/*                                 v1.4    2023/11/23         */
/*                                                            */
/**************************************************************/

SetAssertions(true);
// Alarm(30*60);
// SetQuitOnError(true);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////            AUXILIARY FUNCTIONS FOR THE TESTS             ///////////////  
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// This function, given an element u of p^t(Z/p^s) return the corresponding 
// element of the isomorphic group Z/p^(s-t) dividing u by p^t, for t in [0..s-1]. 
// The sequence T=[t1,...,ts] indicates that the first t1 coordinates are over 
// Z/p^s, the next t2 over p(Z/p^(s-1)), ... , and the last ts over p^(s-1)Zp, 
// and they are divided by 1, p, ..., p^(s-1), respectively.
DivideByP := function(u, p, T)
    pos := 1;
    for i in [1..#T] do
        divisor := p^(i-1);
        for j in [1..T[i]] do
            u[pos] := u[pos] div divisor;
            pos +:= 1;
        end for;
    end for;

    return u;
end function;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////                    TEST FOR ENCODING                     ///////////////  
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/****************************************************************/
print "test 0: Zero code over Z3^3 of type (10; 0,0,0)";

p := 3;
s := 3;
n := 10;
np := n*p^(s-1);
Zps := Integers(p^s);
Vp := VectorSpace(GF(p), np);
type := [0,0,0];

C := ZeroCode(Zps, n);

expectedInfSpaceZps := RSpace(Zps, 0);
expectedInfSpaceGFp := VectorSpace(GF(p), 0);

infVectorZps := expectedInfSpaceZps!0;
infVectorGFp := expectedInfSpaceGFp!0;
codewordZps := C![0^^n];
codewordGFp := Vp![0^^np];

// Information sets for C and Cp = Phi(C)
Iset := [];
Ipset := [];

// Test ZpInformationSpace function
// Systematic
outputV, outputVp, outputf, outputfp := ZpInformationSpace(C);
assert outputV eq expectedInfSpaceZps;
assert outputVp eq expectedInfSpaceGFp;
assert outputf(infVectorZps) eq codewordZps;
assert outputfp(infVectorGFp) eq codewordGFp;
assert codewordZps @@ outputf eq infVectorZps;
assert codewordGFp @@ outputfp eq infVectorGFp;
// Not systematic
outputV, outputVp, outputfE, outputfpE := ZpInformationSpace(C : IsSystematicEncoding := false);
assert outputfE(infVectorZps) eq codewordZps;
assert outputfpE(infVectorGFp) eq codewordGFp;
assert codewordZps @@ outputfE eq infVectorZps;
assert codewordGFp @@ outputfpE eq infVectorGFp;

// Test ZpInformationSet function
outputI, outputIp := ZpInformationSet(C);
assert outputI eq Iset;
assert outputIp eq Ipset;

// Check that the information vector over GF(p) coincides with the codeword
// restricted to the information coordinates, that is, it is sysematic.
assert expectedInfSpaceGFp!Eltseq(codewordGFp)[outputIp] eq infVectorGFp;
infVectorGFp_rand := Random(expectedInfSpaceGFp);
codewordGFp_rand := outputfp(infVectorGFp_rand); 
assert expectedInfSpaceGFp!Eltseq(codewordGFp_rand)[outputIp] eq infVectorGFp_rand;

// Check that given a random information vector over Zps, the corresponding
// information vector over GF(p) under the Gray map coincides with the codeword
// over GF(p) restricting to the information coordinates. The codeword is obtained
// by applying the Gray map to the encoding of the information vector over Zps.
infVectorZps_rand := Random(expectedInfSpaceZps);
grayMapInfo := CarletGrayMap(p, type);
infVectorGFp_rand := grayMapInfo(infVectorZps_rand);
codewordZps_rand := outputf(infVectorZps_rand);
grayMapC := CarletGrayMap(C);
codewordGFp_rand := grayMapC(codewordZps_rand);
assert expectedInfSpaceGFp!Eltseq(codewordGFp_rand)[outputIp] eq infVectorGFp_rand;

// Test IsZpInformationSet function
outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, Iset);
assert outputIsInfoSet_C eq true;
assert outputIsInfoSet_Cp eq true;

outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, Ipset);
assert outputIsInfoSet_C eq true;
assert outputIsInfoSet_Cp eq true;

// Test ZpType function
outputType := ZpType(C);
assert outputType eq type;

/****************************************************************/
print "test 1: Linear code over Z2^3 of type (5; 1,1,1) with generator matrix in standard form";

p := 2;
s := 3;
n := 5;
np := n*p^(s-1);
Zps := Integers(p^s);
Vp := VectorSpace(GF(p), np);
type := [1,1,1];

C := LinearCode<Zps, n | [[1,1,0,6,3],
                          [0,2,2,2,6],
                          [0,0,4,0,4]]>;

expectedInfSpaceZps := RSpace(LinearCode(DiagonalMatrix(Zps, [1,2,4])));
expectedInfSpaceGFp := VectorSpace(GF(p), 6);

infVectorZps := expectedInfSpaceZps![3,6,4];
infVectorGFp := expectedInfSpaceGFp![0,1,1,1,0,1];

codewordZps := C![3,7,4,6,5];
codewordGFp := Vp![0,1,1,0, 1,0,0,1, 1,1,1,1, 1,1,0,0, 1,0,1,0];

codewordZpsE := C![3,1,2,0,7]; // For the non-systematic encoding
codewordGFpE := Vp![0,1,1,0, 0,1,0,1, 0,0,1,1, 0,0,0,0, 1,0,0,1];

// Information sets for C and Cp = Phi(C)
Iset := [1,2,3];
Ipset := [1,2,3,5,7,9];

// Sets of coordinates that are not information sets
K1 := [1,2];
K2 := [1,2,4];
K3 := [1,2,3,4,5,6];
K4 := [1,2,7,8,13,14];

// Test ZpInformationSpace function
// Systematic
outputV, outputVp, outputf, outputfp := ZpInformationSpace(C);
assert outputV eq expectedInfSpaceZps;
assert outputVp eq expectedInfSpaceGFp;
assert outputf(infVectorZps) eq codewordZps;
assert outputfp(infVectorGFp) eq codewordGFp;
assert codewordZps @@ outputf eq infVectorZps;
assert codewordGFp @@ outputfp eq infVectorGFp;
// Not systematic
outputV, outputVp, outputfE, outputfpE := ZpInformationSpace(C : IsSystematicEncoding := false);
assert outputfE(infVectorZps) eq codewordZpsE;
assert outputfpE(infVectorGFp) eq codewordGFpE;
assert codewordZpsE @@ outputfE eq infVectorZps;
assert codewordGFpE @@ outputfpE eq infVectorGFp;

// Test Systematic Encoding function
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, infVectorZps);
assert outputCodewordZps eq codewordZps;
assert outputCodewordGFp eq codewordGFp;
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, infVectorGFp);
assert outputCodewordZps eq codewordZps;
assert outputCodewordGFp eq codewordGFp;
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Eltseq(infVectorZps));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Eltseq(infVectorGFp));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Matrix(infVectorZps));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Matrix(infVectorGFp));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];

// Test Encoding function
outputCodewordZps, outputCodewordGFp := Encoding(C, infVectorZps);
assert outputCodewordZps eq codewordZpsE;
assert outputCodewordGFp eq codewordGFpE;
outputCodewordZps, outputCodewordGFp := Encoding(C, infVectorGFp);
assert outputCodewordZps eq codewordZpsE;
assert outputCodewordGFp eq codewordGFpE;
outputCodewordZps, outputCodewordGFp := Encoding(C, Eltseq(infVectorZps));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Eltseq(infVectorGFp));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(infVectorZps));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(infVectorGFp));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];

LZps := Eltseq(infVectorZps) cat Eltseq(infVectorZps);
LGFp := Eltseq(infVectorGFp) cat Eltseq(infVectorGFp);
outputCodewordZps, outputCodewordGFp := Encoding(C, LZps);
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, LGFp);
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(2, #Iset, LZps));
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(2, #Ipset, LGFp));
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];

// Test ZpInformationSet function
outputI, outputIp := ZpInformationSet(C);
assert outputI eq Iset;
assert outputIp eq Ipset;

// Check that the information vector over GF(p) coincides with the codeword
// restricted to the information coordinates, that is, it is sysematic.
assert expectedInfSpaceGFp!Eltseq(codewordGFp)[outputIp] eq infVectorGFp;
infVectorGFp_rand := Random(expectedInfSpaceGFp);
codewordGFp_rand := outputfp(infVectorGFp_rand); 
assert expectedInfSpaceGFp!Eltseq(codewordGFp_rand)[outputIp] eq infVectorGFp_rand;

// Check that given a random information vector over Zps, the corresponding
// information vector over GF(p) under the Gray map coincides with the codeword
// over GF(p) restricting to the information coordinates. The codeword is obtained
// by applying the Gray map to the encoding of the information vector over Zps.
infVectorZps_rand := Random(expectedInfSpaceZps);
grayMapInfo := CarletGrayMap(p, type);
infVectorGFp_rand := grayMapInfo(infVectorZps_rand);
codewordZps_rand := outputf(infVectorZps_rand);
grayMapC := CarletGrayMap(C);
codewordGFp_rand := grayMapC(codewordZps_rand);
assert expectedInfSpaceGFp!Eltseq(codewordGFp_rand)[outputIp] eq infVectorGFp_rand;

// Check that the non-systematic encoding over Z/p^s corresponds to mutliplying
// the corresponding informetion vector over Z/p^s by a generator matrix of C
assert codewordZpsE eq DivideByP(infVectorZps, p, type)*ZpMinRowsGeneratorMatrix(C); 

// Test IsZpInformationSet function
outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, Iset);
assert outputIsInfoSet_C eq true;
assert outputIsInfoSet_Cp eq false;

outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, Ipset);
assert outputIsInfoSet_C eq false;
assert outputIsInfoSet_Cp eq true;

outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, K1);
assert outputIsInfoSet_C eq false;
assert outputIsInfoSet_Cp eq false;

outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, K2);
assert outputIsInfoSet_C eq false;
assert outputIsInfoSet_Cp eq false;

outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, K3);
assert outputIsInfoSet_C eq false;
assert outputIsInfoSet_Cp eq false;

outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, K4);
assert outputIsInfoSet_C eq false;
assert outputIsInfoSet_Cp eq false;

// Test ZpType function
outputType := ZpType(C);
assert outputType eq type;

/****************************************************************/
print "test 2: Linear code over Z3^3 of type (5; 1,1,1) with generator matrix in standard form";

p := 3;
s := 3;
n := 5;
np := n*p^(s-1);
Zps := Integers(p^s);
Vp := VectorSpace(GF(p), np);
type := [1,1,1];

C := LinearCode<Zps, n | [[1,2,5,8,6],
                          [0,3,6,12,21],
                          [0,0,9,9,18]]>;

expectedInfSpaceZps := RSpace(LinearCode(DiagonalMatrix(Zps, [1,3,9])));
expectedInfSpaceGFp := VectorSpace(GF(p), 6);

infVectorZps := expectedInfSpaceZps![21,12,18];
infVectorGFp := expectedInfSpaceGFp![2,2,0,1,2,2];

codewordZps := C![21,12,18,21,24];
codewordGFp := Vp![2,2,2,0,0,0,1,1,1,
                   1,1,1,2,2,2,0,0,0,
                   2,2,2,2,2,2,2,2,2,
                   2,2,2,0,0,0,1,1,1,
                   2,2,2,1,1,1,0,0,0];

codewordZpsE := C![21,0,12,18,3]; // For the non-systematic encoding
codewordGFpE := Vp![2,2,2,0,0,0,1,1,1,
                    0,0,0,0,0,0,0,0,0,
                    1,1,1,2,2,2,0,0,0,
                    2,2,2,2,2,2,2,2,2,
                    0,0,0,1,1,1,2,2,2];

// Information sets for C and Cp = Phi(C)
Iset := [1,2,3];
Ipset := [1,2,4,10,13,19];

// Sets of coordinates that are not information sets
K1 := [1,2];
K2 := [3,4,6];
K3 := [1,2,3,4,5,6];
K4 := [1,2,7,8,13,14];

// Test ZpInformationSpace function
// Systematic
outputV, outputVp, outputf, outputfp := ZpInformationSpace(C);
assert outputV eq expectedInfSpaceZps;
assert outputVp eq expectedInfSpaceGFp;
assert outputf(infVectorZps) eq codewordZps;
assert outputfp(infVectorGFp) eq codewordGFp;
assert codewordZps @@ outputf eq infVectorZps;
assert codewordGFp @@ outputfp eq infVectorGFp;
// Not systematic
outputV, outputVp, outputfE, outputfpE := ZpInformationSpace(C : IsSystematicEncoding := false);
assert outputfE(infVectorZps) eq codewordZpsE;
assert outputfpE(infVectorGFp) eq codewordGFpE;
assert codewordZpsE @@ outputfE eq infVectorZps;
assert codewordGFpE @@ outputfpE eq infVectorGFp;

// Test Systematic Encoding function
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, infVectorZps);
assert outputCodewordZps eq codewordZps;
assert outputCodewordGFp eq codewordGFp;
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, infVectorGFp);
assert outputCodewordZps eq codewordZps;
assert outputCodewordGFp eq codewordGFp;
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Eltseq(infVectorZps));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Eltseq(infVectorGFp));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Matrix(infVectorZps));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Matrix(infVectorGFp));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];

// Test Encoding function
outputCodewordZps, outputCodewordGFp := Encoding(C, infVectorZps);
assert outputCodewordZps eq codewordZpsE;
assert outputCodewordGFp eq codewordGFpE;
outputCodewordZps, outputCodewordGFp := Encoding(C, infVectorGFp);
assert outputCodewordZps eq codewordZpsE;
assert outputCodewordGFp eq codewordGFpE;
outputCodewordZps, outputCodewordGFp := Encoding(C, Eltseq(infVectorZps));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Eltseq(infVectorGFp));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(infVectorZps));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(infVectorGFp));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];

LZps := Eltseq(infVectorZps) cat Eltseq(infVectorZps);
LGFp := Eltseq(infVectorGFp) cat Eltseq(infVectorGFp);
outputCodewordZps, outputCodewordGFp := Encoding(C, LZps);
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, LGFp);
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(2, #Iset, LZps));
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(2, #Ipset, LGFp));
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];

// Test ZpInformationSet function
outputI, outputIp := ZpInformationSet(C);
assert outputI eq Iset;
assert outputIp eq Ipset;

// Check that the information vector over GF(p) coincides with the codeword
// restricted to the information coordinates, that is, it is sysematic.
assert expectedInfSpaceGFp!Eltseq(codewordGFp)[outputIp] eq infVectorGFp;
infVectorGFp_rand := Random(expectedInfSpaceGFp);
codewordGFp_rand := outputfp(infVectorGFp_rand); 
assert expectedInfSpaceGFp!Eltseq(codewordGFp_rand)[outputIp] eq infVectorGFp_rand;

// Check that given a random information vector over Zps, the corresponding
// information vector over GF(p) under the Gray map coincides with the codeword
// over GF(p) restricting to the information coordinates. The codeword is obtained
// by applying the Gray map to the encoding of the information vector over Zps.
infVectorZps_rand := Random(expectedInfSpaceZps);
grayMapInfo := CarletGrayMap(p, type);
infVectorGFp_rand := grayMapInfo(infVectorZps_rand);
codewordZps_rand := outputf(infVectorZps_rand);
grayMapC := CarletGrayMap(C);
codewordGFp_rand := grayMapC(codewordZps_rand);
assert expectedInfSpaceGFp!Eltseq(codewordGFp_rand)[outputIp] eq infVectorGFp_rand;

// Check that the non-systematic encoding over Z/p^s corresponds to mutliplying
// the corresponding informetion vector over Z/p^s by a generator matrix of C
assert codewordZpsE eq DivideByP(infVectorZps, p, type)*ZpMinRowsGeneratorMatrix(C); 

// Test IsZpInformationSet function
outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, Iset);
assert outputIsInfoSet_C eq true;
assert outputIsInfoSet_Cp eq false;

outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, Ipset);
assert outputIsInfoSet_C eq false;
assert outputIsInfoSet_Cp eq true;

outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, K1);
assert outputIsInfoSet_C eq false;
assert outputIsInfoSet_Cp eq false;

outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, K2);
assert outputIsInfoSet_C eq false;
assert outputIsInfoSet_Cp eq false;

outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, K3);
assert outputIsInfoSet_C eq false;
assert outputIsInfoSet_Cp eq false;

outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, K4);
assert outputIsInfoSet_C eq false;
assert outputIsInfoSet_Cp eq false;

// Test ZpType function
outputType := ZpType(C);
assert outputType eq type;

/****************************************************************/
print "test 3: Linear code over Z3^3 of type (6; 1,2,2) with generator matrix in standard form";

p := 3;
s := 3;
n := 6;
np := n*p^(s-1);
Zps := Integers(p^s);
Vp := VectorSpace(GF(p), np);
type := [1,2,2];

C := LinearCode<Zps, n | [[1,0,1,4,2,8],
                          [0,3,0,0,6,3],
                          [0,0,3,3,0,24],
                          [0,0,0,9,0,0],
                          [0,0,0,0,9,18]]>;

expectedInfSpaceZps := RSpace(LinearCode(DiagonalMatrix(Zps, [1,3,3,9,9])));
expectedInfSpaceGFp := VectorSpace(GF(p), 9);

infVectorZps := expectedInfSpaceZps![13,24,9,0,18];
infVectorGFp := expectedInfSpaceGFp![1,2,2,2,1,1,1,0,2];

codewordZps := C![13,24,10,4,20,23];
codewordGFp := Vp![1,2,0,2,0,1,0,1,2,
                   2,2,2,1,1,1,0,0,0,
                   1,2,0,1,2,0,1,2,0,
                   0,1,2,1,2,0,2,0,1,
                   2,1,0,2,1,0,2,1,0,
                   2,1,0,0,2,1,1,0,2];

codewordZpsE := C![13,24,22,7,11,20]; // For the non-systematic encoding
codewordGFpE := Vp![1,2,0,2,0,1,0,1,2,
                    2,2,2,1,1,1,0,0,0,
                    2,0,1,0,1,2,1,2,0,
                    0,1,2,2,0,1,1,2,0,
                    1,0,2,1,0,2,1,0,2,
                    2,1,0,2,1,0,2,1,0];

// Information sets for C and Cp = Phi(C)
Iset := [1,2,3,4,5];
Ipset := [1,2,4,10,13,19,22,28,37];

// Sets of coordinates that are not information sets
K1 := [1,2];
K2 := [2,3,4,5,7];
K3 := [1,2,3,4,5,6,7,8,9];
K4 := [1,2,7,8,13,14];

// Test ZpInformationSpace function
// Systematic
outputV, outputVp, outputf, outputfp := ZpInformationSpace(C);
assert outputV eq expectedInfSpaceZps;
assert outputVp eq expectedInfSpaceGFp;
assert outputf(infVectorZps) eq codewordZps;
assert outputfp(infVectorGFp) eq codewordGFp;
assert codewordZps @@ outputf eq infVectorZps;
assert codewordGFp @@ outputfp eq infVectorGFp;
// Not systematic
outputV, outputVp, outputfE, outputfpE := ZpInformationSpace(C : IsSystematicEncoding := false);
assert outputfE(infVectorZps) eq codewordZpsE;
assert outputfpE(infVectorGFp) eq codewordGFpE;
assert codewordZpsE @@ outputfE eq infVectorZps;
assert codewordGFpE @@ outputfpE eq infVectorGFp;

// Test Systematic Encoding function
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, infVectorZps);
assert outputCodewordZps eq codewordZps;
assert outputCodewordGFp eq codewordGFp;
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, infVectorGFp);
assert outputCodewordZps eq codewordZps;
assert outputCodewordGFp eq codewordGFp;
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Eltseq(infVectorZps));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Eltseq(infVectorGFp));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Matrix(infVectorZps));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Matrix(infVectorGFp));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];

// Test Encoding function
outputCodewordZps, outputCodewordGFp := Encoding(C, infVectorZps);
assert outputCodewordZps eq codewordZpsE;
assert outputCodewordGFp eq codewordGFpE;
outputCodewordZps, outputCodewordGFp := Encoding(C, infVectorGFp);
assert outputCodewordZps eq codewordZpsE;
assert outputCodewordGFp eq codewordGFpE;
outputCodewordZps, outputCodewordGFp := Encoding(C, Eltseq(infVectorZps));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Eltseq(infVectorGFp));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(infVectorZps));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(infVectorGFp));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];

LZps := Eltseq(infVectorZps) cat Eltseq(infVectorZps);
LGFp := Eltseq(infVectorGFp) cat Eltseq(infVectorGFp);
outputCodewordZps, outputCodewordGFp := Encoding(C, LZps);
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, LGFp);
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(2, #Iset, LZps));
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(2, #Ipset, LGFp));
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];

// Test ZpInformationSet function
outputI, outputIp := ZpInformationSet(C);
assert outputI eq Iset;
assert outputIp eq Ipset;

// Check that the information vector over GF(p) coincides with the codeword
// restricted to the information coordinates, that is, it is sysematic.
assert expectedInfSpaceGFp!Eltseq(codewordGFp)[outputIp] eq infVectorGFp;
infVectorGFp_rand := Random(expectedInfSpaceGFp);
codewordGFp_rand := outputfp(infVectorGFp_rand); 
assert expectedInfSpaceGFp!Eltseq(codewordGFp_rand)[outputIp] eq infVectorGFp_rand;

// Check that given a random information vector over Zps, the corresponding
// information vector over GF(p) under the Gray map coincides with the codeword
// over GF(p) restricting to the information coordinates. The codeword is obtained
// by applying the Gray map to the encoding of the information vector over Zps.
infVectorZps_rand := Random(expectedInfSpaceZps);
grayMapInfo := CarletGrayMap(p, type);
infVectorGFp_rand := grayMapInfo(infVectorZps_rand);
codewordZps_rand := outputf(infVectorZps_rand);
grayMapC := CarletGrayMap(C);
codewordGFp_rand := grayMapC(codewordZps_rand);
assert expectedInfSpaceGFp!Eltseq(codewordGFp_rand)[outputIp] eq infVectorGFp_rand;

// Check that the non-systematic encoding over Z/p^s corresponds to mutliplying
// the corresponding informetion vector over Z/p^s by a generator matrix of C
assert codewordZpsE eq DivideByP(infVectorZps, p, type)*ZpMinRowsGeneratorMatrix(C); 

// Test IsZpInformationSet function
outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, Iset);
assert outputIsInfoSet_C eq true;
assert outputIsInfoSet_Cp eq false;

outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, Ipset);
assert outputIsInfoSet_C eq false;
assert outputIsInfoSet_Cp eq true;

outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, K1);
assert outputIsInfoSet_C eq false;
assert outputIsInfoSet_Cp eq false;

outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, K2);
assert outputIsInfoSet_C eq false;
assert outputIsInfoSet_Cp eq false;

outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, K3);
assert outputIsInfoSet_C eq false;
assert outputIsInfoSet_Cp eq false;

outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, K4);
assert outputIsInfoSet_C eq false;
assert outputIsInfoSet_Cp eq false;

// Test ZpType function
outputType := ZpType(C);
assert outputType eq type;

/****************************************************************/
print "test 4: Simple linear code over Z2^3 of type (3; 1,1,1) with generator matrix not in standard form";

p := 2;
s := 3;
n := 3;
np := n*p^(s-1);
Zps := Integers(p^s);
Vp := VectorSpace(GF(p), np);
type := [1,1,1];

C := LinearCode<Zps, n | [[0,1,0],
                          [0,0,2],
                          [4,0,0]]>;

expectedInfSpaceZps := RSpace(LinearCode(DiagonalMatrix(Zps, [1,2,4])));
expectedInfSpaceGFp := VectorSpace(GF(p), 6);

infVectorZps := expectedInfSpaceZps![1,6,4];
infVectorGFp := expectedInfSpaceGFp![0,1,0,1,0,1];

codewordZps := C![4,1,6];
codewordGFp := Vp![1,1,1,1,
                   0,1,0,1,
                   1,1,0,0];

codewordZpsE := C![4,1,6]; // For the non-systematic encoding
codewordGFpE := Vp![1,1,1,1,
                    0,1,0,1,
                    1,1,0,0];

// Information sets for C and Cp = Phi(C)
Iset := [1,2,3];
Ipset := [1,5,6,7,9,11];

// Sets of coordinates that are not information sets
K1 := [1,2];
K2 := [2,3,4,5,7];
K3 := [1,2,3,4,5,6,7,8,9];
K4 := [1,2,7,8,11,12];

// Test ZpInformationSpace function
// Systematic
outputV, outputVp, outputf, outputfp := ZpInformationSpace(C);
assert outputV eq expectedInfSpaceZps;
assert outputVp eq expectedInfSpaceGFp;
assert outputf(infVectorZps) eq codewordZps;
assert outputfp(infVectorGFp) eq codewordGFp;
assert codewordZps @@ outputf eq infVectorZps;
assert codewordGFp @@ outputfp eq infVectorGFp;
// Not systematic
outputV, outputVp, outputfE, outputfpE := ZpInformationSpace(C : IsSystematicEncoding := false);
assert outputfE(infVectorZps) eq codewordZpsE;
assert outputfpE(infVectorGFp) eq codewordGFpE;
assert codewordZpsE @@ outputfE eq infVectorZps;
assert codewordGFpE @@ outputfpE eq infVectorGFp;

// Test Systematic Encoding function
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, infVectorZps);
assert outputCodewordZps eq codewordZps;
assert outputCodewordGFp eq codewordGFp;
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, infVectorGFp);
assert outputCodewordZps eq codewordZps;
assert outputCodewordGFp eq codewordGFp;
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Eltseq(infVectorZps));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Eltseq(infVectorGFp));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Matrix(infVectorZps));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Matrix(infVectorGFp));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];

// Test Encoding function
outputCodewordZps, outputCodewordGFp := Encoding(C, infVectorZps);
assert outputCodewordZps eq codewordZpsE;
assert outputCodewordGFp eq codewordGFpE;
outputCodewordZps, outputCodewordGFp := Encoding(C, infVectorGFp);
assert outputCodewordZps eq codewordZpsE;
assert outputCodewordGFp eq codewordGFpE;
outputCodewordZps, outputCodewordGFp := Encoding(C, Eltseq(infVectorZps));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Eltseq(infVectorGFp));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(infVectorZps));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(infVectorGFp));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];

LZps := Eltseq(infVectorZps) cat Eltseq(infVectorZps);
LGFp := Eltseq(infVectorGFp) cat Eltseq(infVectorGFp);
outputCodewordZps, outputCodewordGFp := Encoding(C, LZps);
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, LGFp);
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(2, #Iset, LZps));
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(2, #Ipset, LGFp));
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];

// Test ZpInformationSet function
outputI, outputIp := ZpInformationSet(C);
assert Set(outputI) eq Set(Iset);
assert Set(outputIp) eq Set(Ipset);

// Check that the information vector over GF(p) coincides with the codeword
// restricted to the information coordinates, that is, it is sysematic.
assert expectedInfSpaceGFp!Eltseq(codewordGFp)[outputIp] eq infVectorGFp;
infVectorGFp_rand := Random(expectedInfSpaceGFp);
codewordGFp_rand := outputfp(infVectorGFp_rand); 
assert expectedInfSpaceGFp!Eltseq(codewordGFp_rand)[outputIp] eq infVectorGFp_rand;

// Check that given a random information vector over Zps, the corresponding
// information vector over GF(p) under the Gray map coincides with the codeword
// over GF(p) restricting to the information coordinates. The codeword is obtained
// by applying the Gray map to the encoding of the information vector over Zps.
infVectorZps_rand := Random(expectedInfSpaceZps);
grayMapInfo := CarletGrayMap(p, type);
infVectorGFp_rand := grayMapInfo(infVectorZps_rand);
codewordZps_rand := outputf(infVectorZps_rand);
grayMapC := CarletGrayMap(C);
codewordGFp_rand := grayMapC(codewordZps_rand);
assert expectedInfSpaceGFp!Eltseq(codewordGFp_rand)[outputIp] eq infVectorGFp_rand;

// Check that the non-systematic encoding over Z/p^s corresponds to mutliplying
// the corresponding informetion vector over Z/p^s by a generator matrix of C
assert codewordZpsE eq DivideByP(infVectorZps, p, type)*ZpMinRowsGeneratorMatrix(C); 

// Test IsZpInformationSet function
outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, Iset);
assert outputIsInfoSet_C eq true;
assert outputIsInfoSet_Cp eq false;

outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, Ipset);
assert outputIsInfoSet_C eq false;
assert outputIsInfoSet_Cp eq true;

outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, K1);
assert outputIsInfoSet_C eq false;
assert outputIsInfoSet_Cp eq false;

outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, K2);
assert outputIsInfoSet_C eq false;
assert outputIsInfoSet_Cp eq false;

outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, K3);
assert outputIsInfoSet_C eq false;
assert outputIsInfoSet_Cp eq false;

outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, K4);
assert outputIsInfoSet_C eq false;
assert outputIsInfoSet_Cp eq false;

// Test ZpType function
outputType := ZpType(C);
assert outputType eq type;

/****************************************************************/
print "test 5: Linear code over Z5^3 of type (10; 1,2,1) with generator matrix not in standard form";

p := 5;
s := 3;
n := 10;
np := n*p^(s-1);
Zps := Integers(p^s);
Vp := VectorSpace(GF(p), np);
type := [1,2,1];

C := LinearCode<Zps, n | [[1,49,24,14,2,4,58,2,105,44],
                          [0,115,100,5,5,0,40,30,95,85],
                          [0,70,105,5,0,5,25,120,115,55],
                          [0,0,50,25,0,0,75,25,50,75]]>;

expectedInfSpaceZps := RSpace(LinearCode(DiagonalMatrix(Zps, [1,5,5,25])));
expectedInfSpaceGFp := VectorSpace(GF(p), 8);

infVectorZps := expectedInfSpaceZps![89,30,40,0];
infVectorGFp := expectedInfSpaceGFp![3,2,0,1,2,1,4,0];

codewordZps := C![89,31,41,101,3,111,112,23,60,71];
codewordGFp := Vp![3,2,1,0,4,0,4,3,2,1,2,1,0,4,3,4,3,2,1,0,1,0,4,3,2,
                   1,2,3,4,0,2,3,4,0,1,3,4,0,1,2,4,0,1,2,3,0,1,2,3,4,
                   1,2,3,4,0,4,0,1,2,3,2,3,4,0,1,0,1,2,3,4,3,4,0,1,2,
                   4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,
                   0,3,1,4,2,0,3,1,4,2,0,3,1,4,2,0,3,1,4,2,0,3,1,4,2,
                   4,0,1,2,3,1,2,3,4,0,3,4,0,1,2,0,1,2,3,4,2,3,4,0,1,
                   4,1,3,0,2,1,3,0,2,4,3,0,2,4,1,0,2,4,1,3,2,4,1,3,0,
                   0,3,1,4,2,4,2,0,3,1,3,1,4,2,0,2,0,3,1,4,1,4,2,0,3,
                   2,2,2,2,2,4,4,4,4,4,1,1,1,1,1,3,3,3,3,3,0,0,0,0,0,
                   2,3,4,0,1,1,2,3,4,0,0,1,2,3,4,4,0,1,2,3,3,4,0,1,2];

codewordZpsE := C![89,11,21,76,98,16,72,38,5,66];
codewordGFpE := Vp![3,2,1,0,4,0,4,3,2,1,2,1,0,4,3,4,3,2,1,0,1,0,4,3,2,
                    0,1,2,3,4,2,3,4,0,1,4,0,1,2,3,1,2,3,4,0,3,4,0,1,2,
                    0,1,2,3,4,4,0,1,2,3,3,4,0,1,2,2,3,4,0,1,1,2,3,4,0,
                    3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,
                    3,1,4,2,0,2,0,3,1,4,1,4,2,0,3,0,3,1,4,2,4,2,0,3,1,
                    0,1,2,3,4,3,4,0,1,2,1,2,3,4,0,4,0,1,2,3,2,3,4,0,1,
                    2,4,1,3,0,1,3,0,2,4,0,2,4,1,3,4,1,3,0,2,3,0,2,4,1,
                    1,4,2,0,3,3,1,4,2,0,0,3,1,4,2,2,0,3,1,4,4,2,0,3,1,
                    0,0,0,0,0,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,
                    2,3,4,0,1,0,1,2,3,4,3,4,0,1,2,1,2,3,4,0,4,0,1,2,3];

// Information sets for C and Cp = Phi(C)
Iset := [1,2,3,5];
Ipset := [1,2,6,26,31,51,56,101];

// Sets of coordinates that are not information sets
K1 := [1,2];
K2 := [2,3,4,5,7];
K3 := [1,2,3,4,5,6,7,8,9];
K4 := [1,2,7,8,11,12,100];

// Test ZpInformationSpace function
// Systematic
outputV, outputVp, outputf, outputfp := ZpInformationSpace(C);
assert outputV eq expectedInfSpaceZps;
assert outputVp eq expectedInfSpaceGFp;
assert outputf(infVectorZps) eq codewordZps;
assert outputfp(infVectorGFp) eq codewordGFp;
assert codewordZps @@ outputf eq infVectorZps;
assert codewordGFp @@ outputfp eq infVectorGFp;
// Not systematic
outputV, outputVp, outputfE, outputfpE := ZpInformationSpace(C : IsSystematicEncoding := false);
assert outputfE(infVectorZps) eq codewordZpsE;
assert outputfpE(infVectorGFp) eq codewordGFpE;
assert codewordZpsE @@ outputfE eq infVectorZps;
assert codewordGFpE @@ outputfpE eq infVectorGFp;

// Test Systematic Encoding function
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, infVectorZps);
assert outputCodewordZps eq codewordZps;
assert outputCodewordGFp eq codewordGFp;
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, infVectorGFp);
assert outputCodewordZps eq codewordZps;
assert outputCodewordGFp eq codewordGFp;
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Eltseq(infVectorZps));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Eltseq(infVectorGFp));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Matrix(infVectorZps));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Matrix(infVectorGFp));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];

// Test Encoding function
outputCodewordZps, outputCodewordGFp := Encoding(C, infVectorZps);
assert outputCodewordZps eq codewordZpsE;
assert outputCodewordGFp eq codewordGFpE;
outputCodewordZps, outputCodewordGFp := Encoding(C, infVectorGFp);
assert outputCodewordZps eq codewordZpsE;
assert outputCodewordGFp eq codewordGFpE;
outputCodewordZps, outputCodewordGFp := Encoding(C, Eltseq(infVectorZps));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Eltseq(infVectorGFp));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(infVectorZps));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(infVectorGFp));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];

LZps := Eltseq(infVectorZps) cat Eltseq(infVectorZps);
LGFp := Eltseq(infVectorGFp) cat Eltseq(infVectorGFp);
outputCodewordZps, outputCodewordGFp := Encoding(C, LZps);
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, LGFp);
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(2, #Iset, LZps));
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(2, #Ipset, LGFp));
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];

// Test ZpInformationSet function
outputI, outputIp := ZpInformationSet(C);
assert Set(outputI) eq Set(Iset);
assert Set(outputIp) eq Set(Ipset);

// Check that the information vector over GF(p) coincides with the codeword
// restricted to the information coordinates, that is, it is sysematic.
assert expectedInfSpaceGFp!Eltseq(codewordGFp)[outputIp] eq infVectorGFp;
infVectorGFp_rand := Random(expectedInfSpaceGFp);
codewordGFp_rand := outputfp(infVectorGFp_rand); 
assert expectedInfSpaceGFp!Eltseq(codewordGFp_rand)[outputIp] eq infVectorGFp_rand;

// Check that given a random information vector over Zps, the corresponding
// information vector over GF(p) under the Gray map coincides with the codeword
// over GF(p) restricting to the information coordinates. The codeword is obtained
// by applying the Gray map to the encoding of the information vector over Zps.
infVectorZps_rand := Random(expectedInfSpaceZps);
grayMapInfo := CarletGrayMap(p, type);
infVectorGFp_rand := grayMapInfo(infVectorZps_rand);
codewordZps_rand := outputf(infVectorZps_rand);
grayMapC := CarletGrayMap(C);
codewordGFp_rand := grayMapC(codewordZps_rand);
assert expectedInfSpaceGFp!Eltseq(codewordGFp_rand)[outputIp] eq infVectorGFp_rand;

// Check that the non-systematic encoding over Z/p^s corresponds to mutliplying
// the corresponding informetion vector over Z/p^s by a generator matrix of C
assert codewordZpsE eq DivideByP(infVectorZps, p, type)*ZpMinRowsGeneratorMatrix(C); 

// Test IsZpInformationSet function
outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, Iset);
assert outputIsInfoSet_C eq true;
assert outputIsInfoSet_Cp eq false;

outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, Ipset);
assert outputIsInfoSet_C eq false;
assert outputIsInfoSet_Cp eq true;

outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, K1);
assert outputIsInfoSet_C eq false;
assert outputIsInfoSet_Cp eq false;

outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, K2);
assert outputIsInfoSet_C eq false;
assert outputIsInfoSet_Cp eq false;

outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, K3);
assert outputIsInfoSet_C eq false;
assert outputIsInfoSet_Cp eq false;

outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, K4);
assert outputIsInfoSet_C eq false;
assert outputIsInfoSet_Cp eq false;

// Test ZpType function
outputType := ZpType(C);
assert outputType eq type;

/****************************************************************/
print "test 6: Hadamard code over Z3^3 of type (27; 1,1,1) with generator matrix not in standard form";

p := 3;
s := 3;
n := 27;
np := n*p^(s-1);
Zps := Integers(p^s);
Vp := VectorSpace(GF(p), np);
type := [1,1,1];

C := ZpHadamardCode(p,type);

expectedInfSpaceZps := RSpace(LinearCode(DiagonalMatrix(Zps, [1,3,9])));
expectedInfSpaceGFp := VectorSpace(GF(p), 6);

infVectorZps := expectedInfSpaceZps![3,12,0];
infVectorGFp := expectedInfSpaceGFp![0,0,1,1,2,0];

codewordZps := C![3,12,21,3,12,21,3,12,21,3,12,21,3,12,21,3,12,21,3,12,21,3,12,21,3,12,21];
codewordGFp := Vp![0,0,0,1,1,1,2,2,2, 1,1,1,2,2,2,0,0,0, 2,2,2,0,0,0,1,1,1, 
                   0,0,0,1,1,1,2,2,2, 1,1,1,2,2,2,0,0,0, 2,2,2,0,0,0,1,1,1, 
                   0,0,0,1,1,1,2,2,2, 1,1,1,2,2,2,0,0,0, 2,2,2,0,0,0,1,1,1, 
                   0,0,0,1,1,1,2,2,2, 1,1,1,2,2,2,0,0,0, 2,2,2,0,0,0,1,1,1,
                   0,0,0,1,1,1,2,2,2, 1,1,1,2,2,2,0,0,0, 2,2,2,0,0,0,1,1,1, 
                   0,0,0,1,1,1,2,2,2, 1,1,1,2,2,2,0,0,0, 2,2,2,0,0,0,1,1,1, 
                   0,0,0,1,1,1,2,2,2, 1,1,1,2,2,2,0,0,0, 2,2,2,0,0,0,1,1,1, 
                   0,0,0,1,1,1,2,2,2, 1,1,1,2,2,2,0,0,0, 2,2,2,0,0,0,1,1,1,
                   0,0,0,1,1,1,2,2,2, 1,1,1,2,2,2,0,0,0, 2,2,2,0,0,0,1,1,1];

codewordZpsE := C![3,15,0,12,24,9,21,6,18,3,15,0,12,24,9,21,6,18,3,15,0,12,24,9,21,6,18];
codewordGFpE := Vp![0,0,0,1,1,1,2,2,2, 1,1,1,0,0,0,2,2,2, 0,0,0,0,0,0,0,0,0,
                    1,1,1,2,2,2,0,0,0, 2,2,2,1,1,1,0,0,0, 1,1,1,1,1,1,1,1,1,
                    2,2,2,0,0,0,1,1,1, 0,0,0,2,2,2,1,1,1, 2,2,2,2,2,2,2,2,2,
                    0,0,0,1,1,1,2,2,2, 1,1,1,0,0,0,2,2,2, 0,0,0,0,0,0,0,0,0,
                    1,1,1,2,2,2,0,0,0, 2,2,2,1,1,1,0,0,0, 1,1,1,1,1,1,1,1,1,
                    2,2,2,0,0,0,1,1,1, 0,0,0,2,2,2,1,1,1, 2,2,2,2,2,2,2,2,2,
                    0,0,0,1,1,1,2,2,2, 1,1,1,0,0,0,2,2,2, 0,0,0,0,0,0,0,0,0,
                    1,1,1,2,2,2,0,0,0, 2,2,2,1,1,1,0,0,0, 1,1,1,1,1,1,1,1,1,
                    2,2,2,0,0,0,1,1,1, 0,0,0,2,2,2,1,1,1, 2,2,2,2,2,2,2,2,2];

// Information sets for C and Cp = Phi(C)
Iset := [1,2,10];
Ipset := [1,2,4,10,13,82];

// Sets of coordinates that are not information sets
K1 := [1,2];
K2 := [2,3,4,5,7];
K3 := [1,2,3,4,5,6,7,8,9];
K4 := [1,2,7,8,11,12,100];

// Test ZpInformationSpace function
// Systematic
outputV, outputVp, outputf, outputfp := ZpInformationSpace(C);
assert outputV eq expectedInfSpaceZps;
assert outputVp eq expectedInfSpaceGFp;
assert outputf(infVectorZps) eq codewordZps;
assert outputfp(infVectorGFp) eq codewordGFp;
assert codewordZps @@ outputf eq infVectorZps;
assert codewordGFp @@ outputfp eq infVectorGFp;
// Not systematic
outputV, outputVp, outputfE, outputfpE := ZpInformationSpace(C : IsSystematicEncoding := false);
assert outputfE(infVectorZps) eq codewordZpsE;
assert outputfpE(infVectorGFp) eq codewordGFpE;
assert codewordZpsE @@ outputfE eq infVectorZps;
assert codewordGFpE @@ outputfpE eq infVectorGFp;

// Test Systematic Encoding function
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, infVectorZps);
assert outputCodewordZps eq codewordZps;
assert outputCodewordGFp eq codewordGFp;
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, infVectorGFp);
assert outputCodewordZps eq codewordZps;
assert outputCodewordGFp eq codewordGFp;
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Eltseq(infVectorZps));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Eltseq(infVectorGFp));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Matrix(infVectorZps));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Matrix(infVectorGFp));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];

// Test Encoding function
outputCodewordZps, outputCodewordGFp := Encoding(C, infVectorZps);
assert outputCodewordZps eq codewordZpsE;
assert outputCodewordGFp eq codewordGFpE;
outputCodewordZps, outputCodewordGFp := Encoding(C, infVectorGFp);
assert outputCodewordZps eq codewordZpsE;
assert outputCodewordGFp eq codewordGFpE;
outputCodewordZps, outputCodewordGFp := Encoding(C, Eltseq(infVectorZps));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Eltseq(infVectorGFp));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(infVectorZps));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(infVectorGFp));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];

LZps := Eltseq(infVectorZps) cat Eltseq(infVectorZps);
LGFp := Eltseq(infVectorGFp) cat Eltseq(infVectorGFp);
outputCodewordZps, outputCodewordGFp := Encoding(C, LZps);
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, LGFp);
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(2, #Iset, LZps));
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(2, #Ipset, LGFp));
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];

// Test ZpInformationSet function
outputI, outputIp := ZpInformationSet(C);
assert Set(outputI) eq Set(Iset);
assert Set(outputIp) eq Set(Ipset);

// Check that the information vector over GF(p) coincides with the codeword
// restricted to the information coordinates, that is, it is sysematic.
assert expectedInfSpaceGFp!Eltseq(codewordGFp)[outputIp] eq infVectorGFp;
infVectorGFp_rand := Random(expectedInfSpaceGFp);
codewordGFp_rand := outputfp(infVectorGFp_rand); 
assert expectedInfSpaceGFp!Eltseq(codewordGFp_rand)[outputIp] eq infVectorGFp_rand;

// Check that given a random information vector over Zps, the corresponding
// information vector over GF(p) under the Gray map coincides with the codeword
// over GF(p) restricting to the information coordinates. The codeword is obtained
// by applying the Gray map to the encoding of the information vector over Zps.
infVectorZps_rand := Random(expectedInfSpaceZps);
grayMapInfo := CarletGrayMap(p, type);
infVectorGFp_rand := grayMapInfo(infVectorZps_rand);
codewordZps_rand := outputf(infVectorZps_rand);
grayMapC := CarletGrayMap(C);
codewordGFp_rand := grayMapC(codewordZps_rand);
assert expectedInfSpaceGFp!Eltseq(codewordGFp_rand)[outputIp] eq infVectorGFp_rand;

// Check that the non-systematic encoding over Z/p^s corresponds to mutliplying
// the corresponding informetion vector over Z/p^s by a generator matrix of C
assert codewordZpsE eq DivideByP(infVectorZps, p, type)*ZpMinRowsGeneratorMatrix(C); 

// Test IsZpInformationSet function
outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, Iset);
assert outputIsInfoSet_C eq true;
assert outputIsInfoSet_Cp eq false;

outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, Ipset);
assert outputIsInfoSet_C eq false;
assert outputIsInfoSet_Cp eq true;

outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, K1);
assert outputIsInfoSet_C eq false;
assert outputIsInfoSet_Cp eq false;

outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, K2);
assert outputIsInfoSet_C eq false;
assert outputIsInfoSet_Cp eq false;

outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, K3);
assert outputIsInfoSet_C eq false;
assert outputIsInfoSet_Cp eq false;

outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, K4);
assert outputIsInfoSet_C eq false;
assert outputIsInfoSet_Cp eq false;

// Test ZpType function
outputType := ZpType(C);
assert outputType eq type;

/****************************************************************/
print "test 7: Linear code over Z3^4 of type (10; 2,2,2,2) with generator matrix not in standard form";

p := 3;
s := 4;
n := 10;
np := n*p^(s-1);
Zps := Integers(p^s);
Vp := VectorSpace(GF(p), np);
type := [2,2,2,2];

C := LinearCode<Zps, n | [[43,0,0,1,0,3,1,22,0,74],
                          [61,1,2,1,18,8,0,10,1,18],
                          [ 0,6,3,0,15,6,0,3,0,57],
                          [36,6,0,3,0,0,0,21,0,0],
                          [72,9,0,0,18,0,0,0,0,18],
                          [72,0,0,0,18,9,0,9,0,54],
                          [ 0,0,0,0,27,0,0,0,0,54],
                          [ 0,0,0,0,0,0,0,27,0,54]]>;

expectedInfSpaceZps := RSpace(LinearCode(DiagonalMatrix(Zps, [1,1,3,3,9,9,27,27])));
expectedInfSpaceGFp := VectorSpace(GF(p), 20);

infVectorZps := expectedInfSpaceZps![73,9,33,9,27,0,27,27];
infVectorGFp := expectedInfSpaceGFp![2,0,2,1, 0,0,0,1, 1,0,1, 0,0,1, 1,1, 0,0, 1, 1];

codewordZps := C![73,9,33,34,9,0,1,46,39,65];
codewordGFp := Vp![2,0,1,2,0,1,2,0,1,1,2,0,1,2,0,1,2,0,0,1,2,0,1,2,0,1,2,
                   0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,
                   1,1,1,0,0,0,2,2,2,1,1,1,0,0,0,2,2,2,1,1,1,0,0,0,2,2,2,
                   1,2,0,0,1,2,2,0,1,1,2,0,0,1,2,2,0,1,1,2,0,0,1,2,2,0,1,
                   0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,
                   1,2,0,1,2,0,1,2,0,0,1,2,0,1,2,0,1,2,2,0,1,2,0,1,2,0,1,
                   1,1,1,2,2,2,0,0,0,2,2,2,0,0,0,1,1,1,0,0,0,1,1,1,2,2,2,
                   2,1,0,2,1,0,2,1,0,0,2,1,0,2,1,0,2,1,1,0,2,1,0,2,1,0,2];

codewordZpsE := C![73,9,51,61,9,72,28,46,3,47];
codewordGFpE := Vp![2,0,1,2,0,1,2,0,1,1,2,0,1,2,0,1,2,0,0,1,2,0,1,2,0,1,2,
                    0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,
                    1,1,1,0,0,0,2,2,2,0,0,0,2,2,2,1,1,1,2,2,2,1,1,1,0,0,0,
                    2,0,1,1,2,0,0,1,2,2,0,1,1,2,0,0,1,2,2,0,1,1,2,0,0,1,2,
                    0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,
                    2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,
                    1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,
                    1,2,0,1,2,0,1,2,0,0,1,2,0,1,2,0,1,2,2,0,1,2,0,1,2,0,1,
                    0,0,0,1,1,1,2,2,2,0,0,0,1,1,1,2,2,2,0,0,0,1,1,1,2,2,2,
                    1,0,2,1,0,2,1,0,2,0,2,1,0,2,1,0,2,1,2,1,0,2,1,0,2,1,0];

// Information sets for C and Cp = Phi(C)
Iset := [1,2,3,5,4,6,8,9];
Ipset := [1,2,4,10,28,29,31,37,55,58,64,109,112,118,82,91,136,145,190,217];

// Sets of coordinates that are not information sets
K1 := [1,2];
K2 := [2,3,4,5,7];
K3 := [1,2,3,4,5,6,7,8,9];
K4 := [1,2,7,8,11,12,100];

// Test ZpInformationSpace function
// Systematic
outputV, outputVp, outputf, outputfp := ZpInformationSpace(C);
assert outputV eq expectedInfSpaceZps;
assert outputVp eq expectedInfSpaceGFp;
assert outputf(infVectorZps) eq codewordZps;
assert outputfp(infVectorGFp) eq codewordGFp;
assert codewordZps @@ outputf eq infVectorZps;
assert codewordGFp @@ outputfp eq infVectorGFp;
// Not systematic
outputV, outputVp, outputfE, outputfpE := ZpInformationSpace(C : IsSystematicEncoding := false);
assert outputfE(infVectorZps) eq codewordZpsE;
assert outputfpE(infVectorGFp) eq codewordGFpE;
assert codewordZpsE @@ outputfE eq infVectorZps;
assert codewordGFpE @@ outputfpE eq infVectorGFp;

// Test Systematic Encoding function
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, infVectorZps);
assert outputCodewordZps eq codewordZps;
assert outputCodewordGFp eq codewordGFp;
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, infVectorGFp);
assert outputCodewordZps eq codewordZps;
assert outputCodewordGFp eq codewordGFp;
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Eltseq(infVectorZps));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Eltseq(infVectorGFp));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Matrix(infVectorZps));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Matrix(infVectorGFp));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];

// Test Encoding function
outputCodewordZps, outputCodewordGFp := Encoding(C, infVectorZps);
assert outputCodewordZps eq codewordZpsE;
assert outputCodewordGFp eq codewordGFpE;
outputCodewordZps, outputCodewordGFp := Encoding(C, infVectorGFp);
assert outputCodewordZps eq codewordZpsE;
assert outputCodewordGFp eq codewordGFpE;
outputCodewordZps, outputCodewordGFp := Encoding(C, Eltseq(infVectorZps));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Eltseq(infVectorGFp));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(infVectorZps));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(infVectorGFp));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];

LZps := Eltseq(infVectorZps) cat Eltseq(infVectorZps);
LGFp := Eltseq(infVectorGFp) cat Eltseq(infVectorGFp);
outputCodewordZps, outputCodewordGFp := Encoding(C, LZps);
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, LGFp);
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(2, #Iset, LZps));
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(2, #Ipset, LGFp));
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];

// Test ZpInformationSet function
outputI, outputIp := ZpInformationSet(C);
assert Set(outputI) eq Set(Iset);
assert Set(outputIp) eq Set(Ipset);

// Check that the information vector over GF(p) coincides with the codeword
// restricted to the information coordinates, that is, it is sysematic.
assert expectedInfSpaceGFp!Eltseq(codewordGFp)[outputIp] eq infVectorGFp;
infVectorGFp_rand := Random(expectedInfSpaceGFp);
codewordGFp_rand := outputfp(infVectorGFp_rand); 
assert expectedInfSpaceGFp!Eltseq(codewordGFp_rand)[outputIp] eq infVectorGFp_rand;

// Check that given a random information vector over Zps, the corresponding
// information vector over GF(p) under the Gray map coincides with the codeword
// over GF(p) restricting to the information coordinates. The codeword is obtained
// by applying the Gray map to the encoding of the information vector over Zps.
infVectorZps_rand := Random(expectedInfSpaceZps);
grayMapInfo := CarletGrayMap(p, type);
infVectorGFp_rand := grayMapInfo(infVectorZps_rand);
codewordZps_rand := outputf(infVectorZps_rand);
grayMapC := CarletGrayMap(C);
codewordGFp_rand := grayMapC(codewordZps_rand);
assert expectedInfSpaceGFp!Eltseq(codewordGFp_rand)[outputIp] eq infVectorGFp_rand;

// Check that the non-systematic encoding over Z/p^s corresponds to mutliplying
// the corresponding informetion vector over Z/p^s by a generator matrix of C
assert codewordZpsE eq DivideByP(infVectorZps, p, type)*ZpMinRowsGeneratorMatrix(C); 

// Test IsZpInformationSet function // TOO SLOW
// outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, Iset);
// assert outputIsInfoSet_C eq true;
// assert outputIsInfoSet_Cp eq false;

// outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, Ipset);
// assert outputIsInfoSet_C eq false;
// assert outputIsInfoSet_Cp eq true;

// outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, K1);
// assert outputIsInfoSet_C eq false;
// assert outputIsInfoSet_Cp eq false;

// outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, K2);
// assert outputIsInfoSet_C eq false;
// assert outputIsInfoSet_Cp eq false;

// outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, K3);
// assert outputIsInfoSet_C eq false;
// assert outputIsInfoSet_Cp eq false;

// outputIsInfoSet_C, outputIsInfoSet_Cp := IsZpInformationSet(C, K4);
// assert outputIsInfoSet_C eq false;
// assert outputIsInfoSet_Cp eq false;

// Test ZpType function
outputType := ZpType(C);
assert outputType eq type;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////                    DIRECT ENCODING                       ///////////////  
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/****************************************************************/
print "test 8: Linear code over Z3^3 of type (5; 1,1,1)";

p := 3;
s := 3;
n := 5;
np := n*p^(s-1);
Zps := Integers(p^s);
Vp := VectorSpace(GF(p), np);
type := [1,1,1];

C := LinearCode<Zps, n | [[1,2,3,23,17],
                          [0,3,3,6,18],
                          [0,0,9,0,9]]>;

expectedInfSpaceZps := RSpace(LinearCode(DiagonalMatrix(Zps, [1,3,9])));
expectedInfSpaceGFp := VectorSpace(GF(p), 6);

infVectorZps := expectedInfSpaceZps![5,9,18];
infVectorGFp := expectedInfSpaceGFp![0,2,1, 1,1, 2];
codewordZps := C![5,10,24,7,13];
codewordGFp := Vp![0,2,1,1,0,2,2,1,0,1,2,0,1,2,0,1,2,0,2,2,2,1,1,
                    1,0,0,0,0,1,2,2,0,1,1,2,0,1,2,0,2,0,1,0,1,2]; 
codewordZpsE := C![5,19,15,25,22];
codewordGFpE := Vp![0,2,1,1,0,2,2,1,0,2,0,1,2,0,1,2,0,1,1,1,1,0,
                    0,0,2,2,2,2,0,1,1,2,0,0,1,2,2,0,1,0,1,2,1,2,0];

// Test ZpInformationSpace function
// Systematic
outputV, outputVp, outputf, outputfp := ZpInformationSpace(C);
assert outputV eq expectedInfSpaceZps;
assert outputVp eq expectedInfSpaceGFp;
assert outputf(infVectorZps) eq codewordZps;
assert outputfp(infVectorGFp) eq codewordGFp;
assert codewordZps @@ outputf eq infVectorZps;
assert codewordGFp @@ outputfp eq infVectorGFp;
// Not systematic
outputV, outputVp, outputfE, outputfpE := ZpInformationSpace(C : IsSystematicEncoding := false);
assert outputfE(infVectorZps) eq codewordZpsE;
assert outputfpE(infVectorGFp) eq codewordGFpE;
assert codewordZpsE @@ outputfE eq infVectorZps;
assert codewordGFpE @@ outputfpE eq infVectorGFp;

// Test Systematic Encoding function
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, infVectorZps);
assert outputCodewordZps eq codewordZps;
assert outputCodewordGFp eq codewordGFp;
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, infVectorGFp);
assert outputCodewordZps eq codewordZps;
assert outputCodewordGFp eq codewordGFp;
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Eltseq(infVectorZps));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Eltseq(infVectorGFp));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Matrix(infVectorZps));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Matrix(infVectorGFp));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];

// Test Encoding function
outputCodewordZps, outputCodewordGFp := Encoding(C, infVectorZps);
assert outputCodewordZps eq codewordZpsE;
assert outputCodewordGFp eq codewordGFpE;
outputCodewordZps, outputCodewordGFp := Encoding(C, infVectorGFp);
assert outputCodewordZps eq codewordZpsE;
assert outputCodewordGFp eq codewordGFpE;
outputCodewordZps, outputCodewordGFp := Encoding(C, Eltseq(infVectorZps));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Eltseq(infVectorGFp));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(infVectorZps));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(infVectorGFp));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];

LZps := Eltseq(infVectorZps) cat Eltseq(infVectorZps);
LGFp := Eltseq(infVectorGFp) cat Eltseq(infVectorGFp);
outputCodewordZps, outputCodewordGFp := Encoding(C, LZps);
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, LGFp);
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(2, 3, LZps));
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(2, 6, LGFp));
assert outputCodewordZps eq [codewordZpsE, codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE, codewordGFpE];

// Test ZpType function
outputType := ZpType(C);
assert outputType eq type;

/****************************************************************/
print "test 9: Linear code over Z5^2,of type (10; 4,1) with non-minrows generator matrix";

p := 5;
s := 2;
n := 10;
np := n*p^(s-1);
Zps := Integers(p^s);
Vp := VectorSpace(GF(p), np);
type := [4,1];

C := LinearCode<Zps, n | [[ 1,0,0,0,3,1,8,16,7,18],
                          [ 0,1,0,0,2,4,18,22,15,3],
                          [ 0,0,1,3,0,0,16,11,20,0],
                          [ 0,0,0,5,0,4,23,6,20,10],
                          [ 0,0,0,0,5,2,19,3,5,15],
                          [ 0,0,0,0,0,5,10,20,0,0]]>;

expectedInfSpaceZps := RSpace(LinearCode(DiagonalMatrix(Zps, [1,1,1,1,5])));
expectedInfSpaceGFp := VectorSpace(GF(p), 9);

infVectorZps := expectedInfSpaceZps![13,16,16,11,5];
infVectorGFp := expectedInfSpaceGFp![2,0, 3,4, 3,4, 2,3, 1];
codewordZps := C![13,16,16,3,21,11,6,12,21,17];
codewordGFp := Vp![2,0,3,1,4,3,4,0,1,2,3,4,0,1,2,0,3,1,4,2,4,0,1,2,3,
                   2,3,4,0,1,1,2,3,4,0,2,4,1,3,0,4,0,1,2,3,3,0,2,4,1];

// Test ZpInformationSpace function
// Not systematic
outputV, outputVp, outputf, outputfp := ZpInformationSpace(C : IsSystematicEncoding := false);
assert outputf(infVectorZps) eq codewordZps;
assert outputfp(infVectorGFp) eq codewordGFp;
assert codewordZps @@ outputf eq infVectorZps;
assert codewordGFp @@ outputfp eq infVectorGFp;

// Test Encoding function
outputCodewordZps, outputCodewordGFp := Encoding(C, infVectorZps);
assert outputCodewordZps eq codewordZps;
assert outputCodewordGFp eq codewordGFp;
outputCodewordZps, outputCodewordGFp := Encoding(C, infVectorGFp);
assert outputCodewordZps eq codewordZps;
assert outputCodewordGFp eq codewordGFp;
outputCodewordZps, outputCodewordGFp := Encoding(C, Eltseq(infVectorZps));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];
outputCodewordZps, outputCodewordGFp := Encoding(C, Eltseq(infVectorGFp));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(infVectorZps));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(infVectorGFp));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];

LZps := Eltseq(infVectorZps) cat Eltseq(infVectorZps);
LGFp := Eltseq(infVectorGFp) cat Eltseq(infVectorGFp);
outputCodewordZps, outputCodewordGFp := Encoding(C, LZps);
assert outputCodewordZps eq [codewordZps, codewordZps];
assert outputCodewordGFp eq [codewordGFp, codewordGFp];
outputCodewordZps, outputCodewordGFp := Encoding(C, LGFp);
assert outputCodewordZps eq [codewordZps, codewordZps];
assert outputCodewordGFp eq [codewordGFp, codewordGFp];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(2, 5, LZps));
assert outputCodewordZps eq [codewordZps, codewordZps];
assert outputCodewordGFp eq [codewordGFp, codewordGFp];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(2, 9, LGFp));
assert outputCodewordZps eq [codewordZps, codewordZps];
assert outputCodewordGFp eq [codewordGFp, codewordGFp];

// Test ZpType function
outputType := ZpType(C);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////          SEQUENCE/MATRIX SYSTEMATIC/DIRECT ENCODING        /////////////  
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/****************************************************************/
print "test 10: Linear code over Z3^3 of type (10; 2,2,2)";

p := 3;
s := 3;
n := 10;
np := n*p^(s-1);
Zps := Integers(p^s);
Vp := VectorSpace(GF(p), np);
type := [2,2,2];

C := LinearCode<Zps, n | [[ 1, 0, 2, 2, 1, 3,23,22, 4,20],
                          [ 0, 1, 2, 0, 1, 6, 7,15, 17, 4],
                          [ 0, 0, 3, 0, 6, 0,15,18,24,24],
                          [ 0, 0, 0, 3, 0, 6,12,21, 9, 0],
                          [ 0, 0, 0, 0, 9, 0, 0, 9, 0,18],
                          [ 0, 0, 0, 0, 0, 9, 0,18, 9, 0]]>;

expectedInfSpaceZps := RSpace(LinearCode(DiagonalMatrix(Zps, [1,1,3,3,9,9])));
expectedInfSpaceGFp := VectorSpace(GF(p), 12);

infVectorZps := expectedInfSpaceZps![9,23,12,9,9,9];
infVectorGFp := expectedInfSpaceGFp![1,1,1, 2,1,0, 1,2, 1,1, 1, 1];
codewordZps := C![9,23,13,9,11,12,23,12,19,26];
codewordGFp := Vp![1,1,1,1,1,1,1,1,1,2,1,0,0,2,1,1,0,2,1,2,0,2,0,1,0,1,2,1,1,1,
                   1,1,1,1,1,1,1,0,2,1,0,2,1,0,2,1,1,1,2,2,2,0,0,0,2,1,0,0,2,1,
                   1,0,2,1,1,1,2,2,2,0,0,0,2,0,1,2,0,1,2,0,1,2,1,0,1,0,2,0,2,1];
codewordZpsE := C![9,23,22,0,11,3,5,3,19,8];
codewordGFpE := Vp![1,1,1,1,1,1,1,1,1,2,1,0,0,2,1,1,0,2,2,0,1,0,1,2,1,2,0,0,0,0,
                    0,0,0,0,0,0,1,0,2,1,0,2,1,0,2,0,0,0,1,1,1,2,2,2,0,2,1,1,0,2,
                    2,1,0,0,0,0,1,1,1,2,2,2,2,0,1,2,0,1,2,0,1,0,2,1,2,1,0,1,0,2];

infSequenceZps := [9,23,12,9,9,9,19,20,21,3,0,18];
infSequenceGFp := [1,1,1, 2,1,0, 1,2, 1,1, 1, 1, 2,0,2, 2,1,2, 2,0, 0,1, 0, 2];
codewordListZps := [C![9,23,13,9,11,12,23,12,19,26],
                    C![19,20,21,5,6,21,25,19,14,4]];
codewordListGFp := [Vp![1,1,1,1,1,1,1,1,1,2,1,0,0,2,1,1,0,2,1,2,0,2,0,1,0,1,2,1,1,1,
                        1,1,1,1,1,1,1,0,2,1,0,2,1,0,2,1,1,1,2,2,2,0,0,0,2,1,0,0,2,1,
                        1,0,2,1,1,1,2,2,2,0,0,0,2,0,1,2,0,1,2,0,1,2,1,0,1,0,2,0,2,1],
                    Vp![2,0,1,2,0,1,2,0,1,2,1,0,2,1,0,2,1,0,2,2,2,0,0,0,1,1,1,0,2,1,
                        1,0,2,2,1,0,0,0,0,2,2,2,1,1,1,2,2,2,0,0,0,1,1,1,2,0,1,1,2,0,
                        0,1,2,2,0,1,2,0,1,2,0,1,1,0,2,2,1,0,0,2,1,0,1,2,1,2,0,2,0,1]];
codewordListZpsE := [C![9,23,22,0,11,3,5,3,19,8],
                     C![19,20,18,14,0,12,19,10,17,7]];
codewordListGFpE := [Vp![1,1,1,1,1,1,1,1,1,2,1,0,0,2,1,1,0,2,2,0,1,0,1,2,1,2,0,0,0,0,
                         0,0,0,0,0,0,1,0,2,1,0,2,1,0,2,0,0,0,1,1,1,2,2,2,0,2,1,1,0,2,
                         2,1,0,0,0,0,1,1,1,2,2,2,2,0,1,2,0,1,2,0,1,0,2,1,2,1,0,1,0,2],
                     Vp![2,0,1,2,0,1,2,0,1,2,1,0,2,1,0,2,1,0,2,2,2,2,2,2,2,2,2,1,0,2,
                         2,1,0,0,2,1,0,0,0,0,0,0,0,0,0,1,1,1,2,2,2,0,0,0,2,0,1,2,0,1,
                         2,0,1,1,2,0,1,2,0,1,2,0,1,0,2,0,2,1,2,1,0,0,1,2,2,0,1,1,2,0]];

// Test Systematic Encoding function
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, infVectorZps);
assert outputCodewordZps eq codewordZps;
assert outputCodewordGFp eq codewordGFp;
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, infVectorGFp);
assert outputCodewordZps eq codewordZps;
assert outputCodewordGFp eq codewordGFp;
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Eltseq(infVectorZps));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Eltseq(infVectorGFp));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Matrix(infVectorZps));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];
outputCodewordZps, outputCodewordGFp := SystematicEncoding(C, Matrix(infVectorGFp));
assert outputCodewordZps eq [codewordZps];
assert outputCodewordGFp eq [codewordGFp];

// Test Encoding functions
outputCodewordZps, outputCodewordGFp := Encoding(C, infVectorZps);
assert outputCodewordZps eq codewordZpsE;
assert outputCodewordGFp eq codewordGFpE;
outputCodewordZps, outputCodewordGFp := Encoding(C, infVectorGFp);
assert outputCodewordZps eq codewordZpsE;
assert outputCodewordGFp eq codewordGFpE;
outputCodewordZps, outputCodewordGFp := Encoding(C, Eltseq(infVectorZps));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Eltseq(infVectorGFp));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(infVectorZps));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(infVectorGFp));
assert outputCodewordZps eq [codewordZpsE];
assert outputCodewordGFp eq [codewordGFpE];

LZps := Eltseq(RSpace(Zps, 12)!infSequenceZps);
LGFp := Eltseq(VectorSpace(GF(p), 24)!infSequenceGFp);
outputCodewordZps, outputCodewordGFp := Encoding(C, LZps);
assert outputCodewordZps eq codewordListZpsE;
assert outputCodewordGFp eq codewordListGFpE;
outputCodewordZps, outputCodewordGFp := Encoding(C, LGFp);
assert outputCodewordZps eq codewordListZpsE;
assert outputCodewordGFp eq codewordListGFpE;
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(Zps, 6, LZps));
assert outputCodewordZps eq codewordListZpsE;
assert outputCodewordGFp eq codewordListGFpE;
outputCodewordZps, outputCodewordGFp := Encoding(C, Matrix(GF(p), 12, LGFp));
assert outputCodewordZps eq codewordListZpsE;
assert outputCodewordGFp eq codewordListGFpE;