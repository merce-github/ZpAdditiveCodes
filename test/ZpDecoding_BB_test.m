/**************************************************************/
/*                                                            */
/* Project name: Z/p^s-additive codes in MAGMA                */
/* Test file name: ZpDecoding_BB_test.m                       */
/*                                                            */
/* Comments: Black-box tests for the functions                */
/*           IsZpPermutationDecodeSet, ZpPermutationDecode    */
/*           and PDSetZpHadamardCode                          */
/*           included in the ZpAdditiveCodes_Decoding.m file  */
/*                                                            */                                                
/* Authors: Adrián Torres and Mercè Villanueva                */
/*                                                            */
/* Revision version and last date: v1.0    2021/09/09         */
/*                                 v2.0    2023/11/12         */
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

// This function returns the theoretical upper bound for the value of r such
// that an r-PD-set of size r+1 may exist for Hadamard codes over Z/p^s
//import "ZpAdditiveCodes_Decoding.m": MaximumSizePDSetHadamard;
function MaximumSizePDSetHadamard(p, T)
    s := #T;
    sumT := &+T;
    return Floor((p^(&+[(s-i+1)*T[i] : i in [1..s]]-s) - sumT)/sumT); 
end function;

// This function checks that the rows of the inverse matrices 
// of the r-PDset of size r+1 obtained from some constructions
// for Hadamard codes over Z/p^s are all different 
procedure testRowsPDSetZpHadamardCode(p, type, seqMatrices)
    s := #type;
    rows := {};
    for Mp in seqMatrices do
        assert IsInvertible(Mp);
        M := Mp^(-1);
        Include(~rows, M[1]);
        nrow := 1;
        typeH := type; 
        typeH[1] -:= 1; 
        for k in [1..s] do
            for i in [nrow+1..nrow+typeH[k]] do
                Include(~rows, M[1]+p^(k-1)*M[i]);
            end for;
            nrow +:= typeH[k];
        end for;
    end for;
    assert #rows eq (#seqMatrices * &+type); 
end procedure; 

// This function returns the corresponding permutation in Sym(p^(s-1)n) 
// from a permutation in Sym(n) 
PermZpsToPermZp := function(permZps, p, s)	
    n := Degree(Parent(permZps));
    permList := &cat[ [p^(s-1)*(i^permZps)-(p^(s-1)-j) : j in [1..p^(s-1)]]
                                                       : i in [1..n]];
    return Sym(p^(s-1)*n)!permList;	
end function;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////                 PD-SETS FOR HADAMARD CODES               ///////////////  
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/****************************************************************/
// The functions tested in this file do not work for the zero code
// print "test 0: Zero code";

/****************************************************************/
print "test 1: Z/27-linear Hadamard code of type (1; 1,0,0)";

p := 3;
type := [1,0,0];
C := ZpHadamardCode(p, type);
I, Ip := ZpInformationSet(C);

expectedSizePDset := 1;
upperBoundPDset := MaximumSizePDSetHadamard(p, type);
assert upperBoundPDset eq expectedSizePDset-1;
outputI, outputS, outputIp, outputSp, outputQ := PDSetZpHadamardCode(p, type);
assert #outputSp eq expectedSizePDset;
assert #outputS eq expectedSizePDset;
assert outputI eq I; 
assert outputIp eq Ip;
// It does not make sense to check since it is a 0-PD-set
//assert IsZpPermutationDecodeSet(C, outputIp, outputSp, expectedSizePDset-1);
//assert IsZpPermutationDecodeSet(C, outputI, outputS, expectedSizePDset-1);
assert #[p : p in outputS | C^p eq C] eq expectedSizePDset;
assert outputSp eq [PermZpsToPermZp(perm, p, #type) : perm in outputS];
testRowsPDSetZpHadamardCode(p, type, outputQ);

/****************************************************************/
print "test 2: Z/8-linear Hadamard code of type (8; 2,0,0)";
print("This test may take some minutes. Please wait...");

p := 2;
type := [2,0,0];
C := ZpHadamardCode(p, type);
I, Ip := ZpInformationSet(C);

expectedSizePDset := 4;
upperBoundPDset := MaximumSizePDSetHadamard(p, type);
assert upperBoundPDset eq expectedSizePDset-1;
outputI, outputS, outputIp, outputSp, outputQ := PDSetZpHadamardCode(p, type);
assert #outputSp eq expectedSizePDset;
assert #outputS eq expectedSizePDset;
assert outputI eq I;
assert outputIp eq Ip;
assert IsZpPermutationDecodeSet(C, outputIp, outputSp, expectedSizePDset-1);
assert IsZpPermutationDecodeSet(C, outputI, outputS, expectedSizePDset-1);
assert #[p : p in outputS | C^p eq C] eq expectedSizePDset;
assert outputSp eq [PermZpsToPermZp(perm, p, #type) : perm in outputS];
testRowsPDSetZpHadamardCode(p, type, outputQ);

p1 := Sym(8)!(1,3,5,7)(2,4,6,8);
p2 := Sym(8)!(1,5)(2,6)(3,7)(4,8);
p3 := Sym(8)!(1,7,5,3)(2,8,6,4);
S := sub<Sym(8) | p1, p2, p3>;
p1 := Sym(32)!(1,9,17,25)(2,10,18,26)(3,11,19,27)(4,12,20,28)(5,13,21,29)(6,14,22,30)(7,15,23,31)(8,16,24,32);
p2 := Sym(32)!(1,17)(2,18)(3,19)(4,20)(5,21)(6,22)(7,23)(8,24)(9,25)(10,26)(11,27)(12,28)(13,29)(14,30)(15,31)(16,32);
p3 := Sym(32)!(1,25,17,9)(2,26,18,10)(3,27,19,11)(4,28,20,12)(5,29,21,13)(6,30,22,14)(7,31,23,15)(8,32,24,16);
Sp := sub<Sym(32) | p1, p2, p3>;
assert Set(outputS) eq Set(S);
assert Set(outputSp) eq Set(Sp);

p1 := Sym(32)!(5,16)(6,15)(7,30)(8,29)(13,24)(14,23)(21,32)(22,31);
p2 := Sym(32)!(5,7)(6,8)(13,31)(14,32)(15,29)(16,30)(21,23)(22,24);
p3 := Sym(32)!(5,6)(7,8)(13,14)(15,16)(21,22)(23,24)(29,30)(31,32);
p4 := Sym(32)!(1,5)(2,22)(3,15)(4,32)(6,18)(7,10)(8,25)(9,24)(11,30)(12,13)(14,27)(16,20)(17,21)(19,31)(23,26)(28,29);
p5 := Sym(32)!(2,18)(3,19)(6,22)(7,23)(10,26)(11,27)(14,30)(15,31);
Sp := sub<Sym(32) | p1, p2, p3, p4, p5>;
assert IsZpPermutationDecodeSet(C, Ip, [s : s in Sp], 7);

/****************************************************************/
print "test 3: Z/27-linear Hadamard code of type (729; 3,0,0)";

p := 3;
type := [3,0,0];
C := ZpHadamardCode(p, type);
I, Ip := ZpInformationSet(C);

expectedSizePDset := 243;
upperBoundPDset := MaximumSizePDSetHadamard(p, type);
assert upperBoundPDset eq expectedSizePDset-1;
outputI, outputS, outputIp, outputSp, outputQ := PDSetZpHadamardCode(p, type);
assert #outputSp eq expectedSizePDset;
assert #outputS eq expectedSizePDset;
assert outputI eq I;
assert outputIp eq Ip;
// It takes too long...
//assert IsZpPermutationDecodeSet(C, outputIp, outputSp, expectedSizePDset-1);
//assert IsZpPermutationDecodeSet(C, outputI, outputS, expectedSizePDset-1);
// It is checked that the elements of outputS belong to the PAut of C
assert #[p : p in outputS | C^p eq C] eq expectedSizePDset;
assert outputSp eq [PermZpsToPermZp(perm, p, #type) : perm in outputS];
// It is checked that the rows of the inverse matrices are all different 
testRowsPDSetZpHadamardCode(p, type, outputQ);

/****************************************************************/
print "test 4: Z/25-linear Hadamard code of type (625; 3,0)";

p := 5;
type := [3,0];
C := ZpHadamardCode(p, type);
I, Ip := ZpInformationSet(C);

expectedSizePDset := 208;
upperBoundPDset := MaximumSizePDSetHadamard(p, type);
assert upperBoundPDset eq expectedSizePDset-1;
outputI, outputS, outputIp, outputSp, outputQ := PDSetZpHadamardCode(p, type);
assert #outputSp eq expectedSizePDset;
assert #outputS eq expectedSizePDset;
assert outputI eq I;
assert outputIp eq Ip;
// It takes too long...
//assert IsZpPermutationDecodeSet(C, outputIp, outputSp, expectedSizePDset-1);
//assert IsZpPermutationDecodeSet(C, outputI, outputS, expectedSizePDset-1);
// It is checked that the elements of outputS belong to the PAut of C
assert #[p : p in outputS | C^p eq C] eq expectedSizePDset;
assert outputSp eq [PermZpsToPermZp(perm, p, #type) : perm in outputS];
// It is checked that the rows of the inverse matrices are all different 
testRowsPDSetZpHadamardCode(p, type, outputQ);

/****************************************************************/
print "test 5: Z/8-linear Hadamard code of type (4; 1,1,0)";

p := 2;
type := [1,1,0];
C := ZpHadamardCode(p, type);
I, Ip := ZpInformationSet(C);

expectedSizePDset := 2;
upperBoundPDset := MaximumSizePDSetHadamard(p, type);
assert upperBoundPDset eq expectedSizePDset-1;
outputI, outputS, outputIp, outputSp, outputQ := PDSetZpHadamardCode(p, type);
assert #outputSp eq expectedSizePDset;
assert #outputS eq expectedSizePDset;
assert outputI eq I;
assert outputIp eq Ip;
assert IsZpPermutationDecodeSet(C, outputIp, outputSp, expectedSizePDset-1);
assert IsZpPermutationDecodeSet(C, outputI, outputS, expectedSizePDset-1);
assert #[p : p in outputS | C^p eq C] eq expectedSizePDset;
assert outputSp eq [PermZpsToPermZp(perm, p, #type) : perm in outputS];
testRowsPDSetZpHadamardCode(p, type, outputQ);

p1 := Sym(4)!(1,3)(2,4);
S := [Sym(4)!1, p1];
p1 := Sym(16)!(1,9)(2,10)(3,11)(4,12)(5,13)(6,14)(7,15)(8,16);
Sp := [Sym(16)!1, p1];
assert Set(outputS) eq Set(S);
assert Set(outputSp) eq Set(Sp);

/****************************************************************/
print "test 6: Z/8-linear Hadamard code of type (64; 1,3,0)";

p := 2;
type := [1,3,0];
C := ZpHadamardCode(p, type);
I, Ip := ZpInformationSet(C);

expectedSizePDset := 16;
upperBoundPDset := MaximumSizePDSetHadamard(p, type);
assert upperBoundPDset eq expectedSizePDset-1;
outputI, outputS, outputIp, outputSp, outputQ := PDSetZpHadamardCode(p, type);
assert #outputSp eq expectedSizePDset;
assert #outputS eq expectedSizePDset;
assert outputI eq I;
assert outputIp eq Ip;
// It takes too long...
//assert IsZpPermutationDecodeSet(C, outputIp, outputSp, expectedSizePDset-1);
//assert IsZpPermutationDecodeSet(C, outputI, outputS, expectedSizePDset-1);
// It is checked that the elements of outputS belong to the PAut of C
assert #[p : p in outputS | C^p eq C] eq expectedSizePDset;
assert outputSp eq [PermZpsToPermZp(perm, p, #type) : perm in outputS];
// It is checked that the rows of the inverse matrices are all different 
testRowsPDSetZpHadamardCode(p, type, outputQ);

/****************************************************************/
print "test 7: Z/8-linear Hadamard code of type (32; 1,0,5)";
print("This test may take some minutes. Please wait...");

p := 2;
type := [1,0,5];
C := ZpHadamardCode(p, type);
I, Ip := ZpInformationSet(C);

expectedSizePDset := 5;
upperBoundPDset := MaximumSizePDSetHadamard(p, type);
assert upperBoundPDset eq expectedSizePDset-1;
outputI, outputS, outputIp, outputSp, outputQ := PDSetZpHadamardCode(p, type);
assert #outputSp eq expectedSizePDset;
assert #outputS eq expectedSizePDset;
assert outputI eq I;
assert outputIp eq Ip;
assert IsZpPermutationDecodeSet(C, outputIp, outputSp, expectedSizePDset-1);
assert IsZpPermutationDecodeSet(C, outputI, outputS, expectedSizePDset-1);
assert #[p : p in outputS | C^p eq C] eq expectedSizePDset;
assert outputSp eq [PermZpsToPermZp(perm, p, #type) : perm in outputS];
testRowsPDSetZpHadamardCode(p, type, outputQ);

p1 := Sym(32)!(1,19,13,16,20,8,24,31,18,17,28,4,15,25,32,27,9,23,22,10,30,14,5,12,11,2,26,21,3,6)(7,29);
p2 := Sym(32)!(1,19,18,28,17,30,5,11,26,29,3,31,15,2,21,6,13,14,12,32,9,22,4,25,27,23,10,20,24,16,8);
p3 := Sym(32)!(1,19,15,26,7,3,16,24,14,11,21,31,25,9,10,8,13,5,32,23,4,2,29,6,18,30,12,27,22,17,20);
p4 := Sym(32)!(1,19,25,23,17,8,18,20,13,12,9,4,26,3,14,32,22,28,30,11,29,31,2,7,6,15,21,16)(5,27,10,24);
S := [Sym(32)!1, p1, p2, p3, p4];
assert Set(outputS) eq Set(S);

/****************************************************************/
print "test 8: Z/8-linear Hadamard code of type (128; 2,1,2) using recursivity from type (2,0,0)";
print("This test may take some minutes. Please wait...");

p := 2;
type := [2,1,2];
C := ZpHadamardCode(p, type);
I, Ip := ZpInformationSet(C);

expectedSizePDset := 4;
upperBoundPDset := MaximumSizePDSetHadamard(p, type);
assert upperBoundPDset ge expectedSizePDset-1;
outputI, outputS, outputIp, outputSp, outputQ := PDSetZpHadamardCode(p, type);
assert #outputSp eq expectedSizePDset;
assert #outputS eq expectedSizePDset;
assert outputI eq I;
assert outputIp eq Ip;
assert IsZpPermutationDecodeSet(C, outputIp, outputSp, expectedSizePDset-1);
assert IsZpPermutationDecodeSet(C, outputI, outputS, expectedSizePDset-1);
testRowsPDSetZpHadamardCode(p, type, outputQ);

p1 := Sym(128)!(1,7,5,3)(2,8,6,4)(9,15,13,11)(10,16,14,12)(17,23,21,19)(18,24,22,20)(25,31,29,27)
               (26,32,30,28)(33,39,37,35)(34,40,38,36)(41,47,45,43)(42,48,46,44)(49,55,53,51)
               (50,56,54,52)(57,63,61,59)(58,64,62,60)(65,71,69,67)(66,72,70,68)(73,79,77,75)
               (74,80,78,76)(81,87,85,83)(82,88,86,84)(89,95,93,91)(90,96,94,92)(97,103,101,99)
               (98,104,102,100)(105,111,109,107)(106,112,110,108)(113,119,117,115)(114,120,118,116)
               (121,127,125,123)(122,128,126,124);
p2 := Sym(128)!(1,5)(2,6)(3,7)(4,8)(9,13)(10,14)(11,15)(12,16)(17,21)(18,22)(19,23)(20,24)(25,29)
               (26,30)(27,31)(28,32)(33,37)(34,38)(35,39)(36,40)(41,45)(42,46)(43,47)(44,48)(49,53)
               (50,54)(51,55)(52,56)(57,61)(58,62)(59,63)(60,64)(65,69)(66,70)(67,71)(68,72)(73,77)
               (74,78)(75,79)(76,80)(81,85)(82,86)(83,87)(84,88)(89,93)(90,94)(91,95)(92,96)(97,101)
               (98,102)(99,103)(100,104)(105,109)(106,110)(107,111)(108,112)(113,117)(114,118)(115,119)
               (116,120)(121,125)(122,126)(123,127)(124,128);
p3 := Sym(128)!(1,3,5,7)(2,4,6,8)(9,11,13,15)(10,12,14,16)(17,19,21,23)(18,20,22,24)(25,27,29,31)
               (26,28,30,32)(33,35,37,39)(34,36,38,40)(41,43,45,47)(42,44,46,48)(49,51,53,55)
               (50,52,54,56)(57,59,61,63)(58,60,62,64)(65,67,69,71)(66,68,70,72)(73,75,77,79)
               (74,76,78,80)(81,83,85,87)(82,84,86,88)(89,91,93,95)(90,92,94,96)(97,99,101,103)
               (98,100,102,104)(105,107,109,111)(106,108,110,112)(113,115,117,119)(114,116,118,120)
               (121,123,125,127)(122,124,126,128);
S := [Sym(128)!1, p1, p2, p3];
assert Set(outputS) eq Set(S);

/****************************************************************/
print "test 9: Z/16-linear Hadamard code of type (128; 1,1,1,2) using recursivity from type (1,1,0,0)";
print("This test may take some minutes. Please wait...");

p := 2;
type := [1,1,1,2];
C := ZpHadamardCode(p, type);
I, Ip := ZpInformationSet(C);

expectedSizePDset := 4;
upperBoundPDset := MaximumSizePDSetHadamard(p, type);
assert upperBoundPDset ge expectedSizePDset-1;
outputI, outputS, outputIp, outputSp, outputQ := PDSetZpHadamardCode(p, type);
assert #outputSp eq expectedSizePDset;
assert #outputS eq expectedSizePDset;
assert outputI eq I;
assert outputIp eq Ip;
assert IsZpPermutationDecodeSet(C, outputIp, outputSp, expectedSizePDset-1);
assert IsZpPermutationDecodeSet(C, outputI, outputS, expectedSizePDset-1);
testRowsPDSetZpHadamardCode(p, type, outputQ);

// The permutations in the 3-PD-set S are the same as in test 8
//S := [Sym(128)!1, p1, p2, p3];
assert Set(outputS) eq Set(S);

/****************************************************************/
print "test 10: Z/16-linear Hadamard code of type (256; 1,0,2,4) using recursivity from type (1,1,0,0)";

p := 2;
type := [1,0,2,4];
C := ZpHadamardCode(p, type);
I, Ip := ZpInformationSet(C);

expectedSizePDset := 5;
upperBoundPDset := MaximumSizePDSetHadamard(p, type);
assert upperBoundPDset ge expectedSizePDset-1;
outputI, outputS, outputIp, outputSp, outputQ := PDSetZpHadamardCode(p, type);
assert #outputSp eq expectedSizePDset;
assert #outputS eq expectedSizePDset;
assert outputI eq I;
assert outputIp eq Ip;
// It takes too long...
//assert IsZpPermutationDecodeSet(C, outputIp, outputSp, expectedSizePDset-1);
//assert IsZpPermutationDecodeSet(C, outputI, outputS, expectedSizePDset-1);
testRowsPDSetZpHadamardCode(p, type, outputQ);

p1 := Sym(256)!(1,12,13,8,9,4,5,16)(2,15,14,11,10,7,6,3)(17,28,29,24,25,20,21,32)(18,31,30,27,26,23,22,19)
         (33,44,45,40,41,36,37,48)(34,47,46,43,42,39,38,35)(49,60,61,56,57,52,53,64)(50,63,62,59,58,55,54,51)
         (65,76,77,72,73,68,69,80)(66,79,78,75,74,71,70,67)(81,92,93,88,89,84,85,96)(82,95,94,91,90,87,86,83)
         (97,108,109,104,105,100,101,112)(98,111,110,107,106,103,102,99)(113,124,125,120,121,116,117,128)
         (114,127,126,123,122,119,118,115)(129,140,141,136,137,132,133,144)(130,143,142,139,138,135,134,131)
         (145,156,157,152,153,148,149,160)(146,159,158,155,154,151,150,147)(161,172,173,168,169,164,165,176)
         (162,175,174,171,170,167,166,163)(177,188,189,184,185,180,181,192)(178,191,190,187,186,183,182,179)
         (193,204,205,200,201,196,197,208)(194,207,206,203,202,199,198,195)(209,220,221,216,217,212,213,224)
         (210,223,222,219,218,215,214,211)(225,236,237,232,233,228,229,240)(226,239,238,235,234,231,230,227)
         (241,252,253,248,249,244,245,256)(242,255,254,251,250,247,246,243);
p2 := Sym(256)!(1,13,11,7)(2,8,12,14)(3,15,9,5)(4,6,10,16)(17,29,27,23)(18,24,28,30)(19,31,25,21)(20,22,26,32)
         (33,45,43,39)(34,40,44,46)(35,47,41,37)(36,38,42,48)(49,61,59,55)(50,56,60,62)(51,63,57,53)(52,54,58,64)
         (65,77,75,71)(66,72,76,78)(67,79,73,69)(68,70,74,80)(81,93,91,87)(82,88,92,94)(83,95,89,85)(84,86,90,96)
         (97,109,107,103)(98,104,108,110)(99,111,105,101)(100,102,106,112)(113,125,123,119)(114,120,124,126)
         (115,127,121,117)(116,118,122,128)(129,141,139,135)(130,136,140,142)(131,143,137,133)(132,134,138,144)
         (145,157,155,151)(146,152,156,158)(147,159,153,149)(148,150,154,160)(161,173,171,167)(162,168,172,174)
         (163,175,169,165)(164,166,170,176)(177,189,187,183)(178,184,188,190)(179,191,185,181)(180,182,186,192)
         (193,205,203,199)(194,200,204,206)(195,207,201,197)(196,198,202,208)(209,221,219,215)(210,216,220,222)
         (211,223,217,213)(212,214,218,224)(225,237,235,231)(226,232,236,238)(227,239,233,229)(228,230,234,240)
         (241,253,251,247)(242,248,252,254)(243,255,249,245)(244,246,250,256);
p3 := Sym(256)!(1,8,5,10)(2,9,16,13)(3,14,7,4)(6,15,12,11)(17,24,21,26)(18,25,32,29)(19,30,23,20)(22,31,28,27)
         (33,40,37,42)(34,41,48,45)(35,46,39,36)(38,47,44,43)(49,56,53,58)(50,57,64,61)(51,62,55,52)(54,63,60,59)
         (65,72,69,74)(66,73,80,77)(67,78,71,68)(70,79,76,75)(81,88,85,90)(82,89,96,93)(83,94,87,84)(86,95,92,91)
         (97,104,101,106)(98,105,112,109)(99,110,103,100)(102,111,108,107)(113,120,117,122)(114,121,128,125)
         (115,126,119,116)(118,127,124,123)(129,136,133,138)(130,137,144,141)(131,142,135,132)(134,143,140,139)
         (145,152,149,154)(146,153,160,157)(147,158,151,148)(150,159,156,155)(161,168,165,170)(162,169,176,173)
         (163,174,167,164)(166,175,172,171)(177,184,181,186)(178,185,192,189)(179,190,183,180)(182,191,188,187)
         (193,200,197,202)(194,201,208,205)(195,206,199,196)(198,207,204,203)(209,216,213,218)(210,217,224,221)
         (211,222,215,212)(214,223,220,219)(225,232,229,234)(226,233,240,237)(227,238,231,228)(230,239,236,235)
         (241,248,245,250)(242,249,256,253)(243,254,247,244)(246,255,252,251);
p4 := Sym(256)!(1,11)(2,12)(3,9)(4,10)(5,15)(6,16)(7,13)(8,14)(17,27)(18,28)(19,25)(20,26)(21,31)(22,32)(23,29)
         (24,30)(33,43)(34,44)(35,41)(36,42)(37,47)(38,48)(39,45)(40,46)(49,59)(50,60)(51,57)(52,58)(53,63)
         (54,64)(55,61)(56,62)(65,75)(66,76)(67,73)(68,74)(69,79)(70,80)(71,77)(72,78)(81,91)(82,92)(83,89)
         (84,90)(85,95)(86,96)(87,93)(88,94)(97,107)(98,108)(99,105)(100,106)(101,111)(102,112)(103,109)
         (104,110)(113,123)(114,124)(115,121)(116,122)(117,127)(118,128)(119,125)(120,126)(129,139)(130,140)
         (131,137)(132,138)(133,143)(134,144)(135,141)(136,142)(145,155)(146,156)(147,153)(148,154)(149,159)
         (150,160)(151,157)(152,158)(161,171)(162,172)(163,169)(164,170)(165,175)(166,176)(167,173)(168,174)
         (177,187)(178,188)(179,185)(180,186)(181,191)(182,192)(183,189)(184,190)(193,203)(194,204)(195,201)
         (196,202)(197,207)(198,208)(199,205)(200,206)(209,219)(210,220)(211,217)(212,218)(213,223)(214,224)
         (215,221)(216,222)(225,235)(226,236)(227,233)(228,234)(229,239)(230,240)(231,237)(232,238)(241,251)
         (242,252)(243,249)(244,250)(245,255)(246,256)(247,253)(248,254);
S := [Sym(256)!1, p1, p2, p3, p4];
assert Set(outputS) eq Set(S);

/****************************************************************/
print "test 11: Z/16-linear Hadamard code of type (4096; 3,1,0,1) using general construction";
print("This test may take some minutes. Please wait...");

p := 2;
type := [3,1,0,1];
C := ZpHadamardCode(p, type);
I, Ip := ZpInformationSet(C);

expectedSizePDset := 512;
upperBoundPDset := MaximumSizePDSetHadamard(p, type);
assert upperBoundPDset ge expectedSizePDset-1;
outputI, outputS, outputIp, outputSp, outputQ := PDSetZpHadamardCode(p, type);
assert #outputSp eq expectedSizePDset;
assert #outputS eq expectedSizePDset;
assert outputI eq I;
assert outputIp eq Ip;
// It takes too long...
//assert IsZpPermutationDecodeSet(C, outputIp, outputSp, expectedSizePDset-1);
//assert IsZpPermutationDecodeSet(C, outputI, outputS, expectedSizePDset-1);
assert #[p : p in outputS | C^p eq C] eq expectedSizePDset;
assert outputSp eq [PermZpsToPermZp(perm, p, #type) : perm in outputS];
testRowsPDSetZpHadamardCode(p, type, outputQ);

/****************************************************************/
print "test 11b: Z/32-linear Hadamard code of type (4096; 1,2,1,0,1) using general construction from type (3,1,0,1)";
print("This test may take some minutes. Please wait...");

p := 2;
type := [1,2,1,0,1];
C := ZpHadamardCode(p, type);
I, Ip := ZpInformationSet(C);

expectedSizePDset := 512;
upperBoundPDset := MaximumSizePDSetHadamard(p, type);
assert upperBoundPDset ge expectedSizePDset-1;
outputI, outputS, outputIp, outputSp, outputQ := PDSetZpHadamardCode(p, type);
assert #outputSp eq expectedSizePDset;
assert #outputS eq expectedSizePDset;
assert outputI eq I;
assert outputIp eq Ip;
// It takes too long...
//assert IsZpPermutationDecodeSet(C, outputIp, outputSp, expectedSizePDset-1);
//assert IsZpPermutationDecodeSet(C, outputI, outputS, expectedSizePDset-1);
assert #[p : p in outputS | C^p eq C] eq expectedSizePDset;
assert outputSp eq [PermZpsToPermZp(perm, p, #type) : perm in outputS];
testRowsPDSetZpHadamardCode(p, type, outputQ);

/****************************************************************/
print "test 12: Z/8-linear Hadamard code of type (64; 2,1,1) using general construction";

p := 2;
type := [2,1,1];
C := ZpHadamardCode(p, type);
I, Ip := ZpInformationSet(C);

expectedSizePDset := 16;
upperBoundPDset := MaximumSizePDSetHadamard(p, type);
assert upperBoundPDset ge expectedSizePDset-1;
outputI, outputS, outputIp, outputSp, outputQ := PDSetZpHadamardCode(p, type);
assert #outputSp eq expectedSizePDset;
assert #outputS eq expectedSizePDset;
assert outputI eq I;
assert outputIp eq Ip;
// It takes too long...
//assert IsZpPermutationDecodeSet(C, outputIp, outputSp, expectedSizePDset-1);
//assert IsZpPermutationDecodeSet(C, outputI, outputS, expectedSizePDset-1);
assert #[p : p in outputS | C^p eq C] eq expectedSizePDset;
assert outputSp eq [PermZpsToPermZp(perm, p, #type) : perm in outputS];
testRowsPDSetZpHadamardCode(p, type, outputQ);

/****************************************************************/
print "test 12b: Z/16-linear Hadamard code of type (64; 1,1,1,1) using general construction from type (2,1,1)";

p := 2;
type := [1,1,1,1];
C := ZpHadamardCode(p, type);
I, Ip := ZpInformationSet(C);

expectedSizePDset := 16;
upperBoundPDset := MaximumSizePDSetHadamard(p, type);
assert upperBoundPDset ge expectedSizePDset-1;
outputI, outputS, outputIp, outputSp, outputQ := PDSetZpHadamardCode(p, type);
assert #outputSp eq expectedSizePDset;
assert #outputS eq expectedSizePDset;
assert outputI eq I;
assert outputIp eq Ip;
// It takes too long...
//assert IsZpPermutationDecodeSet(C, outputIp, outputSp, expectedSizePDset-1);
//assert IsZpPermutationDecodeSet(C, outputI, outputS, expectedSizePDset-1);
assert #[p : p in outputS | C^p eq C] eq expectedSizePDset;
assert outputSp eq [PermZpsToPermZp(perm, p, #type) : perm in outputS];
testRowsPDSetZpHadamardCode(p, type, outputQ);

/****************************************************************/
print "test 13: Z/9-linear Hadamard code of type (729; 3,2) using general construction";

p := 3;
type := [3,2];
C := ZpHadamardCode(p, type);
I, Ip := ZpInformationSet(C);

expectedSizePDset := 135;
upperBoundPDset := MaximumSizePDSetHadamard(p, type);
assert upperBoundPDset ge expectedSizePDset-1;
outputI, outputS, outputIp, outputSp, outputQ := PDSetZpHadamardCode(p, type);
assert #outputSp eq expectedSizePDset;
assert #outputS eq expectedSizePDset;
assert outputI eq I;
assert outputIp eq Ip;
// It takes too long...
//assert IsZpPermutationDecodeSet(C, outputIp, outputSp, expectedSizePDset-1);
//assert IsZpPermutationDecodeSet(C, outputI, outputS, expectedSizePDset-1);
assert #[p : p in outputS | C^p eq C] eq expectedSizePDset;
assert outputSp eq [PermZpsToPermZp(perm, p, #type) : perm in outputS];
testRowsPDSetZpHadamardCode(p, type, outputQ);

/****************************************************************/
print "test 13b: Z/81-linear Hadamard code of type (729; 1,0,2,2) using general construction from type (3,2)";

p := 3;
type := [1,0,2,2];
C := ZpHadamardCode(p, type);
I, Ip := ZpInformationSet(C);

expectedSizePDset := 135;
upperBoundPDset := MaximumSizePDSetHadamard(p, type);
assert upperBoundPDset ge expectedSizePDset-1;
outputI, outputS, outputIp, outputSp, outputQ := PDSetZpHadamardCode(p, type);
assert #outputSp eq expectedSizePDset;
assert #outputS eq expectedSizePDset;
assert outputI eq I;
assert outputIp eq Ip;
// It takes too long...
//assert IsZpPermutationDecodeSet(C, outputIp, outputSp, expectedSizePDset-1);
//assert IsZpPermutationDecodeSet(C, outputI, outputS, expectedSizePDset-1);
assert #[p : p in outputS | C^p eq C] eq expectedSizePDset;
assert outputSp eq [PermZpsToPermZp(perm, p, #type) : perm in outputS];
testRowsPDSetZpHadamardCode(p, type, outputQ);

/****************************************************************/
print "test 14: Z/9-linear Hadamard code of type (27; 2,1) using general construction";

p := 3;
type := [2,1];
C := ZpHadamardCode(p, type);
I, Ip := ZpInformationSet(C);

expectedSizePDset := 6;
upperBoundPDset := MaximumSizePDSetHadamard(p, type);
assert upperBoundPDset ge expectedSizePDset-1;
outputI, outputS, outputIp, outputSp, outputQ := PDSetZpHadamardCode(p, type);
assert #outputSp eq expectedSizePDset;
assert #outputS eq expectedSizePDset;
assert outputI eq I;
assert outputIp eq Ip;
// It takes too long...
//assert IsZpPermutationDecodeSet(C, outputIp, outputSp, expectedSizePDset-1);
assert IsZpPermutationDecodeSet(C, outputI, outputS, expectedSizePDset-1);
assert #[p : p in outputS | C^p eq C] eq expectedSizePDset;
assert outputSp eq [PermZpsToPermZp(perm, p, #type) : perm in outputS];
testRowsPDSetZpHadamardCode(p, type, outputQ);

/****************************************************************/
print "test 14b: Z/81-linear Hadamard code of type (27; 1,0,1,1) using general construction from type (2,1)";

p := 3;
type := [1,0,1,1];
C := ZpHadamardCode(p, type);
I, Ip := ZpInformationSet(C);

expectedSizePDset := 6;
upperBoundPDset := MaximumSizePDSetHadamard(p, type);
assert upperBoundPDset ge expectedSizePDset-1;
outputI, outputS, outputIp, outputSp, outputQ := PDSetZpHadamardCode(p, type);
assert #outputSp eq expectedSizePDset;
assert #outputS eq expectedSizePDset;
assert outputI eq I;
assert outputIp eq Ip;
// It takes too long...
//assert IsZpPermutationDecodeSet(C, outputIp, outputSp, expectedSizePDset-1);
assert IsZpPermutationDecodeSet(C, outputI, outputS, expectedSizePDset-1);
assert #[p : p in outputS | C^p eq C] eq expectedSizePDset;
assert outputSp eq [PermZpsToPermZp(perm, p, #type) : perm in outputS];
testRowsPDSetZpHadamardCode(p, type, outputQ);

/****************************************************************/
print "test 15: Z/4-linear Hadamard code of type (16; 3,0) using a non-deterministic method";

p := 2;
type := [3,0];
C := ZpHadamardCode(p, type);
I, Ip := ZpInformationSet(C);

expectedSizePDset := 5;
upperBoundPDset := MaximumSizePDSetHadamard(p, type);
assert upperBoundPDset ge expectedSizePDset-1;
outputI, outputS, outputIp, outputSp, outputQ := PDSetZpHadamardCode(p, type : 
                                                   AlgMethod := "Nondeterministic");
assert #outputSp eq expectedSizePDset;
assert #outputS eq expectedSizePDset;
assert outputI eq I;
assert outputIp eq Ip;
assert #[p : p in outputS | C^p eq C] eq expectedSizePDset;
assert outputSp eq [PermZpsToPermZp(perm, p, #type) : perm in outputS];
testRowsPDSetZpHadamardCode(p, type, outputQ);

/****************************************************************/
print "test 16: Z/81-linear Hadamard code of type (81; 1,0,2,0) using a non-deterministic method";

p := 3;
type := [1,0,2,0];
C := ZpHadamardCode(p, type);
I, Ip := ZpInformationSet(C);

expectedSizePDset := 27;
upperBoundPDset := MaximumSizePDSetHadamard(p, type);
assert upperBoundPDset ge expectedSizePDset-1;
outputI, outputS, outputIp, outputSp, outputQ := PDSetZpHadamardCode(p, type : 
                                                   AlgMethod := "Nondeterministic");
assert #outputSp eq expectedSizePDset;
assert #outputS eq expectedSizePDset;
assert outputI eq I;
assert outputIp eq Ip;
assert #[p : p in outputS | C^p eq C] eq expectedSizePDset;
assert outputSp eq [PermZpsToPermZp(perm, p, #type) : perm in outputS];
testRowsPDSetZpHadamardCode(p, type, outputQ);

/****************************************************************/
print "test 17: Z/27-linear Hadamard code of type (16; 2,0,1) using a non-deterministic method";
print("This test may take some minutes. Please wait...");

p := 3;
type := [2,0,1];
C := ZpHadamardCode(p, type);
I, Ip := ZpInformationSet(C);

upperBoundPDset := MaximumSizePDSetHadamard(p, type);
outputI, outputS, outputIp, outputSp, outputQ := PDSetZpHadamardCode(p, type : 
                                                   AlgMethod := "Nondeterministic");
assert (#outputSp-1) le upperBoundPDset;
assert (#outputS-1) le upperBoundPDset;
assert outputI eq I;
assert outputIp eq Ip;
assert #[p : p in outputS | C^p eq C] eq #outputS;
assert outputSp eq [PermZpsToPermZp(perm, p, #type) : perm in outputS];
testRowsPDSetZpHadamardCode(p, type, outputQ);

/****************************************************************/
print "test 18: Z/16-linear Hadamard code of type (32; 1,0,2,1) using a non-deterministic method";
print("This test may take some minutes. Please wait...");

p := 2;
type := [1,0,2,1];
C := ZpHadamardCode(p, type);
I, Ip := ZpInformationSet(C);

upperBoundPDset := MaximumSizePDSetHadamard(p, type);
outputI, outputS, outputIp, outputSp, outputQ := PDSetZpHadamardCode(p, type : 
                                                   AlgMethod := "Nondeterministic");
assert (#outputSp-1) le upperBoundPDset;
assert (#outputS-1) le upperBoundPDset;
assert outputI eq I;
assert outputIp eq Ip;
assert #[p : p in outputS | C^p eq C] eq #outputS;
assert outputSp eq [PermZpsToPermZp(perm, p, #type) : perm in outputS];
testRowsPDSetZpHadamardCode(p, type, outputQ);

/****************************************************************/
print "test 19: Z/4-linear Hadamard codes of type (n; t1,t2) using PDSetHadamardCodeZ4 function";

p := 2; t1 := 3; t2 := 0;
type := [t1,t2];
C := ZpHadamardCode(p, type);
I, S, Ibin, Sbin := PDSetZpHadamardCode(p, type);
I_Z4, S_Z4, Ibin_Z4, Sbin_Z4 := PDSetHadamardCodeZ4(t1, 2*t1+t2-1);

assert Set(I) eq Set(I_Z4);
assert Set(Ibin) eq Set(Ibin_Z4);
assert S eq S_Z4;
assert Sbin eq Sbin_Z4;
assert IsZpPermutationDecodeSet(C, I, S, #S-1);
assert IsZpPermutationDecodeSet(C, I_Z4, S_Z4, #S_Z4-1);

/****************************************************************/
p := 2; t1 := 3; t2 := 1;
type := [t1,t2];
C := ZpHadamardCode(p, type);
I, S, Ibin, Sbin := PDSetZpHadamardCode(p, type);
I_Z4, S_Z4, Ibin_Z4, Sbin_Z4 := PDSetHadamardCodeZ4(t1, 2*t1+t2-1);

assert Set(I) eq Set(I_Z4);
assert Set(Ibin) eq Set(Ibin_Z4);
assert IsZpPermutationDecodeSet(C, I, S, #S-1);
assert IsZpPermutationDecodeSet(C, I_Z4, S_Z4, #S_Z4-1);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////                 PERMUTATION DECODING                     ///////////////
///////                  FOR HADAMARD CODES                      ///////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/****************************************************************/
print "test 20: Z/27-linear Hadamard code of type (27; 2,0,0)";

p := 3;
type := [2,0,0];
s := #type;
C := ZpHadamardCode(p, type);
U := RSpace(Integers(p^s), Length(C));
Up := VectorSpace(GF(p), p^(s-1)*Length(C));
grayMap := CarletGrayMap(p, s);

// u in C and up in Phi(C)
u := Random(C);
up := Up!&cat[Eltseq(grayMap(i)) : i in Eltseq(u)];

// Error of weight 12 in positions {1,11}
e1 := U![1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
u1 := u+e1;
u1p := Up!&cat[Eltseq(grayMap(i)) : i in Eltseq(u1)];

// Error of weight 12 in positions {2,22}
e2 := U![0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,13,0,0,0,0,0];
u2 := u+e2;
u2p := Up!&cat[Eltseq(grayMap(i)) : i in Eltseq(u2)];

// Error of weight 18 in positions {3,10,15}
e3 := U![0,0,11,0,0,0,0,0,0,20,0,0,0,0,23,0,0,0,0,0,0,0,0,0,0,0,0];
u3 := u+e3;
u3p := Up!&cat[Eltseq(grayMap(i)) : i in Eltseq(u3)];

// Sequences of additive tuples over Zps and vectors over Zp
Q := [u1, u2, u3];
Qp := [u1p, u2p, u3p];

I, S, Ip, Sp := PDSetZpHadamardCode(p, type); 
r := #S-1;

OutputIsDecoded_u, OutputDecoded_u, OutputDecoded_up := 
                                            ZpPermutationDecode(C, Ip, Sp, r, u1p);
assert OutputIsDecoded_u eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_up eq up;

OutputIsDecoded_u, OutputDecoded_u, OutputDecoded_up := 
                                            ZpPermutationDecode(C, Ip, Sp, r, u2p);
assert OutputIsDecoded_u eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_up eq up;

OutputIsDecoded_u, OutputDecoded_u, OutputDecoded_up := 
                                            ZpPermutationDecode(C, Ip, Sp, r, u3p);
assert OutputIsDecoded_u eq false;
assert OutputDecoded_u eq U!0;
assert OutputDecoded_up eq Up!0;

OutputIsDecoded_u, OutputDecoded_u, OutputDecoded_up := 
                                            ZpPermutationDecode(C, Ip, Sp, r, u1);
assert OutputIsDecoded_u eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_up eq up;

OutputIsDecoded_u, OutputDecoded_u, OutputDecoded_up := 
                                            ZpPermutationDecode(C, Ip, Sp, r, u2);
assert OutputIsDecoded_u eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_up eq up;

OutputIsDecoded_u, OutputDecoded_u, OutputDecoded_up := 
                                            ZpPermutationDecode(C, Ip, Sp, r, u3);
assert OutputIsDecoded_u eq false;
assert OutputDecoded_u eq U!0;
assert OutputDecoded_up eq Up!0;

OutputIsDecoded_Q, OutputDecoded_Q, OutputDecoded_Qp := 
                                            ZpPermutationDecode(C, Ip, Sp, r, Qp);
assert OutputIsDecoded_Q eq [true, true, false];
assert OutputDecoded_Q eq [u, u, U!0];
assert OutputDecoded_Qp eq [up, up, Up!0];

OutputIsDecoded_Q, OutputDecoded_Q, OutputDecoded_Qp := 
                                            ZpPermutationDecode(C, Ip, Sp, r, Q);
assert OutputIsDecoded_Q eq [true, true, false];
assert OutputDecoded_Q eq [u, u, U!0];
assert OutputDecoded_Qp eq [up, up, Up!0];

/****************************************************************/
print "test 21: Z/8-linear Hadamard code of type (64; 1,3,0)";

p := 2;
type := [1,3,0];
s := #type;
C := ZpHadamardCode(p, type);
U := RSpace(Integers(p^s), Length(C));
Up := VectorSpace(GF(p), p^(s-1)*Length(C));
grayMap := CarletGrayMap(p, s);

// u in C and up in Phi(C)
u := Random(C);
up := Up!&cat[Eltseq(grayMap(i)) : i in Eltseq(u)];

// Error of weight <= 15 in positions {2,3,11,12}
e1 := U!0; 
e1[2] := Random(Integers(p^s)); 
e1[3] := Random(Integers(p^s)); 
e1[11] := Random(Integers(p^s));
e1[12] := 3;
u1 := u+e1;
u1p := Up!&cat[Eltseq(grayMap(i)) : i in Eltseq(u1)];

// Error of weight 16 > 15 in positions {2,3,11,12}
e2 := U!0; 
e2[2] := 4; 
e2[3] := 4; 
e2[11] := 4;
e2[12] := 4;
u2 := u+e2;
u2p := Up!&cat[Eltseq(grayMap(i)) : i in Eltseq(u2)];

// Sequences of additive tuples over Zps and vectors over Zp
Q := [u1, u1, u2];
Qp := [u1p, u1p, u2p];

I, S, Ip, Sp := PDSetZpHadamardCode(p, type); 
r := #S-1;

OutputIsDecoded_u, OutputDecoded_u, OutputDecoded_up := 
                                            ZpPermutationDecode(C, Ip, Sp, r, u1p);
assert OutputIsDecoded_u eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_up eq up;

OutputIsDecoded_u, OutputDecoded_u, OutputDecoded_up := 
                                            ZpPermutationDecode(C, Ip, Sp, r, u2p);
assert OutputIsDecoded_u eq false;
assert OutputDecoded_u eq U!0;
assert OutputDecoded_up eq Up!0;

OutputIsDecoded_u, OutputDecoded_u, OutputDecoded_up := 
                                            ZpPermutationDecode(C, Ip, Sp, r, u1);
assert OutputIsDecoded_u eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_up eq up;

OutputIsDecoded_u, OutputDecoded_u, OutputDecoded_up := 
                                            ZpPermutationDecode(C, Ip, Sp, r, u2);
assert OutputIsDecoded_u eq false;
assert OutputDecoded_u eq U!0;
assert OutputDecoded_up eq Up!0;

OutputIsDecoded_Q, OutputDecoded_Q, OutputDecoded_Qp := 
                                            ZpPermutationDecode(C, Ip, Sp, r, Qp);
assert OutputIsDecoded_Q eq [true, true, false];
assert OutputDecoded_Q eq [u, u, U!0];
assert OutputDecoded_Qp eq [up, up, Up!0];

OutputIsDecoded_Q, OutputDecoded_Q, OutputDecoded_Qp := 
                                            ZpPermutationDecode(C, Ip, Sp, r, Q);
assert OutputIsDecoded_Q eq [true, true, false];
assert OutputDecoded_Q eq [u, u, U!0];
assert OutputDecoded_Qp eq [up, up, Up!0];

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////         PD-SETS FOR CODES WHICH ARE NOT HADAMARD         ///////////////
///////               AND PERMUTATION DECODING                   ///////////////  
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/****************************************************************/
print "test 22a: Z/4-linear Kerdock code of type (16; 5,0)";

C := KerdockCode(4); 

I, Ip := ZpInformationSet(C);
tau := Sym(Length(C))!([2..15] cat [1,16]);
S := [tau^5, tau^10, tau^15];
Sp := [PermZpsToPermZp(perm, 2, 2) : perm in S];
r := #S-1;
assert IsZpPermutationDecodeSet(C, I, S, r);
assert IsZpPermutationDecodeSet(C, Ip, Sp, r);

p := 2; s := 2;
U := RSpace(Integers(p^s), Length(C));
Up := VectorSpace(GF(p), p^(s-1)*Length(C));
grayMap := CarletGrayMap(p, s);
// u in C and up in Phi(C)
u := Random(C);
up := Up!&cat[Eltseq(grayMap(i)) : i in Eltseq(u)];
// Error of weight <= 15 in positions {2,3,11,12}
e1 := U!0; 
e1[2] := 1; 
e1[3] := 3; 
u1 := u+e1;
u1p := Up!&cat[Eltseq(grayMap(i)) : i in Eltseq(u1)];

OutputIsDecoded_u, OutputDecoded_u, OutputDecoded_up := 
                                            ZpPermutationDecode(C, Ip, S, r, u1);
assert OutputIsDecoded_u eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_up eq up;

OutputIsDecoded_u, OutputDecoded_u, OutputDecoded_up := 
                                            ZpPermutationDecode(C, Ip, S, r, u1p);
assert OutputIsDecoded_u eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_up eq up;

/****************************************************************/
print "test 22b: Z/4-linear Kerdock code of type (32; 6,0)";

C := KerdockCode(5); 

I, Ip := ZpInformationSet(C);
tau := Sym(Length(C))!(1,32,9,19,25)(2,18,24,15,31)(3,27,23,28,12)
                       (4,8,20,30,26)(5,14,16,21,13)(6,10,17,29,22);
S := [tau^i : i in [1..5]];
Sp := [PermZpsToPermZp(perm, 2, 2) : perm in S];
r := #S-1;
assert IsZpPermutationDecodeSet(C, I, S, r);
assert IsZpPermutationDecodeSet(C, Ip, Sp, r);

/****************************************************************/
print "test 23: Z/4-linear code of type (15; 0,5)";

Z4 := Integers(4);
G := Matrix(Z4,[[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2], 
		        [0,2,0,0,0,2,0,0,2,2,0,2,0,2,2,2],
                [0,0,2,0,0,2,2,0,2,0,2,2,2,2,0,0],
                [0,0,0,2,0,0,2,2,0,2,0,2,2,2,2,0],
                [0,0,0,0,2,0,0,2,2,0,2,0,2,2,2,2]]);
C := LinearCode(G);

I, Ip := ZpInformationSet(C);
p1 := Sym(16)!1;
p2 := Sym(16)!(1, 14, 11, 9, 6, 10, 13, 3, 15, 5, 16, 2, 12, 8)(4, 7);
p3 := Sym(16)!(1, 14, 11, 2, 7, 9, 5, 12, 3, 16, 13, 6)(4, 15, 8, 10);
S := [p1, p2, p3];
Sp := [PermZpsToPermZp(perm, 2, 2) : perm in S];
r := #S-1;
assert IsZpPermutationDecodeSet(C, I, S, r);
assert IsZpPermutationDecodeSet(C, Ip, Sp, r);

p := 2; s := 2;
U := RSpace(Integers(p^s), Length(C));
Up := VectorSpace(GF(p), p^(s-1)*Length(C));
grayMap := CarletGrayMap(p, s);
// u in C and up in Phi(C)
u := Random(C);
up := Up!&cat[Eltseq(grayMap(i)) : i in Eltseq(u)];
// Error of weight <= 15 in positions {2,3,11,12}
e1 := U!0; 
e1[4] := 1; 
e1[7] := 3; 
u1 := u+e1;
u1p := Up!&cat[Eltseq(grayMap(i)) : i in Eltseq(u1)];

// u in C
u := U![0,2,0,0,0,2,0,0,2,2,0,2,0,2,2,2];  
// v not in C, 2 (quaternary) errors in positions {2, 8};
v := U![0,1,0,0,0,2,0,3,2,2,0,2,0,2,2,2];  
// w not in C, 3 (quaternary) errors in positions {2, 8, 6};
w := U![0,1,0,0,0,1,0,3,2,2,0,2,0,2,2,2]; 

// ubin in Cbin
ubin := Up![0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,1,1,1,1,0,0,1,1,0,0,1,1,1,1,1,1];
// vbin not in Cbin, 2 errors in positions {3, 15}
vbin := Up![0,0,0,1,0,0,0,0,0,0,1,1,0,0,1,0,1,1,1,1,0,0,1,1,0,0,1,1,1,1,1,1];
// wbin not in Cbin, 3 errors in positions {3, 15, 11}
wbin := Up![0,0,0,1,0,0,0,0,0,0,0,1,0,0,1,0,1,1,1,1,0,0,1,1,0,0,1,1,1,1,1,1];

/****************************************************************/
print "test 24: (nonlinear) Z/4-linear code of type (15; 6,1),
                Gray map image of a quaternary linear cyclic code";

Z4 := Integers(4);
PR4<y> := PolynomialRing(Z4);   // Cyclic Linear Code C over Z4
C := CyclicCode(15, y^9 + 3*y^7 + y^6 + y^3 + 3*y^2 + 1); 

I, Ip := ZpInformationSet(C);
tau := Sym(15)!(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15);
S := [tau^i : i in [1..15]];
Sp := [PermZpsToPermZp(perm, 2, 2) : perm in S];
assert IsZpPermutationDecodeSet(C, I, S, 2);
assert not IsZpPermutationDecodeSet(C, I, S, 3);
assert IsZpPermutationDecodeSet(C, Ip, Sp, 2);
