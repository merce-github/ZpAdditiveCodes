/**************************************************************/
/*                                                            */
/* Project name: Z/p^s-additive codes in MAGMA                */
/* Test file name: ZpDual_BB_test.m                           */
/*                                                            */
/* Comments: Black-box tests for the function ZpDual          */
/*           included in the ZpAdditiveCodes_Core.m file      */
/*                                                            */
/* Authors: Adrián Torres and Mercè Villanueva                */
/*                                                            */
/* Revision version and last date: v1.0    2021/08/27         */
/*                                 v1.1    2021/09/29         */
/*                                 v1.2    2024/01/12         */
/*                                                            */
/**************************************************************/

SetAssertions(true);
// Alarm(30*60);
// SetQuitOnError(true);

/****************************************************************/
print "test 0: Zero code over Z3^3 of type (10; 0,0,0)";

p := 3;
s := 3;
n := 10;
Zps := Integers(p^s);
C := ZeroCode(Zps, n);

expectedOutputDual := UniverseCode(Zps, n);
outputDual := ZpDual(C);
assert outputDual eq expectedOutputDual;

G := GeneratorMatrix(C);
Gd := GeneratorMatrix(outputDual);
assert NNZEntries(G*Transpose(Gd)) eq 0;

assert ZpType(C) eq [0^^s];
assert ZpTypeDual(C) eq ZpType(outputDual);

/****************************************************************/
print "test 1: Universe code over Z3^3 of type (10; 10,0,0)";

p := 3;
s := 3;
n := 10;
Zps := Integers(p^s);
C := UniverseCode(Zps, n);

expectedOutputDual := ZeroCode(Zps, n);
outputDual := ZpDual(C);
assert outputDual eq expectedOutputDual;

G := GeneratorMatrix(C);
Gd := GeneratorMatrix(outputDual);
assert NNZEntries(G*Transpose(Gd)) eq 0;

assert ZpType(C) eq [10, 0^^(s-1)];
assert ZpTypeDual(C) eq ZpType(outputDual);

/****************************************************************/
print "test 2: Linear code over Z2^3 of type (8; 2,2,2)";

p := 2;
s := 3;
n := 8;
Zps := Integers(p^s);
C := LinearCode<Zps, n | [[1,0,0,0,2,3,4,6],
                          [0,1,0,1,3,3,7,3],
                          [0,0,2,0,2,0,2,2],
                          [0,0,0,2,2,0,2,6],
                          [0,0,0,0,4,0,4,4],
                          [0,0,0,0,0,4,0,0]]>;

expectedOutputDual := Dual(C);
outputDual := ZpDual(C);
assert outputDual eq expectedOutputDual;

G := GeneratorMatrix(C);
Gd := GeneratorMatrix(outputDual);
assert NNZEntries(G*Transpose(Gd)) eq 0;

assert ZpType(C) eq [2,2,2];
assert ZpTypeDual(C) eq ZpType(outputDual);

/****************************************************************/
print "test 3: Linear code over Z3^4 of type (10; 3,1,2,1)";

p := 3;
s := 4;
n := 10;
Zps := Integers(p^s);
C := LinearCode<Zps, n | [[1,0,0,1,0,4,22,49,21,80],
                          [0,1,0,0,6,7,22,53,37,26],
                          [0,0,1,0,4,3,26,17,67,52],
                          [0,0,0,3,3,3,0,12,39,27],
                          [0,0,0,0,9,0,9,36,27,36],
                          [0,0,0,0,0,9,9,72,0,72],
                          [0,0,0,0,0,0,27,54,27,54]]>;

expectedOutputDual := Dual(C);
outputDual := ZpDual(C);
assert outputDual eq expectedOutputDual;

G := GeneratorMatrix(C);
Gd := GeneratorMatrix(outputDual);
assert NNZEntries(G*Transpose(Gd)) eq 0;

assert ZpType(C) eq [3,1,2,1];
assert ZpTypeDual(C) eq ZpType(outputDual);

/****************************************************************/
print "test 4: Simple linear code over Z2^3 of type (3; 1,1,1) with generator matrix not in standard form";

p := 2;
s := 3;
n := 3;
Zps := Integers(p^s);
C := LinearCode<Zps, n | [[0,1,0],
                          [0,0,2],
                          [4,0,0]]>;

expectedOutputDual := Dual(C);
expectedOutputDualManual := LinearCode<Zps, n | [[2,0,0],[0,0,4]]>;
assert expectedOutputDual eq expectedOutputDualManual;
outputDual := ZpDual(C);
assert outputDual eq expectedOutputDual;

G := GeneratorMatrix(C);
Gd := GeneratorMatrix(outputDual);
assert NNZEntries(G*Transpose(Gd)) eq 0;

assert ZpType(C) eq [1,1,1];
assert ZpTypeDual(C) eq ZpType(outputDual);

/****************************************************************/
print "test 5: Linear code over Z5^3 of type (10; 1,2,1) with generator matrix not in standard form";

p := 5;
s := 3;
n := 10;
Zps := Integers(p^s);
C := LinearCode<Zps, n | [[1,49,24,14,2,4,58,2,105,44],
                          [0,115,100,5,5,0,40,30,95,85],
                          [0,70,105,5,0,5,25,120,115,55],
                          [0,0,50,25,0,0,75,25,50,75]]>;

expectedOutputDual := Dual(C);
outputDual := ZpDual(C);
assert outputDual eq expectedOutputDual;

G := GeneratorMatrix(C);
Gd := GeneratorMatrix(outputDual);
assert NNZEntries(G*Transpose(Gd)) eq 0;

assert ZpType(C) eq [1,2,1];
assert ZpTypeDual(C) eq ZpType(outputDual);

/****************************************************************/
print "test 6: Hadamard code over Z3^3 of type (27; 1,1,1) with generator matrix not in standard form";

p := 3;
type := [1,1,1];
C := ZpHadamardCode(p, type);

expectedOutputDual := Dual(C);
outputDual := ZpDual(C);
assert outputDual eq expectedOutputDual;

G := GeneratorMatrix(C);
Gd := GeneratorMatrix(outputDual);
assert NNZEntries(G*Transpose(Gd)) eq 0;

assert ZpType(C) eq type;
assert ZpTypeDual(C) eq ZpType(outputDual);

/****************************************************************/
print "test 7: Linear code over Z3^4 of type (10; 2,2,2,2) with generator matrix not in standard form";

p := 3;
s := 4;
n := 10;
Zps := Integers(p^s);
C := LinearCode<Zps, n | [[43,0,0,1,0,3,1,22,0,74],
                          [61,1,2,1,18,8,0,10,1,18],
                          [ 0,6,3,0,15,6,0,3,0,57],
                          [36,6,0,3,0,0,0,21,0,0],
                          [72,9,0,0,18,0,0,0,0,18],
                          [72,0,0,0,18,9,0,9,0,54],
                          [ 0,0,0,0,27,0,0,0,0,54],
                          [ 0,0,0,0,0,0,0,27,0,54]]>;

expectedOutputDual := Dual(C);
outputDual := ZpDual(C);
assert outputDual eq expectedOutputDual;

G := GeneratorMatrix(C);
Gd := GeneratorMatrix(outputDual);
assert NNZEntries(G*Transpose(Gd)) eq 0;

assert ZpType(C) eq [2,2,2,2];
assert ZpTypeDual(C) eq ZpType(outputDual);

/****************************************************************/
print "test 8: Hadamard code over Z4 of type (256; 4,2) with generator matrix not in standard form";

p := 2;
type := [4,2];
C := ZpHadamardCode(2, type);

expectedOutputDual := Dual(C);
outputDual := ZpDual(C);
outputDualZ4 := DualZ4(C);
assert outputDual eq expectedOutputDual;
assert outputDualZ4 eq expectedOutputDual;

G := GeneratorMatrix(C);
Gd := GeneratorMatrix(outputDual);
assert NNZEntries(G*Transpose(Gd)) eq 0;

assert ZpType(C) eq type;
assert ZpTypeDual(C) eq ZpType(outputDual);