/************************************************************/
/*                                                          */
/* Project name: Z/p^s-additive codes in MAGMA              */
/* Test file name: ZpWeightDistribution_BB_test.m           */
/*                                                          */
/* Comments: Black-box tests for the functions              */
/*           LeeWeight, HomogeneousWeight,                  */
/*           MinimumHomogeneousWeight, and                  */
/*           HomogeneousWeightDistribution                  */
/*           included in ZpAdditiveCodes_Distances.m file   */
/*                                                          */
/* Authors: Mercè Villanueva and Adrián Torres              */
/*                                                          */
/* Revision version and last date: v1.0    2021/11/24       */
/*                                 v2.0    2024/01/06       */
/*                                                          */
/************************************************************/

SetAssertions(true);
// Alarm(30*60);

SetSeed(0);

////////////////////////////////////////////////////////////////////////////////
///////                                                                 ////////
///////                     LEE/HOMOGENEOUS WEIGHT                      ////////  
///////                                                                 ////////
////////////////////////////////////////////////////////////////////////////////

/****************************************************************/
print "Testing Lee and homogeneous weights on individual components";
a := Integers(27)!3;
expectedLeeWeight := 3;
expectedHomogeneousWeight := 6;
outputLeeWeight := LeeWeight(a);
outputHomogenousWeight := HomogeneousWeight(a);
assert outputLeeWeight eq expectedLeeWeight;
assert outputHomogenousWeight eq expectedHomogeneousWeight;

a := Integers(27)!18;
expectedLeeWeight := 9;
expectedHomogeneousWeight := 9;
outputLeeWeight := LeeWeight(a);
outputHomogenousWeight := HomogeneousWeight(a);
assert outputLeeWeight eq expectedLeeWeight;
assert outputHomogenousWeight eq expectedHomogeneousWeight;

a := Integers(125)!50;
expectedLeeWeight := 50;
expectedHomogeneousWeight := 25;
outputLeeWeight := LeeWeight(a);
outputHomogenousWeight := HomogeneousWeight(a);
assert outputLeeWeight eq expectedLeeWeight;
assert outputHomogenousWeight eq expectedHomogeneousWeight;

/****************************************************************/
print "Testing Lee and homogeneous weights on vectors";
u := RSpace(Integers(81),4)![27,9,3,0];
expectedLeeWeight := 39;
expectedHomogeneousWeight := 63;
outputLeeWeight := LeeWeight(u);
outputHomogenousWeight := HomogeneousWeight(u);
assert outputLeeWeight eq expectedLeeWeight;
assert outputHomogenousWeight eq expectedHomogeneousWeight;

u := RSpace(Integers(125),4)![75,125,63,0];
expectedLeeWeight := 112;
expectedHomogeneousWeight := 45;
outputLeeWeight := LeeWeight(u);
outputHomogenousWeight := HomogeneousWeight(u);
assert outputLeeWeight eq expectedLeeWeight;
assert outputHomogenousWeight eq expectedHomogeneousWeight;

////////////////////////////////////////////////////////////////////////////////
///////                                                                 ////////
///////                         MINIMUM WEIGHT                          ////////  
///////                                                                 ////////
////////////////////////////////////////////////////////////////////////////////

print "Testing minimum homogeneous weights and distribution on linear codes over Z/p^s";

/****************************************************************/
print "test 1: Zero code over Z/16 of length 10";

C := ZeroCode(Integers(16), 10);  
expectedMinWeight := 160;
expectedWeightDist := [<0, 1>];
outputMinWeight := MinimumHomogeneousWeight(C);
outputWeightDist := HomogeneousWeightDistribution(C);
assert outputMinWeight eq expectedMinWeight;
assert outputWeightDist eq expectedWeightDist;

/****************************************************************/
print "test 2: Universe code over Z8 of length 3";

C := UniverseCode(Integers(8), 3);
expectedMinWeight := 2;
expectedWeightDist := [ <0, 1>, <2, 18>, <4, 111>, <6, 252>, 
                        <8, 111>, <10, 18>, <12, 1> ];
outputMinWeight := MinimumHomogeneousWeight(C);
outputWeightDist := HomogeneousWeightDistribution(C);
assert outputMinWeight eq expectedMinWeight;
assert outputWeightDist eq expectedWeightDist;

/****************************************************************/
print "test 3: Repetition code over Z25 of length 6";

C := RepetitionCode(Integers(25), 6);
expectedMinWeight := 24;
expectedWeightDist := [ <0, 1>, <24, 20>, <30, 4> ];
outputMinWeight := MinimumHomogeneousWeight(C);
outputWeightDist := HomogeneousWeightDistribution(C);
assert outputMinWeight eq expectedMinWeight;
assert outputWeightDist eq expectedWeightDist;

/****************************************************************/
print "test 4: Linear code over Z/125 of type [2,1,0]";
C := LinearCode<Integers(125),20 | 
     [[1,0,0,11,79,124,15,7,28,74,71,66,75,89,43,19,66,91,2,20], // Bon exemple
      [0,1,2,35,99,98,44,12,23,68,107,4,120,2,86,123,98,31,103,49],
      [0,0,5,50,10,80,0,100,0,60,35,20,35,85,95,10,10,70,20,60]]>;
expectedMinWeight := 300;
outputMinWeight := MinimumHomogeneousWeight(C);
assert outputMinWeight eq expectedMinWeight;

/****************************************************************/
print "test 5: Linear generalized Hadamard code over Z/8 of type [2,1,0]";

C := ZpHadamardCode(2, [2,1,0]);
expectedMinWeight := 64;
expectedWeightDist := [ <0, 1>, <64, 254>, <128, 1> ];
outputMinWeight := MinimumHomogeneousWeight(C);
outputWeightDist := HomogeneousWeightDistribution(C);
assert outputMinWeight eq expectedMinWeight;
assert outputWeightDist eq expectedWeightDist;

/****************************************************************/
print "test 6: Linear generalized Hadamard code over Z/9 of type [2,1]";

C := ZpHadamardCode(3, [2,1]);
expectedMinWeight := 54;
expectedWeightDist := [ <0, 1>, <54, 240>, <81, 2> ];
outputMinWeight := MinimumHomogeneousWeight(C);
outputWeightDist := HomogeneousWeightDistribution(C);
assert outputMinWeight eq expectedMinWeight;
assert outputWeightDist eq expectedWeightDist;

/****************************************************************/
print "test 7: Linear simplex alpha code over Z/27 of type [2,0,0]";

C := ZpSimplexAlphaCode(3, 3, 2);
expectedMinWeight := 4374;
expectedWeightDist := [ <0, 1>, <4374, 728> ];
outputMinWeight := MinimumHomogeneousWeight(C);
outputWeightDist := HomogeneousWeightDistribution(C);
assert outputMinWeight eq expectedMinWeight;
assert outputWeightDist eq expectedWeightDist;

/****************************************************************/
print "test 8: Linear simplex beta code over Z/27 of type [2,0,0]";

C := ZpSimplexBetaCode(3, 3, 2);
expectedMinWeight := 216;
expectedWeightDist := [ <0, 1>, <216, 720>, <243, 8> ];
outputMinWeight := MinimumHomogeneousWeight(C);
outputWeightDist := HomogeneousWeightDistribution(C);
assert outputMinWeight eq expectedMinWeight;
assert outputWeightDist eq expectedWeightDist;

/****************************************************************/
print "test 9: Linear MacDonald alpha code over Z/27 of type [2,0,0]";

C := ZpMacDonaldAlphaCode(3, 3, 2, 1);
expectedMinWeight := 4212;
expectedWeightDist := [ <0, 1>, <4212, 702>, <4374, 26> ];
outputMinWeight := MinimumHomogeneousWeight(C);
outputWeightDist := HomogeneousWeightDistribution(C);
assert outputMinWeight eq expectedMinWeight;
assert outputWeightDist eq expectedWeightDist;

/****************************************************************/
print "test 10: Linear MacDonald beta code over Z/27 of type [2,0,0]";

C := ZpMacDonaldBetaCode(3, 3, 2, 1);
expectedMinWeight := 207;
expectedWeightDist := [ <0, 1>, <207, 48>, <210, 648>, <216, 24>, <234, 6>, <243, 2> ];
outputMinWeight := MinimumHomogeneousWeight(C);
outputWeightDist := HomogeneousWeightDistribution(C);
assert outputMinWeight eq expectedMinWeight;
assert outputWeightDist eq expectedWeightDist;