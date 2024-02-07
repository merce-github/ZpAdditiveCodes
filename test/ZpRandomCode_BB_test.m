/************************************************************/
/*                                                          */
/* Project name: Z/p^s-additive codes in MAGMA              */
/* Test file name: ZpRandomCode_BB_test.m                   */
/*                                                          */
/* Comments: Black-box tests for the functions              */
/*           RandomZpAdditiveCode included in the           */
/*           ZpAdditiveCodes_Constructions.m file           */
/*                                                          */                                                
/* Authors: Noam von Rotberg and Mercè Villanueva           */
/*                                                          */
/* Revision version and last date: v1.0    2021/08/28       */
/*                                 v1.1    2023/02/12       */
/*                                                          */
/************************************************************/

SetAssertions(true);
//Alarm(30*60);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////                   ZERO CODES OVER Z/p^s                         ////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/****************************************************************/
/*                                                              */
/* Function name: TestRandomCode_Zero                           */
/* Parameters: p, n, s, testNumber                              */
/* Description: Given a prime p, two non-negative integers n    */
/*   and s, and an integer with the testNumber, tests whether   */
/*   RandomZpAdditiveCode(p, n, [0,...,0]) is the zero code.    */
/*   In particular, it checks the following:                    */         
/*   - output code is a zero code over Z/p^s of length n        */
/*   - output code is generated by the output matrix            */ 
/*   - output code has correct type                             */
/*   - output code has correct length                           */
/*   - output code has correct base ring                        */
/*   - output matrix has correct base ring                      */
/*   - output matrix is in standard form                        */
/* Input parameters description:                                */
/*   - p: a prime                                               */
/*   - n: a non-negative integer                                */
/*   - s: an integer > 1 describing the length of [0,...,0]     */
/*   - testNumber: an integer being printed as test number      */
/* Output parameters description:                               */
/*   - testNumber: the by one increased integer                 */
/*                                                              */
/****************************************************************/
TestRandomCode_Zero := function(p, n, s, testNumber)
    Zps := Integers(p^s);
    type := [0^^s];
    expectedOutputCode := ZeroCode(Zps, n);
    OutputCode, OutputMatrix := RandomZpAdditiveCode(p, n, type);

    print "Test 1 .", testNumber,": Zero code over Z /", p, "^", s, "of length", n;
    assert OutputCode eq expectedOutputCode;
    assert OutputCode eq LinearCode(OutputMatrix);
    assert ZpType(OutputCode) eq type; 
    assert Length(OutputCode) eq n;
    assert BaseRing(OutputCode) eq Zps;
    assert BaseRing(OutputMatrix) eq Zps;
    assert IsStandardFormMatrix(OutputMatrix);

    return testNumber + 1;
end function;

/****************************************************************/
// counter for test number
testNumber := 1; 

p := 2;
n := 10;
s := 3;

// Test RandomZpAdditiveCode function on this zero code
testNumber := TestRandomCode_Zero(p, n, s, testNumber);

/****************************************************************/
p := 3;
n := 9;
s := 3;

// Test RandomZpAdditiveCode function on this zero code
testNumber := TestRandomCode_Zero(p, n, s, testNumber);

/****************************************************************/
p := 5;
n := 13;
s := 4;

// Test RandomZpAdditiveCode function on this zero code
testNumber := TestRandomCode_Zero(p, n, s, testNumber);

/****************************************************************/
p := 2;
n := 20;
s := 4;

// Test RandomZpAdditiveCode function on this zero code
testNumber := TestRandomCode_Zero(p, n, s, testNumber);

/****************************************************************/
p := 5;
n := 3;
s := 8;

// Test RandomZpAdditiveCode function on this zero code
testNumber := TestRandomCode_Zero(p, n, s, testNumber);

/****************************************************************/
p := 131;
n := 30;
s := 11;

// Test RandomZpAdditiveCode function on this zero code
testNumber := TestRandomCode_Zero(p, n, s, testNumber);

/****************************************************************/
p := 11;
n := 2;
s := 2;

// Test RandomZpAdditiveCode function on this zero code
testNumber := TestRandomCode_Zero(p, n, s, testNumber);

/****************************************************************/
p := 7;
n := 5;
s := 3;

// Test RandomZpAdditiveCode function on this zero code
testNumber := TestRandomCode_Zero(p, n, s, testNumber);

/****************************************************************/
p := 3;
n := 7;
s := 5;

// Test RandomZpAdditiveCode function on this zero code
testNumber := TestRandomCode_Zero(p, n, s, testNumber);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////                  RANDOM CODES OVER Z/p^s                        ////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/****************************************************************/
/*                                                              */
/* Function name: TestRandomCode_General                        */
/* Parameters: p, n, type, testNumber                           */
/* Description: Given a prime p, a non-negative integer n, a    */
/*   sequence of integers and an integer with the testNumber,   */
/*   test whether RandomZpAdditiveCode(p, n, type) gives a      */
/*   reasonable output by checking the following:               */
/*   - output code is generated by the output matrix            */ 
/*   - output code has correct type                             */
/*   - output code has correct length                           */
/*   - output code has correct base ring                        */
/*   - output matrix has correct base ring                      */ 
/*   - output matrix is in standard form                        */
/* Input parameters description:                                */
/*   - p: a prime                                               */
/*   - n: a non-negative integer                                */
/*   - type: a non-empty list containing non-negative integers  */
/*   - testNumber: an integer being printed as test number      */
/* Output parameters description:                               */
/*   - testNumber: the by one increased integer                 */
/*                                                              */
/****************************************************************/
TestRandomCode_General := function(p, n, type, testNumber)
    Zps := Integers(p^#type);
    
    print "Test 2 .", testNumber, ": Random code over Z /", p, "^", #type, "of length", n, "and type", type;
    OutputCode, OutputMatrix := RandomZpAdditiveCode(p, n, type);

    assert OutputCode eq LinearCode(OutputMatrix);
    assert ZpType(OutputCode) eq type;
    assert Length(OutputCode) eq n;
    assert BaseRing(OutputCode) eq Zps;
    assert BaseRing(OutputMatrix) eq Zps;
    assert IsStandardFormMatrix(OutputMatrix); 

    k := &+type;
    print "Test 3 .", testNumber, ": Random code over Z /", p, "^", #type, "of length", n, "and pseudo-dimension at least", k;
    OutputCode, OutputMatrix := RandomZpAdditiveCode(p, n, #type, k);
    assert OutputCode eq LinearCode(OutputMatrix);
    assert &+ZpType(OutputCode) eq k;
    assert Length(OutputCode) eq n;
    assert BaseRing(OutputCode) eq Zps;
    assert BaseRing(OutputMatrix) eq Zps;
    assert IsStandardFormMatrix(OutputMatrix); 

    return testNumber + 1;
end function;

/****************************************************************/
// counter for printing increasing test number
testNumber := 1; 

/****************************************************************/
p := 2;
n := 10;
type := [1, 1, 1];

// Test RandomZpAdditiveCode function on this random code
testNumber := TestRandomCode_General(p, n, type, testNumber);

/****************************************************************/
p := 3;
n := 12;
type := [4, 4];

// Test RandomZpAdditiveCode function on this random code
testNumber := TestRandomCode_General(p, n, type, testNumber);

/****************************************************************/
p := 3;
n := 9;
type := [4, 2, 3];

// Test RandomZpAdditiveCode function on this random code
testNumber := TestRandomCode_General(p, n, type, testNumber);

/****************************************************************/
p := 5;
n := 13;
type := [1, 2, 3, 4];

// Test RandomZpAdditiveCode function on this random code
testNumber := TestRandomCode_General(p, n, type, testNumber);

/****************************************************************/
p := 103;
n := 68;
type := [1, 0, 10, 3, 5, 2, 8, 7, 0, 4];

// Test RandomZpAdditiveCode function on this random code
testNumber := TestRandomCode_General(p, n, type, testNumber);

/****************************************************************/
p := 7;
n := 25;
type := [8, 15];

// Test RandomZpAdditiveCode function on this random code
testNumber := TestRandomCode_General(p, n, type, testNumber);

/****************************************************************/
p := 2999;
n := 555;
type := [8, 7, 6, 5, 4, 3, 2, 1];

// Test RandomZpAdditiveCode function on this random code
testNumber := TestRandomCode_General(p, n, type, testNumber);

/****************************************************************/
p := 3;
n := 8;
type := [0, 7, 1, 0];

// Test RandomZpAdditiveCode function on this random code
testNumber := TestRandomCode_General(p, n, type, testNumber);

/****************************************************************/
p := 5;
n := 34;
type := [0, 0, 0, 9, 14, 8];

// Test RandomZpAdditiveCode function on this random code
testNumber := TestRandomCode_General(p, n, type, testNumber); 
