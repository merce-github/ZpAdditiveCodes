/************************************************************/
/*                                                          */
/* Project name: Z/p^s-additive codes in MAGMA              */
/* Test file name: ZpResidueTorsionCode_BB_test.m           */
/*                                                          */
/* Comments: Black-box tests for the functions              */
/*           ZpResidueCode and ZpTorsionCode included in    */
/*           ZpAdditiveCodes_Constructions.m file           */
/*                                                          */                                                
/* Authors: Abdullah Irfan Basheer and Mercè Villanueva     */
/*                                                          */
/* Revision version and last date: v1.0    2023/06/20       */
/*                                 v1.1    2023/07/23       */
/*                                 v1.2    2023/08/02       */
/*                                                          */
/************************************************************/

SetAssertions(true);
//Alarm(30*60);

/****************************************************************/
/*                                                              */
/* Function name: TestZpResidueTorsionCode                      */
/* Parameters: C, checkDefinition                               */
/* Description: Given a linear code C over Z/p^s of type (n; t1,*/
/*   ...,ts) and a boolean checkDefinition, the procedure checks*/
/*   that the residue and torsion codes given by functions      */
/*   ZpResidueCode and ZpTorsionCode satisfies some properties. */
/*   It is also checked that they satisfy the definition if     */
/*   checkDefinition is set to true.                            */ 
/* Input parameters description:                                */
/*   - C : a linear code over Z/p^s of length n                 */
/*   - checkDefinition: a boolean that indicates whether it is  */
/*                checked that the codes satisfy its definition */
/*                                                              */
/****************************************************************/
TestZpResidueTorsionCode := procedure(C, checkDefinition) 
    Zps:= Alphabet(C);
    p := Factorization(#Zps)[1][1];
    Z := Integers();
    s := Valuation(#Zps, p);
    n := Length(C); 
    Vp := VectorSpace(GF(p), n);
    V := RSpace(Zps, n);
    type := ZpType(C);

    // Test for the ZpResidueCode(C) function
    outputResidueCode := ZpResidueCode(C);
    assert Dimension(outputResidueCode) eq type[1];
    assert Length(outputResidueCode) eq n;
    assert Alphabet(outputResidueCode) eq GF(p);

    if checkDefinition then 
        expectedResidueCodewords := [Vp!v : v in C];
        expectedResidueCode := LinearCode(Matrix(expectedResidueCodewords));
        assert outputResidueCode eq expectedResidueCode;
    end if;  
    
    // Test for the ZpTorsionCode(C, i) function
    for i in [1..#type] do
        outputTorsionCode := ZpTorsionCode(C, i);
        
        // when i=1 check that the codes is equal to the ZpResidueCode(C) 
        if i eq 1 then 
            assert outputTorsionCode eq outputResidueCode;
        end if;

        // check that the codes form a chain of subcodes    
        if i ge 2 then 
            assert ZpTorsionCode(C, i-1) subset outputTorsionCode;
        end if;

        // check that it is a linear code over GF(2) of dimension t1+...+ti    
        assert Dimension(outputTorsionCode) eq &+type[1..i];
        assert Length(outputTorsionCode) eq n;
        assert Alphabet(outputTorsionCode) eq GF(p);

        // check the code according to the definition of a Torsion code
        // this test can take too much time depending on the size of the code
        if checkDefinition then 
            expectedTorsionCodewords := [ Vp![Z!w[j] div p^(i-1) : j in [1..n]] : 
                                           w in {v : v in C | p^(s-i+1)*v eq C!0 }];
            // another option that also works, but it is slower 
            //expectedTorsionCodewords := [ Vp!Vector(Z, c) : c in Generic(C) | (p^(i-1))*c in C ];
            expectedTorsionCode := LinearCode(Matrix(expectedTorsionCodewords));
            assert outputTorsionCode eq expectedTorsionCode;
        end if;  
    end for;
end procedure;

/****************************************************************/
print "test 1: Zero code over Z16 of length 10";

Zps := Integers(16);
n := 10;
C := ZeroCode(Zps, 10);  

checkDefinition := false;
TestZpResidueTorsionCode(C, checkDefinition);

/****************************************************************/
print "test 2: Universe code over Z8 of length 3";

Zps := Integers(8);
n := 3;
C := UniverseCode(Zps, n);

checkDefinition := true;
TestZpResidueTorsionCode(C, checkDefinition);

/****************************************************************/
print "test 3: Repetition code over Z25 of length 6";

Zps := Integers(25);
n := 6;
C := RepetitionCode(Zps, n);

checkDefinition := false;
TestZpResidueTorsionCode(C, checkDefinition);

/****************************************************************/
print "test 4: Linear code over Z9 of type (3; 1,2)";

C := LinearCode(Matrix(Integers(9), 3, 3, [1, 2, 2,
                                           0, 3, 0,
                                           0, 0, 3]));
checkDefinition := true;
TestZpResidueTorsionCode(C, checkDefinition);

/****************************************************************/
print "test 5: Linear code over Z25 of type (3; 0,1)";

C := LinearCode(Matrix(Integers(25), 1, 3, [5, 0, 0]));

checkDefinition := true;
TestZpResidueTorsionCode(C, checkDefinition);

/****************************************************************/
print "test 6: Linear code over Z343 of type (3; 0,0,2)";

C := LinearCode(Matrix(Integers(343), 2, 3, [ 49,   0,  49,
                                               0,  49, 245]));

checkDefinition := true;
TestZpResidueTorsionCode(C, checkDefinition);

/****************************************************************/
print "test 7: Linear code over Z27 of type (5; 2,0,0)";

C := LinearCode(Matrix(Integers(27), 2, 5, [ 1, 0,  8,  8, 20,
                                             0, 1, 26, 24,  6 ]));

checkDefinition := true;
TestZpResidueTorsionCode(C, checkDefinition);

/****************************************************************/
print "test 8: Linear code over Z243 of type (5; 2,0,0,0,0)";

C := LinearCode(Matrix(Integers(243), 2, 5, [1, 0, 212,  85,  82,
                                             0, 1,  85, 239, 122]));

checkDefinition := false;
TestZpResidueTorsionCode(C, checkDefinition);

/****************************************************************/
print "test 9: Linear code over Z125 of type (3; 0,1,0)";

C := LinearCode(Matrix(Integers(125), 1, 3, [5, 85, 65]));

checkDefinition := true;
TestZpResidueTorsionCode(C, checkDefinition);

/****************************************************************/
print "test 10: Linear code over Z25 of type (3; 2,1)";

C := LinearCode(Matrix(Integers(25), 3, 3, [1, 0, 599,
                                            0, 1, 560,
                                            0, 0,  5]));
checkDefinition := true;
TestZpResidueTorsionCode(C, checkDefinition);

/****************************************************************/
print "test 11: Linear code over Z27 of type (8; 3,3,2)";

// Generator matrix of Linear Code over Z3^3 in Non-Standard Form.
C := LinearCode(Matrix(Integers(27), 8, 8, [ 1, 0, 0, 2, 0, 6, 1, 3,
                                             0, 1, 1, 2, 0, 5, 2, 2,
                                             0, 0, 3, 0, 0, 3, 0, 3,
                                             0, 0, 0, 3, 0, 6, 0, 3,
                                             0, 0, 0, 0, 1, 6, 0, 4,
                                             0, 0, 0, 0, 0, 9, 0, 0,
                                             0, 0, 0, 0, 0, 0, 3, 6,
                                             0, 0, 0, 0, 0, 0, 0, 9]));
checkDefinition := false; // it is too big to check that satisfies definition
TestZpResidueTorsionCode(C, checkDefinition);

/****************************************************************/
print "test 12: Linear code over Z4 of type (6; 2,1)";

C := ZpHadamardCode(2, [3,1]);
checkDefinition := true;
TestZpResidueTorsionCode(C, checkDefinition);
assert ZpResidueCode(C) eq BinaryResidueCode(C);
assert ZpTorsionCode(C, 1) eq BinaryResidueCode(C);
assert ZpTorsionCode(C, 2) eq BinaryTorsionCode(C);
