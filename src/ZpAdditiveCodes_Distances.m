////////////////////////////////////////////////////////////////////////////////
/////////       Copyright 2021-2023 Adrián Torres and Mercè Villanueva  ////////
/////////                                                               ////////
/////////       This program is distributed under the terms of GNU      ////////
/////////               General Public License                          ////////
/////////                                                               ////////
////////////////////////////////////////////////////////////////////////////////

/*******************************************************************************
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful, 
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*******************************************************************************/


/*************************************************************/
/*                                                           */
/* Project name: Z/p^s-additive codes in MAGMA               */
/* File name: ZpAdditiveCodes_Distances.m                    */
/*                                                           */
/* Comment: Package developed within the CCSG group          */
/*                                                           */
/* Authors: Adrián Torres and Mercè Villanueva               */
/*                                                           */
/* Revision version and last date: v1.0   05-11-2021         */
/*                                                           */
/*************************************************************/
//Uncomment freeze when package finished
//freeze

intrinsic ZpAdditiveCodes_Distances_version() -> SeqEnum
{Return the current version of this package.}
    
    version := [1, 0];
    return version;

end intrinsic;

/****************************************************************
    GLOBAL VARIABLES
*****************************************************************/

import "ZpAdditiveCodes_Core.m": IsLinearCodeOverZps;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////          HOMOGENEOUS WEIGHT AND DISTRIBUTION                    ////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/******************************************************************************/
/*                                                                            */
/* Function name: HomogeneousWeight                                           */
/* Parameters:  a                                                             */
/* Function description: The homogeneous weight of the element a of Z/p^s     */
/*   with p prime and s>=2.                                                   */
/* Input parameters description:                                              */
/*   - a : An element of Z/p^s                                                */
/* Output parameters description:                                             */
/*   - The homogeneous weight of a                                            */
/*                                                                            */
/* Function developed by Adrián Torres                                        */
/*                                                                            */
/* Signature: (<RngIntResElt> a) -> RngIntElt                                 */
/*                                                                            */
/******************************************************************************/
intrinsic HomogeneousWeight(a::RngIntResElt) -> RngIntElt
{
The homogeneous weight of the element a of Z/p^s with p prime and s>=2.

If C is over Z4, this function coincides with function LeeWeight(a).
}
    ps := Factorization(#Parent(a))[1];
    p := ps[1];
    s := ps[2]; 
    require s ge 2 : "The element must be of Z/p^s with p prime and s>=2";
	require Parent(a) eq Integers(p^s) : 
	                 "The element must be of Z/p^s with p prime and s>=2";

    if a eq 0 then
        return 0;
    elif Integers()!a mod p^(s-1) eq 0 then
        return p^(s-1);
    else 
        return (p-1)*p^(s-2);
    end if;              
    
end intrinsic;

/******************************************************************************/
/*                                                                            */
/* Function name: HomogeneousWeight                                           */
/* Parameters:  v                                                             */
/* Function description: The homogeneous weight of the vector v in (Z/p^s)^n  */
/*   with p prime and s>=2.                                                   */
/* Input parameters description:                                              */
/*   - v : An element of (Z/p^s)^n                                            */
/* Output parameters description:                                             */
/*   - The homogeneous weight of v                                            */
/*                                                                            */
/* Function developed by Adrián Torres                                        */
/*                                                                            */
/* Signature: (<ModTupRngElt> v) -> RngIntElt                                 */
/*                                                                            */
/******************************************************************************/
intrinsic HomogeneousWeight(v::ModTupRngElt) -> RngIntElt
{
The homogeneous weight of the vector v in (Z/p^s)^n with p prime and s>=2.	

If C is over Z4, this function coincides with function LeeWeight(v).
}
    n := OverDimension(v);
    require n gt 0: "The length of the vector must be greater than 0";
    ps := Factorization(#Parent(v[1]))[1];
    p := ps[1];
    s := ps[2];
    require s ge 2 : "The vector must be over Z/p^s with p prime and s>=2";
	require Parent(v[1]) eq Integers(p^s) : 
	                 "The vector must be over Z/p^s with p prime and s>=2";

    w := 0;
    for i in [1..n] do
        w +:= HomogeneousWeight(v[i]);
    end for;

    return w;                
    
end intrinsic;

/******************************************************************************/
/*                                                                            */
/* Function name: HomogeneousDistance                                         */
/* Parameters:  u, v                                                          */
/* Function description: The homogeneous distance of the vectors u and v in   */
/*   (Z/p^s)^n with p prime and s>=2. This is defined to be the homogeneous   */
/*   weight of u-v.                                                           */
/* Input parameters description:                                              */
/*   - u : An element of (Z/p^s)^n                                            */
/*   - v : An element of (Z/p^s)^n                                            */
/* Output parameters description:                                             */
/*   - The homogeneous distance of u and v                                    */
/*                                                                            */
/* Function developed by Merce Villanueva                                     */
/*                                                                            */
/* Signature: (<ModTupRngElt> u, <ModTupRngElt> v) -> RngIntElt               */
/*                                                                            */
/******************************************************************************/
intrinsic HomogeneousDistance(u::ModTupRngElt, v::ModTupRngElt) -> RngIntElt
{
The homogeneous distance between the vectors u and v in (Z/p^s)^n with p prime and s>=2.
This is defined to be the homogeneous weight of u-v. 

If C is over Z4, this function coincides with function LeeDistance(u, v).
}
    n := OverDimension(v);
    require n gt 0: "The length of the vectors must be greater than 0";
	require n eq OverDimension(u): "Both vectors must be of the same length";
    ps := Factorization(#Parent(v[1]))[1];
    p := ps[1];
    s := ps[2];
    require s ge 2 : "The vector must be over Z/p^s with p prime and s>=2";
	require Parent(u[1]) eq Integers(p^s) : 
	                 "The first vector must be over Z/p^s with p prime and s>=2";
	require Parent(v[1]) eq Integers(p^s) : 
	                 "The second vector must be over Z/p^s with p prime and s>=2";
	
    return HomogeneousWeight(u-v);

end intrinsic;

/******************************************************************************/
/*                                                                            */
/* Function name: MinimumHomogeneousWeight                                    */
/* Parameters:  C                                                             */
/* Function description: The minimum homogeneous weight of the linear code C  */
/*   over Z/p^s with p prime and s>=2.                                        */
/* Input parameters description:                                              */
/*   - C : A linear code over Z/p^s                                           */
/* Output parameters description:                                             */
/*   - The minimum homogeneous weight of C                                    */
/*                                                                            */
/* Function developed by Adrián Torres                                        */
/*                                                                            */
/* Signature: (<CodeLinRng> C) -> RngIntElt                                   */
/*                                                                            */
/******************************************************************************/
intrinsic MinimumHomogeneousWeight(C::CodeLinRng) -> RngIntElt
{
The minimum homogeneous weight of the linear code C over Z/p^s with p prime and s>=2.

If C is over Z4, this function coincides with function MinumumLeeWeight(C) and 
MinumumLeeDistance(C), but the former may perform less efficiently because they 
are implemented just by using a brute force approach.
}
    isOverZps, p, s := IsLinearCodeOverZps(C); 
    require isOverZps and (s ge 1): "The code C must be over Z/p^s with s>=2";

    minWeight := #BaseRing(C)*Length(C);
    for c in C do
        newWeight := HomogeneousWeight(c);
        if (newWeight lt minWeight) and (newWeight ne 0) then 
            minWeight := newWeight; 
        end if;
    end for;

    return minWeight;               
    
end intrinsic;

intrinsic MinimumHomogeneousDistance(C::CodeLinRng) -> RngIntElt 
{
The minimum homogeneous weight of the linear code C over Z/p^s with p prime and s>=2.

If C is over Z4, this function coincides with functions MinumumLeeWeight(C) and 
MinumumLeeDistance(C), but the former may perform less efficiently because they 
are implemented just by using a brute force approach.
}
    return MinimumHomogeneousWeight(C);

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: AddWeightDistribution                         */
/* Parameters:  weightSeq, weightDistribution                   */
/* Function description: Given a sequence which has the number  */
/*   of codewords of each Lee weight, update the sequence by    */
/*   adding the values given  by weightDistribution             */  
/* Input parameters description:                                */
/*   - weightSeq : Sequence of integers                         */
/*   - weightDistribution : Sequence of tuples.                 */   
/*                                                              */
/****************************************************************/ 
procedure AddWeightDistribution(~weightSeq, weightDistribution)  
    for tuple in weightDistribution do 
        weightSeq[tuple[1] + 1] := weightSeq[tuple[1] + 1] + tuple[2]; 
    end for; 
end procedure; 

/******************************************************************************/
/*                                                                            */
/* Function name: HomogeneousWeightDistribution                               */
/* Parameters: C                                                              */
/* Function description: The homogeneous weight distribution of the linear    */
/*   code C over Z/p^s with p prime and s>=2. The distribution is returned in */
/*   the form of a sequence of tuples, where the ith tuple contains the ith   */
/*   homogeneous weight, wi say, and the number of codewords having weight wi.*/
/* Input parameters description:                                              */
/*   - C: A linear code over Z/p^s                                            */
/* Output parameters description:                                             */
/*   - Sequence of tuples <Lee weight, number of codewords>                   */
/*                                                                            */
/* Signature: (<CodeLinRng> C) -> SeqEnum                                     */
/*                                                                            */
/******************************************************************************/
intrinsic HomogeneousWeightDistribution(C::CodeLinRng) -> SeqEnum
{
The homogeneous weight distribution of the linear code C over Z/p^s with p prime and s>=2. 
The distribution is returned in the form of a sequence of tuples, where the ith tuple 
contains the ith homogeneous weight, wi say, and the number of codewords having weight wi.

This function coincides with function LeeWeightDistribution(C) for linear codes C over Z4.
}
	isOverZps, p, s := IsLinearCodeOverZps(C); 
    require isOverZps and (s ge 1): "The code C must be over Z/p^s with s>=2";

    if (PseudoDimension(C) eq 0) then 
        return [<0, 1>];
    else
        n := Length(C)*p^(s-1) + 1; 
        weightHomSeq := [0^^n];

        for c in C do 
            weightHomCodeword := HomogeneousWeight(c); 
            weightHomSeq[weightHomCodeword + 1] +:= 1;
        end for; 

        return [<i-1, weightHomSeq[i]> : i in [1..n] | not IsZero(weightHomSeq[i])];
    end if;

end intrinsic;

/******************************************************************************/
/*                                                                            */
/* Function name: DualHomogeneousWeightDistribution                           */
/* Parameters: C                                                              */
/* Function description: The homogeneous weight distribution of the dual of   */
/*   the linear code C over Z/p^s with p prime and s>=2                       */
/*   (see HomogeneousWeightDistribution). The distribution is returned in     */
/*   the form of a sequence of tuples, where the ith tuple contains the ith   */
/*   homogeneous weight, wi say, and the number of codewords having weight wi.*/
/* Input parameters description:                                              */
/*   - C: A linear code over Z/p^s                                            */
/* Output parameters description:                                             */
/*   - Sequence of tuples <Lee weight, number of codewords>                   */
/*                                                                            */
/* Signature: (<CodeLinRng> C) -> SeqEnum                                     */
/*                                                                            */
/******************************************************************************/
intrinsic DualHomogeneousWeightDistribution(C::CodeLinRng) -> SeqEnum
{
The homogeneous weight distribution of the dual of the linear code C over Z/p^s with p prime 
and s>=2 (see HomogeneousWeightDistribution). The distribution is returned in the form of 
a sequence of tuples, where the i-th tuple contains the i-th homogeneous weight, wi say, 
and the number of codewords having weight wi. 

This function coincides with function DualLeeWeightDistribution(C) for linear codes C over Z4.
}
    isOverZps, p, s := IsLinearCodeOverZps(C); 
    require isOverZps and (s ge 1): "The code C must be over Z/p^s with s>=2";

    //if (ZpInformationRate(C) gt 0.5) then
        return HomogeneousWeightDistribution(ZpDual(C));
    //else
    //    type := ZpType(C);
    //    dimension := PseudoDimension(C);
    //    pLength := Length(C)*s;     
    //    W := HomogeneousWeightDistribution(C);
    //    return MacWilliamsTransform(pLength, dimension, p, W);
    //end if;

end intrinsic;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////                LEE WEIGHT AND DISTRIBUTION                      ////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/******************************************************************************/
/*                                                                            */
/* Function name: LeeWeight                                                   */
/* Parameters:  a                                                             */
/* Function description: The generalized Lee weight of the element a in Z/p^s */
/*   with p prime and s>=2.                                                   */
/* Input parameters description:                                              */
/*   - a : An element of Z/p^s                                                */
/* Output parameters description:                                             */
/*   - The Lee weight of a                                                    */
/*                                                                            */
/* Function developed by Adrián Torres                                        */
/*                                                                            */
/* Signature: (<RngIntResElt> a) -> RngIntElt                                 */
/*                                                                            */
/******************************************************************************/
intrinsic LeeWeight(a::RngIntResElt) -> RngIntElt
{
The generalized Lee weight of the element a in Z/p^s with p prime and s>=2.
}
    ps := Factorization(#Parent(a))[1];
    p := ps[1];
    s := ps[2]; 
    require s ge 2 : "The element must be of Z/p^s with p prime and s>=2";
	require Parent(a) eq Integers(p^s) : 
	                 "The element must be of Z/p^s with p prime and s>=2";
    
    try
        Z := Integers();
        w := Min([Z!a, Z!(-a)]);
    catch e
        error "The element must be in Z/p^s with p prime and s>=1";
    end try;

    return w;                
    
end intrinsic;

/******************************************************************************/
/*                                                                            */
/* Function name: LeeWeight                                                   */
/* Parameters:  v                                                             */
/* Function description: The generalized Lee weight of the vector v in        */
/*   (Z/p^s)^n.                                                               */
/* Input parameters description:                                              */
/*   - v : An element of (Z/p^s)^n                                            */
/* Output parameters description:                                             */
/*   - The Lee weight of v                                                    */
/*                                                                            */
/* Function developed by Adrián Torres                                        */
/*                                                                            */
/* Signature: (<ModTupRngElt> v) -> RngIntElt                                 */
/*                                                                            */
/******************************************************************************/
intrinsic LeeWeight(v::ModTupRngElt) -> RngIntElt
{
The generalized Lee weight of the vector v in (Z/p^s)^n with p prime and s>=2.	
}
	n := OverDimension(v);
    require n gt 0: "The length of the vector must be greater than 0";
    ps := Factorization(#Parent(v[1]))[1];
    p := ps[1];
    s := ps[2];
    require s ge 2 : "The vector must be over Z/p^s with p prime and s>=2";
	require Parent(v[1]) eq Integers(p^s) : 
	                 "The vector must be over Z/p^s with p prime and s>=2";
    
	w := 0;
    for i in [1..n] do
        try
            w +:= LeeWeight(v[i]);
        catch e
            error "The vector must be an element of (Z/p^s)^n";
        end try;
    end for;

    return w;                
    
end intrinsic;