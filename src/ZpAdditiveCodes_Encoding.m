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
/* File name: ZpAdditiveCodes_Encoding.m                     */
/*                                                           */
/* Comment: Package developed within the CCSG group          */
/*                                                           */
/* Authors: Adrián Torres-Martín and Mercè Villanueva        */
/*                                                           */
/* Revision version and last date: v1.0   10-02-2023         */
/*  (moved from ZpAdditiveCodes_Extension v2.0)              */
/*                                 v2.0   12-10-2023         */
/*                                                           */
/*************************************************************/
//Uncomment freeze when package finished
freeze

intrinsic ZpAdditiveCodes_Encoding_version() -> SeqEnum
{Return the current version of this package.}
    
    version := [2, 0];
    return version;

end intrinsic;

/****************************************************************
    GLOBAL VARIABLES
*****************************************************************/

declare verbose EncodingFlag, 1;

// Values used in function IsZpInformationSet to decide which method to use
// in order to check whether I is an information set for C_p or not
FACTOR_RESTRICTED := 7*10^(-5);
EXP_RESTRICTED := 1.15;
FACTOR_ORIGINAL := 2.45*10^(-5);
EXP_ORIGINAL := 1.00;

import "ZpAdditiveCodes_Core.m": ParyExpansion;
import "ZpAdditiveCodes_Core.m": ParyComposition;
import "ZpAdditiveCodes_Core.m": PhiZp;  
import "ZpAdditiveCodes_Core.m": IsLinearCodeOverZps;
import "ZpAdditiveCodes_Core.m": ZpsPhiInfo;
import "ZpAdditiveCodes_Core.m": ZpsPhiInverseInfo;
import "ZpAdditiveCodes_Core.m": MultiplyByP;
import "ZpAdditiveCodes_Core.m": DivideByP;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////          INFORMATION SPACES AND INFORMATION SETS                ////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/****************************************************************/
/*                                                              */
/* Function name: ReducePsi                                     */
/* Parameters:  u, p, s, t                                      */
/* Function description: Given a vector u over Z/p^s, return a  */
/*   vector over Z/p^t which is the image of u under function   */
/*   Psi_{s,t} coordinate-wise. For t=s, Psi_{s,s}(u)=u. For    */
/*   t=s-1, Psi_{s,s-t}=(u-i)/p where i=u mod p and 0<=i<p. In  */
/*   general, if [u_0,...,u_{s-1}] is the p-ary expansion of u, */
/*   then Psi_{s,t} is the element of Z/p^t such that its p-ary */
/*   expansion is [u_{s-t},...,u_{s-1}]. We have that           */
/*   Psi_{s,t}=Psi_{t+1,t}·Psi_{t+2,t+1}·...·Psi_{s,s-1}, where */
/*   · denotes composition of functions.                        */
/* Input parameters description:                                */
/*   - u : A vector over Z/p^s                                  */
/*   - p : A prime number                                       */
/*   - s : An integer greater than or equal to 2                */
/*   - t : An integer between 1 and s                           */
/* Output parameters description:                               */
/*   - A vector over Z/p^t, but given as a vector over Z/p^s    */
/*                                                              */ 
/* Function developed by Adrián Torres                          */
/*                                                              */
/****************************************************************/
ReducePsi := function(u, p, s, t)
    firstCoord := s-t+1;
    uOutput := [];
    for j in [1..OverDimension(u)] do
        uOutput := uOutput cat 
                   [ParyComposition(ParyExpansion(u[j])[firstCoord..s])];
    end for;
    return Vector(Integers(p^s), uOutput);
end function;

/****************************************************************/
/*                                                              */
/* Function name: SigmaTransform                                */
/* Parameters:  b, G, T, p, s                                   */
/* Function description: Given a vector b in the submodule      */
/*   (Z/p^s)^t1 x (Z/p^(s-1))^t2 x ··· x Zp^ts, return a        */
/*   vector over the same submodule such that the encoding is   */
/*   systematic (after multiplying by a generator matrix G in   */
/*   standard form and applying Carlet's Gray map).             */
/* Input parameters description:                                */
/*   - b : A vector of (Z/p^s)^t1 x (Z/p^(s-1))^t2 x ··· x Zp^ts*/
/*   - G : A generator matrix in standard form                  */
/*   - T : A sequence with the type t1,...,ts                   */
/*   - p : A prime number                                       */
/*   - s : An integer greater than or equal to 2                */
/* Output parameters description:                               */
/*   - A vector of (Z/p^s)^t1 x (Z/p^(s-1))^t2 x ··· x Zp^ts    */
/*                                                              */
/* Reference: "Systematic encoding and permutation decoding for */
/*   Z/p^s-linear codes", by A. Torres and M. Villanueva, IEEE  */
/*   Trans. Inf. Theory, 68(7), pp 4435-4443, 2022.             */
/*                                                              */
/* Function developed by Adrián Torres                          */
/*                                                              */
/****************************************************************/
SigmaTransform := function(b, G, T, p, s) 
    pos := 1;
    V := Parent(b);
    bs := [];
    for k in [1..s] do
        Vs := RSpace(Integers(p^s), T[k]);
        Vk := RSpace(Integers(p^(s-k+1)), T[k]);
        Vsum := RSpace(Integers(p^s), pos-1);

        prod := Vsum!bs * Submatrix(G, 1, pos, pos-1, T[k]);
        sum := Vk!Eltseq(b)[pos..(pos+T[k]-1)] - Vk!ReducePsi(prod, p, s, s-k+1);
        sumZps := Vs!Eltseq(sum);

        bs := bs cat Eltseq(sumZps);
        pos := pos + T[k]; 
    end for;
    return V!bs;
end function;

/****************************************************************/
/*                                                              */
/* Function name: InverseSigmaTransform                         */
/* Parameters:  b, G, T, p, s                                   */
/* Function description: Given a vector b in the submodule      */
/*   (Z/p^s)^t1 x (Z/p^(s-1))^t2 x ··· x Zp^ts, return a        */
/*   vector over the same submodule such that is the inverse    */
/*   of the SigmaTransform function.                            */
/* Input parameters description:                                */
/*   - b : A vector of (Z/p^s)^t1 x (Z/p^(s-1))^t2 x ··· x Zp^ts*/
/*   - G : A generator matrix in standard form                  */
/*   - T : A sequence with the type t1,...,ts                   */
/*   - p : A prime number                                       */
/*   - s : An integer greater than or equal to 2                */
/* Output parameters description:                               */
/*   - A vector of (Z/p^s)^t1 x (Z/p^(s-1))^t2 x ··· x Zp^ts    */
/*                                                              */
/* Reference: "Systematic encoding and permutation decoding for */
/*   Z/p^s-linear codes", by A. Torres and M. Villanueva, IEEE  */
/*   Trans. Inf. Theory, 68(7), pp 4435-4443, 2022.             */
/*                                                              */
/* Function developed by Adrián Torres                          */
/*                                                              */
/****************************************************************/
InverseSigmaTransform := function(b, G, T, p, s)
    pos := 1;
    V := Parent(b);
    bs := [];
    for k in [1..s] do
        Vs := RSpace(Integers(p^s), T[k]);
        Vk := RSpace(Integers(p^(s-k+1)), T[k]);
        Vsum := RSpace(Integers(p^s), pos-1);

        prod := Vsum!Eltseq(b)[1..(pos-1)] * Submatrix(G, 1, pos, pos-1, T[k]);
        sum := Vk!Eltseq(b)[pos..(pos+T[k]-1)] + Vk!ReducePsi(prod, p, s, s-k+1);
        sumZps := Vs!Eltseq(sum);

        bs := bs cat Eltseq(sumZps);
        pos := pos + T[k]; 
    end for;
    return V!bs;
end function;

/****************************************************************/
/*                                                              */
/* Function name: InformationFromCodeword                       */
/* Parameters:  y, Ip, encodingForward                          */
/* Function description: Given a vector over GF(p) of the       */
/*   ambient space of a systematic code C_p over GF(p), an      */
/*   information set Ip with respect to which C_p is systematic,*/
/*   and a systematic encoding function, return the information */
/*   associated to the codeword y if y belongs to C_p and an    */
/*   error message otherwise.                                   */
/* Input parameters description:                                */
/*   - A vector over GF(p)                                      */
/*   - A sequence with coordinate positions                     */
/*   - A systematic encoding function for code C_p              */ 
/* Output parameters description:                               */
/*   - A sequence with the coordinates of a vector in the       */
/*     information space of C_p if y belongs to C_p             */
/*                                                              */
/* Function developed by Adrián Torres                          */
/*                                                              */
/****************************************************************/
function InformationFromCodeword(y, Ip, encodingForward)
    infoY := Eltseq(y)[Ip];
    // check whether y belongs to the systematic code C_p, by encoding the information
    if y eq encodingForward(infoY) then
        return infoY;
    else
        error "Element is not in the codomain of the map where the inverse is defined";
    end if;
end function;

/******************************************************************************/
/*                                                                            */
/* Function name: ZpInformationSpace                                          */
/* Parameters:  C, IsSystematicEncoding                                       */
/* Function description: Given a linear code C over Z/p^s of type (n; t1,..., */
/*   ts), return V=(Z/p^s)^t1 x (Z/p^(s-1))^t2 x ··· x Zp^ts, that is, the    */
/*   space of information vectors for C. Note that V is a Z/p^s-submodule of  */
/*   (Z/p^s)^(t1+...+ts) whose first t1 coordinates of V are of order p^s,    */
/*   the next t2 coordinates of order p^(s-1), and so on, until the last ts   */
/*   coordinates of order p. The function also returns the space of           */
/*   information vectors for the corresponding code C_p=Phi_s(C) over GF(p),  */
/*   where Phi_s is Carlet's generalized Gray map; that is, the k-dimensional */
/*   vector space over GF(p), where k=(st1+(s-1)t2+...+ts). Finally, for the  */
/*   encoding process, it returns two isomorphisms f and f_p from these spaces*/
/*   of information vectors, V and V_p, onto C and C_p, respectively.         */
/*   The map f is given as a bijective map from V to C. Nevertheless, f_p is  */
/*   given as an injective map from V_p to GF(p)^(np^(s−1)), where np^(s−1) is*/
/*   the length of C_p, having the inverse only defined for the elements in   */
/*   C_p subseteq GF(p)^(np^(s−1)). These two maps are related to each other  */
/*   in the sense that Phi_s * f = f_p * Phi_s|Jp, where Phi_s|Jp denotes the */
/*   projection to the set of k coordinates Jp=Phi_s(J) subset [1, ...,       */
/*   (t1 +...+ ts)ps−1] as defined in below reference, for the set J=[1, ..., */
/*   t1 +...+ ts].                                                            */
/*   The parameter IsSystematicEncoding specifies whether the map f_p         */
/*   corresponds to a systematic encodings for the code C_p or not. It is set */
/*   to true by default.                                                      */
/* Input parameters description:                                              */
/*   - C : A linear code over Z/p^s                                           */
/* Output parameters description:                                             */
/*   - The information space for C                                            */
/*   - The information space for Phi_s(C)                                     */
/*   - An injective (bijective) map from V to C. Encoding for C               */
/*   - An injective map from Vp to C_p. Encoding for Phi_s(C)                 */
/*                                                                            */
/* Reference: "Systematic encoding and permutation decoding for Z/p^s-linear  */
/*   codes", by A. Torres and M. Villanueva, IEEE Trans. Inf. Theory, 68(7),  */
/*   pp 4435-4443, 2022.                                                      */
/*                                                                            */
/* Function developed by Adrián Torres                                        */
/*                                                                            */
/* Signature: (<CodeLinRng> C : <BoolElt> IsSystematicEncoding)               */
/*                             -> ModTupRng, ModTupFld, Map, Map              */
/*                                                                            */
/******************************************************************************/
intrinsic ZpInformationSpace(C::CodeLinRng : IsSystematicEncoding := true) 
                                        -> ModTupRng, ModTupFld, Map, Map 
{
Given a linear code C over Z/p^s of type (n; t1,..., ts), return V=(Z/p^s)^t1 
x (Z/p^(s-1))^t2 x ··· x Zp^ts, that is, the space of information vectors for C.
Note that V is a Z/p^s-submodule of (Z/p^s)^(t1+...+ts) whose first t1 coordinates 
of V are of order p^s, the next t2 coordinates of order p^(s-1), and so on, 
until the last ts coordinates of order p. The function also returns the space of
information vectors for the corresponding code C_p=Phi_s(C) over GF(p), where Phi_s
is Carlet's generalized Gray map; that is, the k-dimensional vector space over GF(p), 
where k=(st1+(s-1)t2+...+ts). Finally, for the encoding process, it returns two
isomorphisms f and f_p from these spaces of information vectors, V and V_p, onto
C and C_p, respectively. 

The map f is given as a bijective map from V to C. Nevertheless, f_p is given as 
an injective map from V_p to GF(p)^(np^(s−1)), where np^(s−1) is the length of C_p,
having the inverse only defined for the elements in C_p subseteq GF(p)^(np^(s−1)).
These two maps are related to each other in the sense that Phi_s * f = f_p * Phi_s|Jp, 
where Phi_s|Jp denotes the restriction to the set of k coordinates Jp=Phi_s(J)
subset [1, ..., (t1 +...+ ts)ps−1] as defined in the below reference, for the set
J=[1, ..., t1 +...+ ts]. 

The parameter IsSystematicEncoding specifies whether the map f_p corresponds to 
a systematic encodings for the code C_p or not. It is set to true by default. 
In this case, it returns a systematic encoding f_p with respect to the information
set I_p given by function ZpInformationSet(C). Indeed, f_p is such that the diagram
given in the manual commutes for the encoding f given in the below reference. 
Otherwise, it returns an encoding f_p, which may not be systematic. In particular, 
f_p is such that the diagram given in the manual commutes for the encoding f which 
corresponds to multiplying by the generator matrix given by function 
ZpMinRowsGeneratorMatrix(C).

If C is a linear code over Z4 of type (n; t1, t2), then InformationSpace(C) can also 
be used. The second output parameter coincides in both functions. However, 
InformationSpace(C) instead of returning Z4^t1 × pZ2^t2, it returns pZ2^t2 x Z4^t1, 
and both isomorphisms f and f_p from the spaces of information vectors, V and V_p, 
onto C and C_p, are also different.

Reference: "Systematic encoding and permutation decoding for Z/p^s-linear codes", 
by A. Torres and M. Villanueva, IEEE Trans. Inf. Theory, 68(7), pp 4435-4443, 2022. 
}  
    isOverZps, p, s := IsLinearCodeOverZps(C); 
    require isOverZps and (s ge 2): "The code C must be over Z/p^s with s>1";
    require Type(IsSystematicEncoding) eq BoolElt: 
                 "The optional parameter IsSystematicEncoding must be boolean";

    Zps := BaseRing(C);
    n := Length(C);
    np := p^(s-1)*n;
    V := VectorSpace(GF(p), np);
    
    // Special case for the zero code
    if PseudoDimension(C) eq 0 then
        InfRSpace := RSpace(Zps, 0);
        InfVSpace := VectorSpace(GF(p), 0); 
        encodingZps := map<InfRSpace -> C | x :-> C![0^^n],
                                            y :-> InfRSpace!0 >;
        encodingZp := map<InfVSpace -> V | x :-> V![0^^np],
                                           y :-> InfVSpace!0 >;

        return InfRSpace, InfVSpace, encodingZps, encodingZp;
    end if;

    T := ZpType(C);
    R := RSpace(Zps, &+T);

    diagonal := [1 : i in [1..T[1]]];
    for j in [2..#T] do
        diagonal := diagonal cat [p^(j-1) : i in [1..T[j]]];
    end for;
    InfCode := LinearCode(DiagonalMatrix(Zps, diagonal));
    InfRSpace := RSpace(InfCode);   
    InfVSpace := VectorSpace(GF(p), &+[T[i]*(s-i+1) : i in [1..#T]]);
    
    grayMap := CarletGrayMap(C);
    
    if IsSystematicEncoding then
        S, f, G := ZpStandardForm(C);
        encodingZps := map<InfRSpace -> C | 
            // first x is encoded using the standard form, and then it is sent to the code by f 
            x :-> (S!(SigmaTransform(DivideByP(x, p, T), G, T, p, s)*G)) @@ f,
            // first y is sent to the standard form by f, and then the information is obtained
            y :-> MultiplyByP(InverseSigmaTransform(Solution(G, f(y)), G, T, p, s), p, T) >;

        _, Ip := ZpInformationSet(C);
        encodingForward := func < x | grayMap((SigmaTransform(R!ZpsPhiInverseInfo(x, T), G, T, p, s)*G) @@ f) >;
        encodingZp := map<InfVSpace -> V | 
            // systematic encoding with respect to the information set Ip 
            x :-> encodingForward(x),
            // since it is systematic, it would be enough to return Eltseq(y)[Ip]
            // but since the codomain is the whole ambient space, instead of C_p, 
            // function InformationFromCodeword return an error message in case
            // y does not belong to C_p, or the information Eltseq(y)[Ip] otherwise    
            y :-> InfVSpace!InformationFromCodeword(y, Ip, encodingForward) >;
    else
        G := ZpMinRowsGeneratorMatrix(C);
        encodingZps :=  map<InfRSpace -> C | x :-> (C!(DivideByP(x, p, T)*G)),
                                             y :-> MultiplyByP(Solution(G, y), p, T) >;
        encodingZp := map<InfVSpace -> V | x :-> V!grayMap((R!ZpsPhiInverseInfo(x, T)*G)),
                                           y :-> ZpsPhiInfo(Solution(G, y @@ grayMap), T) >;
    end if;

    return InfRSpace, InfVSpace, encodingZps, encodingZp;

end intrinsic;

/******************************************************************************/
/*                                                                            */
/* Function name: ZpInformationSet                                            */
/* Parameters:  C                                                             */
/* Function description: Given a linear code C over Z/p^s of type (n;t1,...,  */
/*   ts), return an information set I subseteq [1,...,n] for C. Moreover, it  */
/*   also returns an information set I_p for the corresponding code C_p =     */
/*   Phi_s(C) over GF(p), where Phi_s is Carlet's generalized Gray map. These */
/*   information sets I and I_p are returned as a sequence of t1+t2+...+ts and*/
/*   st1+(s-1)t2+...+ts coordinate positions, respectively. The information   */
/*   set I_p coincides with Phi(I) as defined in the below reference, and the */
/*   encoding map f_p given by function ZpInformationSpace(C) is systematic   */ 
/*   with respect to the information set I_p.                                 */
/* Input parameters description:                                              */
/*   - C : A linear code over Z/p^s                                           */
/* Output parameters description:                                             */
/*   - An information set for C                                               */
/*   - An information set for Phi_s(C)                                        */
/*                                                                            */
/* Reference: "Systematic encoding and permutation decoding for Z/p^s-linear  */
/*   codes", by A. Torres and M. Villanueva, IEEE Trans. Inf. Theory, 68(7),  */
/*   pp 4435-4443, 2022.                                                      */
/*                                                                            */
/* Function developed by Adrián Torres                                        */
/*                                                                            */
/* Signature: (<CodeLinRng> C) -> SeqEnum[RngIntElt], SeqEnum[RngIntElt]      */
/*                                                                            */
/******************************************************************************/
intrinsic ZpInformationSet(C::CodeLinRng) -> SeqEnum[RngIntElt], SeqEnum[RngIntElt]
{
Given a linear code C over Z/p^s of type (n;t1,...,ts), return an information 
set I subseteq [1,...,n] for C. Moreover, it also returns an information set I_p 
for the corresponding code C_p=Phi_s(C) over GF(p), where Phi_s is Carlet's 
generalized Gray map. These information sets I and I_p are returned as a sequence
of t1+t2+...+ts and st1+(s-1)t2+...+ts coordinate positions, respectively. The
information set I_p coincides with Phi(I) as defined in the below reference, and 
the encoding map f_p given by function ZpInformationSpace(C) is systematic with 
respect to the information set I_p. 

An information set I for C is an ordered set of t1+t2+...+ts coordinate 
positions such that |C_I|=|C|, where C_I=[v_I : v in C] and v_I is the vector
v restricted to the I coordinates. An information set I_p for C_p is an ordered 
set of st1+(s-1)t2+...+ts coordinate positions such that |(C_p)_(I_p)|=|C_p|=|C|.

If C is a linear code over Z4, then function InformationSet(C) can also
be used even though both output parameters may be different.

Reference: "Systematic encoding and permutation decoding for Z/p^s-linear codes", 
by A. Torres and M. Villanueva, IEEE Trans. Inf. Theory, 68(7), pp 4435-4443, 2022. 
}
    isOverZps, p, s := IsLinearCodeOverZps(C); 
    require isOverZps and (s ge 2): "The code C must be over Z/p^s with s>1";

    n := Length(C);
    T := ZpType(C);
    _, _, _, perm := ZpStandardForm(C);

    Ip := [];
    blocks := RSpace(Integers(), n)![p^(s-1)*i + 1 : i in [0..(n-1)]];
    pos := blocks^(perm);
    coord := 1;
    for k in [1..#T] do
        for j in [1..T[k]] do
            J := [pos[coord] + i*p^(k-1) : i in [0..(T[k]*p^(s-k)-1)]];
            L := [1] cat [1 + p^(i-1) : i in [1..(s-k)]];
            Ip := Ip cat [J[i] : i in L];
            coord +:= 1;
        end for;
    end for;

    I := [i : i in [1..&+T]];

    return I^(perm^(-1)), Ip;

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: ProjectedCarletGrayMapImage                   */
/* Parameters: C, Ip                                            */
/* Function description: Given a linear code C over Z/p^s and an*/
/*   information set Ip for the code C_p=Phi_s(C), where Phi_s  */
/*   is the generalized Carlet's Gray map, return a sequence    */
/*   containing all codewords in Phi_s(C), projected to the     */
/*   coordinates in Ip. The projected codewords are given as a  */
/*   sequence of elements of GF(p).                             */ 
/*   Instead of performing the Gray map to the codewords, and   */
/*   then select the coordinates in Ip, the Gray map is only    */
/*   performed to the coordinates corresponding to Ip.          */
/* Input parameters description:                                */
/*   - C : A linear code over Z/p^s                             */
/*   - Ip : A sequence of integers in [1..np^(s-1)]             */
/* Output parameters description:                               */
/*   - A sequence of codewords (given as a sequence of elements */
/*     of GF(p)) of Phi_s(C) projected to the coordinates of Ip */
/*                                                              */
/* Function developed by Adrián Torres                          */
/*                                                              */
/****************************************************************/
function ProjectedCarletGrayMapImage(C, Ip)
    _, p, s := IsLinearCodeOverZps(C);

    // Coordinates in Z/p^s that correspond to each element in Ip after the Gray map
    coordinateZps := [(u-1) div p^(s-1) + 1 : u in Ip];
    // Columns of Y corresponding to each element in Ip
    columnsY := [(u-1) mod p^(s-1) + 1 : u in Ip];
    // The matrix used in the Gray Map
    Y := Matrix([x : x in VectorSpace(GF(p), s-1)]);
       
    return [&cat[Eltseq(PhiZp(c[coordinateZps[i]], Transpose(Matrix(Y[columnsY[i]])))) 
                    : i in [1..#Ip]] : c in C];
end function;

/******************************************************************************/
/*                                                                            */
/* Function name: IsZpInformationSet                                          */
/* Parameters: C, I                                                           */
/* Function description: Given a linear code C over Z/p^s of type (n;t1,...,  */
/*   ts) and a sequence I subseteq [1,...,n] or I subseteq [1,...,np^(s-1)],  */
/*   return true if and only if I subseteq [1,...,n] is an information set    */
/*   for C. This function also returns another boolean, which is true if an   */
/*   only if I subseteq [1,...,np^(s-1)] is an information set for the        */
/*   corresponding code C_p=Phi_s(C) over GF(p), where Phi_s is Carlet's      */
/*   generalized Gray map.                                                    */
/* Input parameters description:                                              */
/*   - C : A linear code over Z/p^s                                           */
/*   - I : A sequence of integers in [1..n] or [1..np^(s-1)]                  */
/* Output parameters description:                                             */
/*   - Boolean, true iff I is an information set for C                        */
/*   - Boolean, true iff I is an information set for C_p                      */
/*                                                                            */
/* Function initially developed by Roland Barrolleta and                      */
/*                    generalized by Adrián Torres                            */
/*                                                                            */
/* Signature: (<CodeLinRng> C, <[RngIntElt]> I) -> BoolElt, BoolElt           */
/*                                                                            */
/******************************************************************************/
intrinsic IsZpInformationSet(C::CodeLinRng, I::[RngIntElt]) -> BoolElt, BoolElt
{
Given a linear code C over Z/p^s of type (n;t1,t2,...,ts) and a sequence I subseteq
[1,...,n] or I subseteq [1,...,np^(s-1)], return true if and only if I subseteq 
[1,...,n] is an information set for C. This function also returns another boolean,
which is true if an only if I subseteq [1,...,np^(s-1)] is an information set for 
the corresponding code C_p=Phi_s(C) over GF(p), where Phi_s is Carlet's generalized 
Gray map.

An information set I for C is an ordered set of t1+t2+...+ts coordinate 
positions such that |C_I|=|C|, where C_I=[v_I : v in C] and v_I is the vector
v restricted to the I coordinates. An information set I_p for C_p is an ordered 
set of st1+(s-1)t2+...+ts coordinate positions such that |(C_p)^(I_p)|=|C_p|=|C|.

If C is over Z4, function IsZpInformationSet(C, I) coincides with function 
IsInformationSet(C, I), which works only for linear codes over Z4, but the former 
may perform less efficiently when I subset [1,...,2n].
}
    isOverZps, p, s := IsLinearCodeOverZps(C); 
    require isOverZps and (s ge 2): "The code C must be over Z/p^s with s>1";

    I := Set(I);  // Eliminate repeated coordinate positions in I
    k := #I;
    n := Length(C);
    np := p^(s-1)*n;
    require (I subset [1..n]) or (I subset [1..np]): 
          "Argument 2 should be a subset of", [1..n], "or a subset of", [1..np];

    T := ZpType(C);
    M := #C;
    sumT := &+T;

    // Special case for the zero code
    if (sumT eq 0) then
        if (k eq 0) then
            return true, true;
        else 
            return false, false;
        end if;
    end if;

    // Check whether I is an information set for C or not
    if (I subset [1..n]) and (sumT eq k) then
        checkSet := {1..n} diff I;
        isInfSetZps := M eq #PunctureCode(C, checkSet);
    else
        isInfSetZps := false;
    end if;
    
    // Check whether I is an information set for C_p or not
    if &+[T[i]*(s-i+1) : i in [1..s]] eq k then
        timeProjection := Log(p, FACTOR_RESTRICTED) + k*EXP_RESTRICTED;
        timeOriginal := Log(p, FACTOR_ORIGINAL) + Log(p, n) + k*EXP_ORIGINAL;
        if timeProjection lt timeOriginal then
            codewordsZpRest := Set(ProjectedCarletGrayMapImage(C, Sort(Setseq(I))));
        else
            codewordsZp := CarletGrayMapImage(C);
            codewordsZpRest := {[c[i] : i in I ] : c in codewordsZp};
        end if;
        isInfSetZp := M eq #codewordsZpRest;
    else
        isInfSetZp := false;
    end if;

    return isInfSetZps, isInfSetZp;

end intrinsic;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////                ENCODING FOR CODES OVER Z/p^s                    ////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/******************************************************************************/
/*                                                                            */
/* Function name: Encoding                                                    */
/* Parameters:  C, L                                                          */
/* Function description: Given a linear code C over Z/p^s of type (n; t1,t2,..*/
/*   ts) and a sequence L of elements from Z/p^s or GF(p), return a sequence  */
/*   of codewords from C, and also the corresponding sequence of codewords    */
/*   from C_p = Phi_s(C), given by an injective map, which corresponds to an  */
/*   encoding of the elements of L. The encodings for C and C_p coincide with */
/*   the ones provided by function ZpInformationSpace(C : IsSystematicEncoding*/
/*   := false).                                                               */
/*   Depending on the elements of L, the function automatically selects the   */
/*   appropriate encoding over Z/p^s or GF(p), by using function Encoding(C,  */
/*   v). If it detects that L contains more than one information vector, then */
/*   it computes the encoding for each one of them and returns a sequence of  */
/*   codewords. If it is necessary, zeros are added at the end of the         */
/*   sequence L to complete the last information vector before encoding.      */
/* Input parameters description:                                              */
/*   - C : A linear code over Z/p^s                                           */
/*   - L : A sequence of elements from Z/p^s or GF(p)                         */
/* Output parameters description:                                             */
/*   - A sequence of codewords of C                                           */
/*   - A sequence of codewords of Phi_s(C)                                    */
/*                                                                            */
/* Functions developed by Adrián Torres                                       */
/*                                                                            */
/* Signature: (<CodeLinRng> C, <SeqEnum[RngIntResElt]> L) -> [ModTupRngElt],  */
/*                                                           [ModTupFldElt]   */
/*            (<CodeLinRng> C, <SeqEnum[FldFinElt]> L) -> [ModTupRngElt],     */
/*                                                        [ModTupFldElt]      */
/*                                                                            */
/******************************************************************************/ 
intrinsic Encoding(C::CodeLinRng, L::SeqEnum[RngIntResElt]) -> SeqEnum, SeqEnum
{
Given a linear code C over Z/p^s of type (n;t1,t2,...,ts) and a sequence L
of elements from Z/p^s, return a sequence of codewords from C, and
also the corresponding sequence of codewords from C_p = Phi_s(C), given by
an injective map, which corresponds to an encoding of the elements of L.
The encodings for C and C_p coincide with the ones provided by function
ZpInformationSpace(C : IsSystematicEncoding := false). 

If it detects that L contains more than one information vector, then it
computes the encoding for each one of them and returns a sequence of
codewords. If it is necessary, zeros are added at the end of the sequence
L to complete the last information vector before encoding.
}
    isOverZps, p, s := IsLinearCodeOverZps(C);
    require isOverZps and (s ge 2): "The code C must be over Z/p^s with s>1";
    require PseudoDimension(C) gt 0: "The code C cannot be the zero code";

    V, Vp, f, fp := ZpInformationSpace(C : IsSystematicEncoding := false);

    Zps := Integers(p^s);
    Zp := Integers(p);
    Fp := GF(p);
    lengthL := #L;

    if lengthL eq 0 then
        return [], [];
    end if;

    if Universe(L) eq Zps then
        lengthV := OverDimension(V);
        // add zeros to complete the sequence having a multiple of lengthV
        L cat:=  [0^^((-lengthL) mod lengthV)];
        return Encoding(C, Matrix(Zps, lengthV, L));
    elif (Universe(L) eq Zp) then
        lengthVp := OverDimension(Vp);
        // add zeros to complete the sequence having a multiple of lengthVp
        L cat:=  [0^^((-lengthL) mod lengthVp)];
        return Encoding(C, Matrix(Fp, lengthVp, L));
    else
        error "The sequence must be defined over Z/p^s with p prime and s>=1.";
    end if;

end intrinsic;

intrinsic Encoding(C::CodeLinRng, L::SeqEnum[FldFinElt]) -> SeqEnum, SeqEnum
{
Given a linear code C over Z/p^s of type (n;t1,t2,...,ts) and a sequence L
of elements from GF(p), return a sequence of codewords from C, and
also the corresponding sequence of codewords from C_p = Phi_s(C), given by
an injective map, which corresponds to an encoding of the elements of L.
The encodings for C and C_p coincide with the ones provided by function
ZpInformationSpace(C : IsSystematicEncoding := false). 

If it detects that L contains more than one information vector, then it
computes the encoding for each one of them and returns a sequence of
codewords. If it is necessary, zeros are added at the end of the sequence
L to complete the last information vector before encoding.
}
    isOverZps, p, s := IsLinearCodeOverZps(C);
    require isOverZps and (s ge 2): "The code C must be over Z/p^s with s>1";
    require PseudoDimension(C) gt 0: "The code C cannot be the zero code";

    V, Vp, f, fp := ZpInformationSpace(C : IsSystematicEncoding := false);

    Fp := GF(p);
    lengthL := #L;

    if lengthL eq 0 then
        return [], [];
    end if;

    if Universe(L) eq Fp then
        lengthVp := OverDimension(Vp);
        // add zeros to complete the sequence having a multiple of lengthVp
        L cat:=  [0^^((-lengthL) mod lengthVp)];
        return Encoding(C, Matrix(Fp, lengthVp, L));
    else
        error "The sequence must be defined over a finite field GF(p) with p prime.";
    end if;

end intrinsic;

/******************************************************************************/
/*                                                                            */
/* Function name: Encoding                                                    */
/* Parameters:  C, M                                                          */
/* Function description: Given a linear code C over Z/p^s of type (n; t1,t2,..*/
/*   ts) and a matrix M over Z/p^s or GF(p), return a sequence of codewords   */
/*   from C, and also the corresponding sequence of codewords from C_p =      */
/*   Phi_s(C), given by an injective map, which corresponds to an encoding of */
/*   the rows of M. The encodings for C and C_p coincide with the ones        */
/*   provided by function ZpInformationSpace(C:IsSystematicEncoding := false).*/
/*   The matrix can contain either information vectors over Z/p^s or GF(p).   */
/*   Depending on the length of the rows of M and its entries, the function   */
/*   automatically selects the appropriate encoding over Z/p^s or GF(p), by   */
/*   using function Encoding(C, v).                                           */
/* Input parameters description:                                              */
/*   - C : A linear code over Z/p^s                                           */
/*   - M : A matrix with elements from Z/p^s or GF(p)                         */
/* Output parameters description:                                             */
/*   - A sequence of codewords of C                                           */
/*   - A sequence of codewords of Phi_s(C)                                    */
/*                                                                            */
/* Functions developed by Adrián Torres                                       */
/*                                                                            */
/* Signature: (<CodeLinRng> C, <ModMatRngElt> M) -> [ModTupRngElt],           */
/*                                                  [ModTupFldElt]            */
/*            (<CodeLinRng> C, <ModMatFldElt> M) -> [ModTupRngElt],           */
/*                                                  [ModTupFldElt]            */
/*                                                                            */
/******************************************************************************/
intrinsic Encoding(C::CodeLinRng, M::ModMatRngElt) -> SeqEnum, SeqEnum
{
Given a linear code C over Z/p^s of type (n;t1,...,ts) and a matrix M
over Z/p^s, returns a sequence of codewords from C, and also the
corresponding sequence of codewords from C_p = Phi_s(C), given by an
injective map, which corresponds to an encoding of the rows of M. The
encodings for C and C_p coincide with the ones provided by function
ZpInformationSpace(C : IsSystematicEncoding := false).
}
    isOverZps, p, s := IsLinearCodeOverZps(C);
    require isOverZps and (s ge 2): "The code C must be over Z/p^s with s>1";
    require PseudoDimension(C) gt 0: "The code C cannot be the zero code";
    
    V, Vp, f, fp := ZpInformationSpace(C : IsSystematicEncoding := false);

    nRows := Nrows(M);

    if nRows eq 0 then
        return [], [];
    end if;

    if (OverDimension(M[1]) eq OverDimension(V)) and (BaseRing(M[1]) eq Integers(p^s)) 
                                          and (M[1] in V) then
        codewords := [f(V!M[i]) : i in [1..nRows]];
        grayMap := CarletGrayMap(C);
        return codewords, [grayMap(codeword): codeword in codewords]; 
    elif (OverDimension(M[1]) eq OverDimension(Vp)) and (BaseRing(M[1]) eq Integers(p)) then
        return Encoding(C, ChangeRing(M, GF(p)));
    else
        error "The matrix must be over Z/p^s with p prime and s>=1.";
    end if;

end intrinsic;

/******************************************************************************/ 
intrinsic Encoding(C::CodeLinRng, M::ModMatFldElt) -> SeqEnum, SeqEnum
{
Given a linear code C over Z/p^s of type (n;t1,...,ts) and a matrix M
over GF(p), returns a sequence of codewords from C, and also the
corresponding sequence of codewords from C_p = Phi_s(C), given by an
injective map, which corresponds to an encoding of the rows of M. The
encodings for C and C_p coincide with the ones provided by function
ZpInformationSpace(C : IsSystematicEncoding := false).  
}
    isOverZps, p, s := IsLinearCodeOverZps(C);
    require isOverZps and (s ge 2): "The code C must be over Z/p^s with s>1";
    require PseudoDimension(C) gt 0: "The code C cannot be the zero code";

    V, Vp, f, fp := ZpInformationSpace(C : IsSystematicEncoding := false);

    nRows := Nrows(M);
    
    if nRows eq 0 then
        return [], [];
    end if;

    if (OverDimension(M[1]) eq OverDimension(Vp)) and (BaseRing(M[1]) eq GF(p)) then
        T := ZpType(C);
        return [f(V!MultiplyByP(ZpsPhiInverseInfo(M[i], T), p, T)) : i in [1..nRows]],
               [fp(Vp!M[i]) : i in [1..nRows]];
    else
        error "The matrix must be over a finite field GF(p) with p prime.";
    end if;

end intrinsic;

/******************************************************************************/
/*                                                                            */
/* Function name: Encoding                                                    */
/* Parameters:  C, v                                                          */
/* Function description: Given a linear code C over Z/p^s of type (n; t1,..., */
/*   ts) and an element v from the space Vp = GF(p)^k, where k=st1+(s−1)t2+...*/
/*   +ts, of information vectors for C_p = Phi_s(C), where Phi_s is Carlet’s  */
/*   generalized Gray map (or an element v from the space V=(Z/p^s)^t1 x      */
/*   (Z/p^(s-1))^t2 x ··· x Zp^ts of information vectors for C, given as a    */
/*   vector of length t1 +...+ ts over Z/p^s), return the codewords c in C and*/
/*   cp = Phi_s(c) in C_p corresponding to an encoding of (Phi_s|Jp)^(−1)(v)  */
/*   in V and v in Vp, respectively (or v in V and (Phi_s|Jp)(v) in Vp if v   */
/*   is given from V), where Phi_s|Jp is the projection of the image of Phi_s */
/*   onto the coordinates from Jp = Phi([1,..., t1+...+ts]).                  */ 
/*   This encoding for C, denoted by f, corresponds to multiplying an element */
/*   in V by the generator matrix given by function ZpMinRowsGeneratorMatrix  */
/*   (C) and the encoding for C_p corresponds to the map fp that makes diagram*/
/*   given in the manual commute for the encoding f. Note that                */
/*   f((Phi_s|Jp )^(−1)(v)) = c and fp(v) = cp (or f(v) = c and               */
/*   fp((Phi_s|Jp)(v)) = cp if v is given from V). The encodings f and fp     */
/*   coincide with the ones provided by function ZpInformationSpace(C :         */
/*   IsSystematicEncoding := false).                                          */
/* Input parameters description:                                              */
/*   - C : A linear code over Z/p^s                                           */
/*   - v : A vector over Z/p^s or GF(p)                                       */
/* Output parameters description:                                             */
/*   - A codeword of C                                                        */
/*   - A codeword of Phi_s(C)                                                 */
/*                                                                            */
/* Functions developed by Adrián Torres                                       */
/*                                                                            */
/* Signature: (<CodeLinRng> C, <ModTupRngElt> v) -> ModTupRngElt, ModTupFldElt*/
/*            (<CodeLinRng> C, <ModTupFldElt> v) -> ModTupRngElt, ModTupFldElt*/
/*                                                                            */
/******************************************************************************/
intrinsic Encoding(C::CodeLinRng, v::ModTupRngElt) -> ModTupRngElt, 
                                                      ModTupFldElt
{
Given a linear code C over Z/p^s of type (n; t1,...,ts) and an element v 
from the space V=(Z/p^s)^t1 x (Z/p^(s-1))^t2 x ··· x Zp^ts of information
vectors for C, given as a vector of length t1 +...+ ts over Z/p^s), return the
codewords c in C and cp = Phi_s(c) in C_p corresponding to an encoding of 
v in V and (Phi_s|Jp)(v) in Vp, where Phi_s|Jp is the projection of the image of 
Phi_s onto the coordinates from Jp = Phi([1,..., t1+...+ts]).

This encoding for C, denoted by f, corresponds to multiplying an element in V 
by the generator matrix given by function ZpMinRowsGeneratorMatrix(C), and the 
encoding for C_p corresponds to the map fp that makes diagram given in the 
manual commute for the encoding f. Note that f(v) = c and fp((Phi_s|Jp)(v)) = cp. 
The encodings f and fp coincide with the ones provided by function 
ZpInformationSpace(C : IsSystematicEncoding := false).
}
    isOverZps, p, s := IsLinearCodeOverZps(C);
    require isOverZps and (s ge 2): "The code C must be over Z/p^s with s>1";
    require PseudoDimension(C) gt 0: "The code C cannot be the zero code";

    V, Vp, f, fp := ZpInformationSpace(C : IsSystematicEncoding := false);

    if (OverDimension(v) eq OverDimension(V)) and (BaseRing(v) eq Integers(p^s)) then
        codeword := f(V!v); 
        return codeword, CarletGrayMap(C)(codeword);
    elif (OverDimension(v) eq OverDimension(Vp)) and (BaseRing(v) eq Integers(p)) then
        return Encoding(C, Vp!v);
    else
        error "The vector must be over Z/p^s with p prime and s>=1.";
    end if;

end intrinsic;

/******************************************************************************/ 
intrinsic Encoding(C::CodeLinRng, vp::ModTupFldElt) -> ModTupRngElt, 
                                                       ModTupFldElt
{
Given a linear code C over Z/p^s of type (n; t1,...,ts) and an element v from 
the space Vp = GF(p)^k, where k=st1+(s−1)t2+...+ts, of information vectors for 
C_p = Phi_s(C), where Phi_s is Carlet’s generalized Gray map, return the
codewords c in C and cp = Phi_s(c) in C_p corresponding to an encoding of
(Phi_s|Jp)^(−1)(v) in V and v in Vp, respectively, where Phi_s|Jp is the 
projection of the image of Phi_s onto the coordinates from 
Jp = Phi([1,..., t1+...+ts]).

This encoding for C, denoted by f, corresponds to multiplying an element in V 
by the generator matrix given by function ZpMinRowsGeneratorMatrix(C), and the 
encoding for C_p corresponds to the map fp that makes diagram given in the 
manual commute for the encoding f. Note that f((Phi_s|Jp )^(−1)(v)) = c and 
fp(v) = cp. The encodings f and fp coincide with the ones provided by function 
ZpInformationSpace(C : IsSystematicEncoding := false).
}
    isOverZps, p, s := IsLinearCodeOverZps(C);
    require isOverZps and (s ge 2): "The code C must be over Z/p^s with s>1";
    require PseudoDimension(C) gt 0: "The code C cannot be the zero code";

    V, Vp, f, fp := ZpInformationSpace(C : IsSystematicEncoding := false);

    if (OverDimension(vp) eq OverDimension(Vp)) and (BaseRing(vp) eq GF(p)) then
        T := ZpType(C);
        return f(V!MultiplyByP(ZpsPhiInverseInfo(vp, T), p, T)), fp(Vp!vp);
    else
        error "The vector must be over a finite field GF(p) with p prime.";
    end if;

end intrinsic;

/******************************************************************************/
/*                                                                            */
/* Function name: SystematicEncoding                                          */
/* Parameters:  C, L                                                          */
/* Function description: Given a linear code C over Z/p^s of type (n; t1,t2,..*/
/*   ts) and a sequence L of elements from Z/p^s or GF(p), return a sequence  */
/*   of codewords from C, and also the corresponding sequence of codewords    */
/*   from C_p = Phi_s(C). The encodings for C and C_p coincide with the ones  */
/*   provided by function ZpInformationSpace(C).Unlike function Encoding(C,L),*/
/*   in this case, the given encoding for C_p is systematic with respect to   */
/*   the information set Ip given by function ZpInformationSet(C).            */ 
/*   Depending on the elements of L, the function automatically selects the   */
/*   appropriate encoding over Z/p^s or GF(p), by using function Systematic   */
/*   Encoding(C, v). If it detects that L contains more than one information  */
/*   vector, then it computes the encoding for each one of them and returns a */
/*   sequence of  codewords. If it is necessary, zeros are added at the end of*/
/*   the sequence L to complete the last information vector before encoding.  */
/* Input parameters description:                                              */
/*   - C : A linear code over Z/p^s                                           */
/*   - L : A sequence of elements from Z/p^s or GF(p)                         */
/* Output parameters description:                                             */
/*   - A sequence of codewords of C                                           */
/*   - A sequence of codewords of Phi_s(C)                                    */
/*                                                                            */
/* Functions developed by Adrián Torres                                       */
/*                                                                            */
/* Signature: (<CodeLinRng> C, <SeqEnum[RngIntResElt]> L) -> [ModTupRngElt],  */
/*                                                           [ModTupFldElt]   */
/*            (<CodeLinRng> C, <SeqEnum[FldFinElt]> L) -> [ModTupRngElt],     */
/*                                                        [ModTupFldElt]      */
/*                                                                            */
/******************************************************************************/ 
intrinsic SystematicEncoding(C::CodeLinRng, L::SeqEnum[RngIntResElt]) -> SeqEnum, SeqEnum
{
Given a linear code C over Z/p^s of type (n;t1,t2,...,ts) and a sequence L
of elements from Z/p^s, return a sequence of codewords from C, and
also the corresponding sequence of codewords from C_p = Phi_s(C).
The encodings for C and C_p coincide with the ones provided by function
ZpInformationSpace(C). Unlike function Encoding(C, L), in this case, the given 
encoding for C_p is systematic with respect to the information set Ip given by 
function ZpInformationSet(C).

If it detects that L contains more than one information vector, then it
computes the encoding for each one of them and returns a sequence of
codewords. If it is necessary, zeros are added at the end of the sequence
L to complete the last information vector before encoding.
}  
    isOverZps, p, s := IsLinearCodeOverZps(C);
    require isOverZps and (s ge 2): "The code C must be over Z/p^s with s>1";
    require PseudoDimension(C) gt 0: "The code C cannot be the zero code";

    V, Vp, f, fp := ZpInformationSpace(C);

    Zps := Integers(p^s);
    Zp := Integers(p);
    Fp := GF(p);
    lengthL := #L;

    if lengthL eq 0 then
        return [], [];
    end if;

    if Universe(L) eq Zps then
        lengthV := OverDimension(V);
        // add zeros to complete the sequence having a multiple of lengthV
        L cat:=  [0^^((-lengthL) mod lengthV)];
        return SystematicEncoding(C, Matrix(Zps, lengthV, L));
    elif (Universe(L) eq Zp) then
        lengthVp := OverDimension(Vp);
        // add zeros to complete the sequence having a multiple of lengthVp
        L cat:=  [0^^((-lengthL) mod lengthVp)];
        return SystematicEncoding(C, Matrix(Fp, lengthVp, L));
    else
        error "The sequence must be defined over Z/p^s with p prime and s>=1.";
    end if;

end intrinsic;
 
intrinsic SystematicEncoding(C::CodeLinRng, L::SeqEnum[FldFinElt]) -> SeqEnum, SeqEnum
{
Given a linear code C over Z/p^s of type (n;t1,t2,...,ts) and a sequence L
of elements from GF(p), return a sequence of codewords from C, and
also the corresponding sequence of codewords from C_p = Phi_s(C).
The encodings for C and C_p coincide with the ones provided by function
ZpInformationSpace(C). Unlike function Encoding(C, L), in this case, the given 
encoding for C_p is systematic with respect to the information set Ip given by 
function ZpInformationSet(C).

If it detects that L contains more than one information vector, then it
computes the encoding for each one of them and returns a sequence of
codewords. If it is necessary, zeros are added at the end of the sequence
L to complete the last information vector before encoding.
}
    isOverZps, p, s := IsLinearCodeOverZps(C);
    require isOverZps and (s ge 2): "The code C must be over Z/p^s with s>1";
    require PseudoDimension(C) gt 0: "The code C cannot be the zero code";

    V, Vp, f, fp := ZpInformationSpace(C);

    Fp := GF(p);
    lengthL := #L;

    if lengthL eq 0 then
        return [], [];
    end if;

    if Universe(L) eq Fp then
        lengthVp := OverDimension(Vp);
        // add zeros to complete the sequence having a multiple of lengthVp
        L cat:=  [0^^((-lengthL) mod lengthVp)];
        return SystematicEncoding(C, Matrix(Fp, lengthVp, L));
    else
        error "The sequence must be defined over a finite field GF(p) with p prime.";
    end if;

end intrinsic;

/******************************************************************************/
/*                                                                            */
/* Function name: SystematicEncoding                                          */
/* Parameters:  C, M                                                          */
/* Function description: Given a linear code C over Z/p^s of type (n; t1,t2,..*/
/*   ts) and a matrix M over Z/p^s or GF(p), return a sequence of codewords   */
/*   from C, and also the corresponding sequence of codewords from C_p =      */
/*   Phi_s(C). The encodings for C and C_p coincide with the ones provided by */
/*   function ZpInformationSpace(C). Unlike function Encoding(C, M), in this  */
/*   case, the given encoding for C_p is systematic with respect to the       */
/*   information set Ip given by function ZpInformationSet(C).                */
/*   The matrix can contain either information vectors over Z/p^s or GF(p).   */
/*   Depending on the length of the rows of M and its entries, the function   */
/*   automatically selects the appropriate encoding over Z/p^s or GF(p), by   */
/*   using function SystematicEncoding(C, v).                                 */
/* Input parameters description:                                              */
/*   - C : A linear code over Z/p^s                                           */
/*   - M : A matrix with elements from Z/p^s or GF(p)                         */
/* Output parameters description:                                             */
/*   - A sequence of codewords of C                                           */
/*   - A sequence of codewords of Phi_s(C)                                    */
/*                                                                            */
/* Functions developed by Adrián Torres                                       */
/*                                                                            */
/* Signature: (<CodeLinRng> C, <ModMatRngElt> M) -> [ModTupRngElt],           */
/*                                                  [ModTupFldElt]            */
/*            (<CodeLinRng> C, <ModMatFldElt> M) -> [ModTupRngElt],           */
/*                                                  [ModTupFldElt]            */
/*                                                                            */
/******************************************************************************/
intrinsic SystematicEncoding(C::CodeLinRng, M::ModMatRngElt) -> SeqEnum, SeqEnum
{
Given a linear code C over Z/p^s of type (n;t1,...,ts) and a matrix M
over Z/p^s, returns a sequence of codewords from C, and also the
corresponding sequence of codewords from C_p = Phi_s(C). The
encodings for C and C_p coincide with the ones provided by function
ZpInformationSpace(C). Unlike function Encoding(C, M), in this case, the 
given encoding for C_p is systematic with respect to the information set Ip 
given by function ZpInformationSet(C).
}
    isOverZps, p, s := IsLinearCodeOverZps(C);
    require isOverZps and (s ge 2): "The code C must be over Z/p^s with s>1";
    require PseudoDimension(C) gt 0: "The code C cannot be the zero code";
    
    V, Vp, f, fp := ZpInformationSpace(C);

    nRows := Nrows(M);

    if nRows eq 0 then
        return [], [];
    end if;

    if (OverDimension(M[1]) eq OverDimension(V)) and (BaseRing(M[1]) eq Integers(p^s)) 
                                          and (M[1] in V) then
        codewords := [f(V!M[i]) : i in [1..nRows]];
        grayMap := CarletGrayMap(C);
        return codewords, [grayMap(codeword): codeword in codewords]; 
    elif (OverDimension(M[1]) eq OverDimension(Vp)) and (BaseRing(M[1]) eq Integers(p)) then
        return SystematicEncoding(C, ChangeRing(M, GF(p)));
    else
        error "The matrix must be over Z/p^s with p prime and s>=1.";
    end if;

end intrinsic;

/******************************************************************************/ 
intrinsic SystematicEncoding(C::CodeLinRng, M::ModMatFldElt) -> SeqEnum, SeqEnum
{
Given a linear code C over Z/p^s of type (n;t1,...,ts) and a matrix M
over GF(p), returns a sequence of codewords from C, and also the
corresponding sequence of codewords from C_p = Phi_s(C). The
encodings for C and C_p coincide with the ones provided by function
ZpInformationSpace(C). Unlike function Encoding(C, M), in this case, the 
given encoding for C_p is systematic with respect to the information set Ip 
given by function ZpInformationSet(C).
}
    isOverZps, p, s := IsLinearCodeOverZps(C);
    require isOverZps and (s ge 2): "The code C must be over Z/p^s with s>1";
    require PseudoDimension(C) gt 0: "The code C cannot be the zero code";

    V, Vp, f, fp := ZpInformationSpace(C);

    nRows := Nrows(M);
    
    if nRows eq 0 then
        return [], [];
    end if;

    if (OverDimension(M[1]) eq OverDimension(Vp)) and (BaseRing(M[1]) eq GF(p)) then
        T := ZpType(C);
        return [f(V!MultiplyByP(ZpsPhiInverseInfo(M[i], T), p, T)) : i in [1..nRows]],
               [fp(Vp!M[i]) : i in [1..nRows]];
    else
        error "The matrix must be over a finite field GF(p) with p prime.";
    end if;

end intrinsic;

/******************************************************************************/
/*                                                                            */
/* Function name: SystematicEncoding                                          */
/* Parameters:  C, v                                                          */
/* Function description: Given a linear code C over Z/p^s of type (n; t1,..., */
/*   ts) and an element v from the space Vp = GF(p)^k, where k=st1+(s−1)t2+...*/
/*   +ts, of information vectors for C_p = Phi_s(C), where Phi_s is Carlet’s  */
/*   generalized Gray map (or an element v from the space V=(Z/p^s)^t1 x      */
/*   (Z/p^(s-1))^t2 x ··· x Zp^ts of information vectors for C, given as a    */
/*   vector of length t1 +...+ ts over Z/p^s), return the codewords c in C and*/
/*   cp = Phi_s(c) in C_p corresponding to an encoding of (Phi_s|Jp)^(−1)(v)  */
/*   in V and v in Vp, respectively (or v in V and (Phi_s|Jp)(v) in Vp if v   */
/*   is given from V), where Phi_s|Jp is the projection of the image of Phi_s */
/*   onto the coordinates from Jp = Phi([1,..., t1+...+ts]). Unlike function  */
/*   Encoding(C, v), in this case, the given encoding for C_p is systematic   */
/*   with respect to the information set Ip given by function                 */
/*   ZpInformationSet(C).                                                     */
/*   This encoding for C_p, denoted by fp, corresponds to the systematic      */
/*   encoding with respect to the information set Ip given by function        */
/*   ZpInformationSet(C), and the encoding for C corresponds to the map f that*/
/*   makes diagram given in the manual commute for the encoding fp. Note that */
/*   f((Phi_s|Jp )^(−1)(v)) = c and fp(v) = cp (or f(v) = c and               */
/*   fp((Phi_s|Jp)(v)) = cp if v is given from V). Moreover, since fp is      */
/*   systematic, (cp)_Ip = v (or (cp)_Ip = (Phi_s|Jp)(v) if v is given from   */
/*   V). The encodings f and fp coincide with the ones provided by function   */
/*   ZpInformationSpace(C).                                                   */
/* Input parameters description:                                              */
/*   - C : A linear code over Z/p^s                                           */
/*   - v : A vector over Z/p^s or GF(p)                                       */
/* Output parameters description:                                             */
/*   - A codeword of C                                                        */
/*   - A codeword of Phi_s(C)                                                 */
/*                                                                            */
/* Functions developed by Adrián Torres                                       */
/*                                                                            */
/* Signature: (<CodeLinRng> C, <ModTupRngElt> v) -> ModTupRngElt, ModTupFldElt*/
/*            (<CodeLinRng> C, <ModTupFldElt> v) -> ModTupRngElt, ModTupFldElt*/
/*                                                                            */
/******************************************************************************/
intrinsic SystematicEncoding(C::CodeLinRng, v::ModTupRngElt) -> ModTupRngElt,
                                                                ModTupFldElt
{
Given a linear code C over Z/p^s of type (n; t1,...,ts) and an element v 
from the space V=(Z/p^s)^t1 x (Z/p^(s-1))^t2 x ··· x Zp^ts of information
vectors for C, given as a vector of length t1 +...+ ts over Z/p^s), return the
codewords c in C and cp = Phi_s(c) in C_p corresponding to an encoding of 
v in V and (Phi_s|Jp)(v) in Vp, where Phi_s|Jp is the projection of the image of 
Phi_s onto the coordinates from Jp = Phi([1,..., t1+...+ts]). Unlike function 
Encoding(C, v), in this case, the given encoding for C_p is systematic with 
respect to the information set Ip given by function ZpInformationSet(C).

This encoding for C_p, denoted by fp, corresponds to the systematic encoding with 
respect to the information set Ip given by function ZpInformationSet(C), and the 
encoding for C corresponds to the map f that makes diagram given in the manual 
commute for the encoding fp. Note that f(v) = c and fp((Phi_s|Jp)(v)) = cp. 
Moreover, since fp is systematic, (cp)_Ip = (Phi_s|Jp)(v). The encodings f and 
fp coincide with the ones provided by function ZpInformationSpace(C).
}
    isOverZps, p, s := IsLinearCodeOverZps(C);
    require isOverZps and (s ge 2): "The code C must be over Z/p^s with s>1";
    require PseudoDimension(C) gt 0: "The code C cannot be the zero code";

    V, Vp, f, fp := ZpInformationSpace(C);

    if (OverDimension(v) eq OverDimension(V)) and (BaseRing(v) eq Integers(p^s)) then
        codeword := f(V!v);
        return codeword, CarletGrayMap(C)(codeword);
    elif (OverDimension(v) eq OverDimension(Vp)) and (BaseRing(v) eq Integers(p)) then
        return SystematicEncoding(C, Vp!v);
    else
        error "The vector must be over Z/p^s with p prime and s>=1.";
    end if;

end intrinsic;

/******************************************************************************/ 
intrinsic SystematicEncoding(C::CodeLinRng, vp::ModTupFldElt) -> ModTupRngElt, 
                                                                 ModTupFldElt
{
Given a linear code C over Z/p^s of type (n; t1,...,ts) and an element v from 
the space Vp = GF(p)^k, where k=st1+(s−1)t2+...+ts, of information vectors for 
C_p = Phi_s(C), where Phi_s is Carlet’s generalized Gray map, return the
codewords c in C and cp = Phi_s(c) in C_p corresponding to an encoding of
(Phi_s|Jp)^(−1)(v) in V and v in Vp, respectively, where Phi_s|Jp is the 
projection of the image of Phi_s onto the coordinates from 
Jp = Phi([1,..., t1+...+ts]). Unlike function Encoding(C, v), in this case, 
the given encoding for C_p is systematic with respect to the information set 
Ip given by function ZpInformationSet(C).

This encoding for C_p, denoted by fp, corresponds to the systematic encoding with 
respect to the information set Ip given by function ZpInformationSet(C), and the 
encoding for C corresponds to the map f that makes diagram given in the manual 
commute for the encoding fp. Note that f(v) = c and fp((Phi_s|Jp)(v)) = cp. 
Moreover, since fp is systematic, (cp)_Ip = (Phi_s|Jp)(v). The encodings f and 
fp coincide with the ones provided by function ZpInformationSpace(C).
}
    isOverZps, p, s := IsLinearCodeOverZps(C);
    require isOverZps and (s ge 2): "The code C must be over Z/p^s with s>1";
    require PseudoDimension(C) gt 0: "The code C cannot be the zero code";

    V, Vp, f, fp := ZpInformationSpace(C);

    if (OverDimension(vp) eq OverDimension(Vp)) and (BaseRing(vp) eq GF(p)) then
        T := ZpType(C);
        return f(V!MultiplyByP(ZpsPhiInverseInfo(vp, T), p, T)), fp(Vp!vp);
    else
        error "The vector must be over a finite field GF(p) with p prime.";
    end if;

end intrinsic;
