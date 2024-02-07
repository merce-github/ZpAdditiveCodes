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
/* File name: ZpAdditiveCodes_Decoding.m                     */
/*                                                           */
/* Comment: Package developed within the CCSG group          */
/*                                                           */
/* Authors: Adrián Torres-Martín and Mercè Villanueva        */
/*                                                           */
/* Revision version and last date: v1.0   07-07-2021         */
/*                                 v2.0   10-02-2023         */
/*                                 v2.1   13-06-2023         */
/*                                 v2.2   22-11-2023         */
/*                                                           */
/*************************************************************/
//Uncomment freeze when package finished
//freeze

intrinsic ZpAdditiveCodes_Decoding_version() -> SeqEnum
{Return the current version of this package.}
    
    version := [2, 2];
    return version;

end intrinsic;

/****************************************************************
    GLOBAL VARIABLES
*****************************************************************/

declare verbose IsPDSetFlag, 3;
declare verbose PDSetHadamardFlag, 2;

Z := Integers();

// Maximum time to find randomly an r-PD-set for the Hadamard code over Z/p^s
MAXTIME := 5.0;
// Maximum number of iterations to try to find randomly an r-PD-set for  
// the Hadamard code over Z/p^s 
NUMITER := 5;

import "ZpAdditiveCodes_Core.m": IsLinearCodeOverZps;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////				PERMUTATION DECODE FUNCTIONS                        ///////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/****************************************************************/
/*                                                              */  
/* Function name: ExistPermutationInPDset                       */
/* Parameters: S, I, v                                          */
/* Function description: Given a PD-set S, an information set I */
/*   and a vector v, return whether there exists a              */
/*   permutation in S such that moves the nonzero coordinates   */
/*   of v out of the information coordinates.                   */  
/*   The PD-set S and information set I are given as sequences. */
/* Input parameters description:                                */
/*   - S : A sequence of elements in Sym(n)                     */
/*   - I : A sequence with coordinate positions                 */
/*   - v : A vector over GF(p)                                  */ 
/* Output parameters description:                               */
/*   - Boolean, true iff there exists a permutation in S such   */
/*     that moves the nonzero coordinats of v out of I          */
/*                                                              */
/* Function developed by Roland Barrolleta                      */
/*                                                              */
/****************************************************************/
ExistPermutationInPDset := function (S, I, v)
    for p in S do
        if IsZero([(v^p)[i] : i in I]) then
            return true;
        end if;
    end for;
    return false;
end function;

/****************************************************************/
/*                                                              */  
/* Function name: IsSubsetOfPAut                                */
/* Parameters: C, S                                             */
/* Function description: Given a linear code C over Z/p^s of    */
/*   length n or a set C with the corresponding codewords, and  */
/*   a sequence S of permutations in Sym(n) or Sym(np^(s-1)),   */
/*   return true whether S is a subset of PAut(C).              */ 
/* Input parameters description:                                */
/*   - C : A linear code over Z/p^s or a set of codewords       */
/*   - S : A sequence of elements in Sym(n)                     */
/* Output parameters description:                               */
/*   - Boolean, true iff S is a subset of PAut(C)               */
/*                                                              */
/* Function developed by Roland Barrolleta                      */
/*                                                              */
/****************************************************************/
IsSubsetOfPAut := function (C, S)
    for p in S do
        if C^p ne C then
            return false;
        end if;
    end for;
    return true;
end function;

/*******************************************************************************/
/*                                                                             */
/* Function name: IsZpPermutationDecodeSet                                     */
/* Parameters: C, I, S, r                                                      */
/* Function description: Given a linear code C over Z/p^s of type (n; t1,...,  */
/*   ts), a sequence I in [1,...,np^(s-1)], a sequence S of elements in the    */
/*   symmetric group Sym(np^(s-1)) of permutations on the set [1,...,np^(s-1)],*/
/*   and an integer r>=1, return true if and only if S is a r-PD-set for       */
/*   Cp=Phi_s(C), where Phi_s is Carlet's generalized Gray map, with respect to*/
/*   the information set I. The parameters I and S can also be given as a      */
/*   sequence I in [1,...,n] and a sequence S of elements in the symmetric     */
/*   group Sym(n) of permutations on the set [1,...,n], respectively. In this  */
/*   case, the function returns true if and only if Phi(S) is a r-PD-set for   */
/*   Cp=Phi_s(C) with respect to the information set Phi(I), where Phi(I) and  */
/*   Phi(S) are the sequences defined as in the manual.                        */
/* Input parameters description:                                               */
/*   - C : A linear code over Z/p^s                                            */
/*   - I : A subset of integers in {1..np^(s-1)} or {1..n} as a sequence       */
/*   - S : A subset of permutations acting on {1..n} or {1..np^(s-1)},         */
/*         respectively, as a sequence                                         */    
/*   - r : An integer in {1..t}, where t is the correcting capability          */
/* Output parameters description:                                              */
/*   - Boolean, true if S is an r-PD-set and false otherwise                   */
/*                                                                             */
/* Function initially developed by Roland Barrolleta and                       */
/*                    generalized by Adrián Torres                             */
/*                                                                             */
/* Signature: (<CodeLinRng> C, <[RngIntElt]> I,                                */
/*              <[GrpPermElt]> S, <RngIntElt> r) -> BoolElt                    */
/*                                                                             */
/*******************************************************************************/
intrinsic IsZpPermutationDecodeSet(C::CodeLinRng, I::[RngIntElt], S::[GrpPermElt], 
                                                        r::RngIntElt) -> BoolElt
{
Given a linear code C over Z/p^s of type (n; t1,...,ts), a sequence I in 
[1,...,np^(s-1)], a sequence S of elements in the symmetric group Sym(np^(s-1)) of 
permutations on the set [1,...,np^(s-1)], and an integer r>=1, return true if and 
only if S is a r-PD-set for Cp=Phi_s(C), where Phi_s is Carlet's generalized Gray 
map, with respect to the information set I.                                      

The parameters I and S can also be given as a sequence I in [1,...,n] and a 
sequence S of elements in the symmetric group Sym(n) of permutations on the set 
[1,...,n], respectively. In this case, the function returns true if and only if 
Phi(S) is a r-PD-set for Cp=Phi_s(C) with respect to the information set Phi(I), 
where Phi(I) and Phi(S) are the sequences defined as in the manual.

Depending on the length of the code C, its type, and the integer r, this function 
could take some time to compute whether S or Phi(S) is a r-PD-set for Cp with 
respect to I or Phi(I), respectively. Specifically, if the function returns true, 
it is necessary to check sum_(i=1)^r (|I| choose i)*((N-|I|) choose (r-i)) r-sets, 
where N=n and |I|=t1+···+ts when I is given as an information set for C, or 
N=np^(s-1) and |I|=s*t1+···+ts when I is given as an information set for Cp=Phi_s(C).

The verbose flag IsPDsetFlag is set to level 0 by default. If it is set to level 
1, the total time used to check the condition is shown. Moreover, the reason why 
the function returns false is also shown, that is, whether I is not an information 
set, S is not a subset of the permutation automorphism group or S is not an r-PD-set. 
If it is set to level 2, the percentage of the computation process performed is also 
printed.

If C is over Z4, function IsZpPermutationDecodeSet(C, I, S, r) coincides with function 
IsPermutationDecodeSet(C, I, S, r), which works only for linear codes over Z4, but the 
former may perform less efficiently when I subset [1,...,2n] because it calls function 
IsZpInformationSet(C, I) instead of IsInformationSet(C, I).
}   
    isOverZps, p, s := IsLinearCodeOverZps(C);
    require isOverZps and (s ge 2): "The code C must be over Z/p^s with s>1"; 
    require PseudoDimension(C) gt 0: "The code C cannot be the zero code";
    requirege r, 1;

    require not IsEmpty(S): "Argument 3 cannot be an empty sequence";   
    n := Length(C);
	np := p^(s-1)*n; 
    numSymG := Degree(Parent(S[1]));
    require ((numSymG eq n) or (numSymG eq np)): 
          "Argument 3 should contain permutations acting on a set of cardinality", 
                                                                    n, "or", np;
    k := #I;
    require (k ge 1): "Argument 2 cannot be an empty sequence";
    if numSymG eq n then
        require (Min(I) ge 1) and (Max(I) le n): 
                                     "Argument 2 should be a subset of", [1..n];
    else  // numSymG = np
        require (Min(I) ge 1) and (Max(I) le np): 
                                  "Argument 2 should be a subset of", [1..np];
    end if;
    
    ///////////////////////////////////////////////////////////////////
	iniTime := Cputime();
	vprintf IsPDSetFlag, 2: "Checking whether I is an information set...\n";
    ///////////////////////////////////////////////////////////////////
    
	isInfSetZps, isInfSetZp := IsZpInformationSet(C, I);

    I := Set(I);  // Eliminate repeated coordinate positions in I
    S := Set(S);  // Eliminate repeated permutations in S

    ///////////////////////////////////////////////////////////////////
	vprintf IsPDSetFlag, 2: 
                  "Checking whether S is in the permutation automorphism group...\n";
    ///////////////////////////////////////////////////////////////////
     
    case numSymG:
        when n:  
            if not isInfSetZps then
                vprintf IsPDSetFlag, 1: "Argument 2 is not an information set for C.\n";
                return false; 
            end if;
            if not IsSubsetOfPAut(C, S) then 
                vprintf IsPDSetFlag, 1: 
                  "Argument 3 is not a subset of the permutation automorphism group of C.\n";   
                return false;     
            end if;
            length := n;

        when np:
            if not isInfSetZp then
                vprintf IsPDSetFlag, 1: "Argument 2 is not an information set for Cp = Phi(C).\n";
                return false; 
            end if;
            if not IsSubsetOfPAut(Set(CarletGrayMapImage(C)), S) then   
                vprintf IsPDSetFlag, 1: 
                  "Argument 3 is not a subset of the permutation automorphism group of Cp = Phi(C).\n";  
                return false;      
            end if;
            length := np;
    end case;
     
    ///////////////////////////////////////////////////////////////////
	vprintf IsPDSetFlag, 2: "Checking whether S is an r-PD-set...\n";
    ///////////////////////////////////////////////////////////////////

    numCheckSets := &+[ Binomial(k, i) * Binomial(length-k, r-i) : i in [1..r]];
    
    ///////////////////////////////////////////////////////////////////
	tenpc := numCheckSets div 10 + 1; //for verbose flag
    ///////////////////////////////////////////////////////////////////
    
    cont := 0;
    V := VectorSpace(GF(2), length);
    checkSet := {1..length} diff I;
    for numErrors in [1..r] do
        allSetsErrorsInfo := Subsets(I, numErrors);
        for errorsSetInfo in allSetsErrorsInfo do
            allSetsErrorsCheck := Subsets(checkSet, r-numErrors);
            for errorsSetCheck in allSetsErrorsCheck do
                errorSet := errorsSetInfo join errorsSetCheck;
                
                errorVec := V!0;
                for i in errorSet do
                    errorVec[i] := 1;
                end for;
                
                cont := cont + 1;
                //////////////////////////////////////////////////////////
                if cont mod tenpc eq 0 then
                    vprintf IsPDSetFlag, 2: "%o %%\n",(cont div tenpc * 10);
                end if;
                //////////////////////////////////////////////////////////
                
                if not ExistPermutationInPDset(S, I, errorVec) then
                
                    /////////////////////////////////////////////////////
                    vprintf IsPDSetFlag, 1: "Argument 3 is not an r-PD-set.\n";
	                vprintf IsPDSetFlag, 1: "Took %o seconds (CPU time).\n", Cputime(iniTime);
                    /////////////////////////////////////////////////////
                
                    return false;
                    
                end if;
            end for;
        end for;
    end for;

    /////////////////////////////////////////////////////
	vprintf IsPDSetFlag, 1: "Took %o seconds (CPU time).\n", Cputime(iniTime);
    /////////////////////////////////////////////////////

	return true;
	 
end intrinsic;

/****************************************************************/
/*                                                              */  
/* Function name: PermZpsToPermZp                               */
/* Parameters: permZps                                          */
/* Function description: Given permutation permZps from Sym(n), */
/*   return the permutation in Sym(np^(s-1)) such that          */
/*   i -> p^(s-1)*permZps((i+x(i))/p^(s-1))-x(i),               */
/*   where x(i)=p^(s-1)-(i mod p^(s-1)).                        */
/* Input parameters description:                                */
/*   - permZps : A permutation in Sym(n)                        */
/*   - p : A prime number                                       */
/*   - s : A positive integer                                   */ 
/* Output parameters description:                               */
/*   - A permutation in Sym(np^(s-1))                           */
/*                                                              */
/* Function initially developed by Roland Barrolleta and        */
/*                    generalized by Adrián Torres              */
/*                                                              */
/****************************************************************/
PermZpsToPermZp := function(permZps, p, s)	
    n := Degree(Parent(permZps));
    permList := &cat[ [p^(s-1)*(i^permZps)-(p^(s-1)-j) : j in [1..p^(s-1)]]
                                                       : i in [1..n]];
    return Sym(p^(s-1)*n)!permList;	
end function;

/****************************************************************/
/*                                                              */  
/* Function name: MatrixToPermZps                               */
/* Parameters: M, G, T, p                                       */
/* Function description: Given a matrix M, a generator matrix   */
/*   G of the Hadamard code over Z/p^s of type (n;t1,...,ts),   */
/*   and a primer number p return the permutation in Sym(n)     */
/*   associated to M.                                           */
/* Input parameters description:                                */
/*   - M : A matrix in GL(t1+...+ts)                            */
/*   - G : A generator matrix of a Hadamard code over Z/p^s     */
/*   - T : A sequence with the type t1,...,ts                   */
/*   - p : A prime number                                       */
/* Output parameters description:                               */
/*   - A permutation in Sym(n)                                  */
/*                                                              */
/* Remark: The transpose of matrix M is multiplied by the       */
/*   generator matrix G. This results in a matrix containing the*/
/*   same column vectors, but permutated. One can identify each */
/*   column with the corresponding original coordinate position */
/*   by adding its elements in a certain linear combination.    */
/*   Then, adding the rows of the permutated matrix with the    */
/*   same linear combination, we obtain a tuple describing a    */
/*   permutation which is equivalent to matrix M.               */
/*                                                              */
/* Reference: "Partial permutation decoding and PD-sets for     */
/*   Z/p^s-linear generalized Hadamard codes," by A. Torres and */
/*   M. Villanueva, Finite Fields Their Appl., 102316, 2023.    */
/*                                                              */
/* Function developed by Adrián Torres                          */
/*                                                              */
/****************************************************************/
MatrixToPermZps := function(M, G, T, p)
    n := Ncols(G);

    permGcol := Matrix(Z, Transpose(M)*G);
    permGcol := Matrix(RealField(), permGcol);
    s := #T;
    positions := permGcol[1];
    T[1] -:= 1;
    row := 1;
    power := 0;
    for k in [1..s] do
        positions +:= (T[k] gt 0) select &+[p^power*p^(s*(i-1)-(k-1)*i)*permGcol[row+i] : i in [1..T[k]]]
                        else 0*permGcol[1];
        row +:= T[k];
        power +:= (s-k+1)*T[k];
    end for;
    positions := Matrix(Z, positions);
    return Sym(n)!Eltseq(positions);   
end function;

/****************************************************************/
/*                                                              */
/* Function name: ExtendPerm                                    */
/* Parameters: perm, p                                          */
/* Function description: Given a permutation perm in Sym(n),    */
/*   returns the permutation in Sym(pn) obtained by applying p  */
/*   in each of {1..n}, {n+1..2n}, ..., {(p-1)n+1..pn}.         */
/* Input parameters description:                                */
/*   - perm : A permutation in Sym(n)                           */
/*   - p : A prime number or a power of a prime number, which   */
/*         represents to apply the function several times       */
/* Output parameters description:                               */
/*   - A permutation in Sym(pn)                                 */
/*                                                              */
/* Function developed by Adrián Torres                          */
/*                                                              */
/****************************************************************/
function ExtendPerm(perm, p)
	pseq := Eltseq(perm);
	n    := #pseq; 
	return Sym(p*n)!(pseq cat &cat[[i+n*j : i in pseq] : j in [1..p-1]] );
end function;

/****************************************************************/
/*                                                              */
/* Function name: TimesExtendPerm                               */
/* Parameters: T, indexTi                                       */
/* Function description: Given a sequence with type T=[t1,..,ts]*/
/*   of a Hadamard code H over Z/p^s, and an integer indexTi in */
/*   {1..s}, returns the number of times is necessary to apply  */
/*   ExtendPerm function to a permutation suitable for a        */
/*   Hadamard code of type (t1,0...0) or (1,0..0,ti,0..0) to    */
/*   obtain a permutation for the Hadamard code of type T.      */
/* Input parameters description:                                */
/*   - T : A sequence with the type t1,...,ts                   */
/*   - indexTi : An integer in {1..s}                           */
/* Output parameters description:                               */
/*   - An integer                                               */
/*                                                              */
/****************************************************************/
function TimesExtendPerm(T, indexTi)
    s := #T;
    T[1] -:= 1;
    return &+[T[i]*(s-i+1) : i in {1..s} diff {indexTi}];
end function;

/****************************************************************/
/*                                                              */  
/* Function name: FindPDSetRandomHadamard                       */
/* Parameters: T, allzero, V, seqGLZps, p, r                    */
/* Function description: Given a sequence T=[t1,...,ts], the    */
/*   all-zero vector of length t1+···+ts, a sequence of vectors */
/*   V for the first row, a sequence containing the general     */
/*   linear groups GL(t1-1,Z/p^s), GL(t2, Zp^(s-1)), ...,       */
/*   GL(ts,Zp), the prime numver p, and an integer r, the       */
/*   function tries to find a sequence of r+1                   */
/*   invertible matrices fulfilling certain conditions. The     */
/*   function returns two parameters: first whether the sequence*/
/*   of matrices is found or not, and then the found sequence.  */
/*   This function is called by the below function              */
/*   PDSetHadamardCodeZps_Random.                               */
/* Input parameters description:                                */
/*   - T : A sequence with the type t1,...,ts                   */
/*   - allzero : The all-zero vector of length t1+...+ts        */
/*   - V : A sequence of vectors to construct the first row     */
/*   - seqGLZps : A sequence with the general linear groups     */
/*   - p : A primer number                                      */
/*   - r : A positive integer                                   */
/* Output parameters description:                               */
/*   - A sequence of invertible matrices in GL(t1+....+ts)      */
/*                                                              */
 /* Function initially developed by Roland Barrolleta and       */
/*                    generalized by Adrián Torres              */
/*                                                              */
/****************************************************************/ 
function FindPDSetRandomHadamard(T, allzero, V, seqGLZps, p, r)    
    s := #T;
    Zps := Integers(p^s);
    sumT := &+T;

    PDSetMatrices := [];
    allStarRows := {};

    // Find a total of (r+1) matrices M from the PAut of the Hadamard code
    repeat

        // Find a matrix M from PAut, having a (star)M with no rows in allStarRows 
        iniTime := Cputime();
        repeat

            // Construct a random matrix from PAut 					
            if (sumT-T[1]) eq 0 then
                M := HorizontalJoin(allzero, VerticalJoin(Random(V), Random(seqGLZps[1])));
            else
                Tminus := T;
                Tminus[1] -:= 1;
                BlockList := [**];
                for i in [1..s] do
                    k := s-i+1;
                    Zpk := Integers(p^k);
                    BlockList := BlockList cat [*
                        HorizontalJoin(<RandomMatrix(Zpk, Tminus[i], Tminus[j]) : j in [1..i-1]> 
                                       cat <Random(seqGLZps[i])> cat
                               <p^(j-i)*RandomMatrix(Zpk, Tminus[i], Tminus[j]) : j in [i+1..s]>)
                    *];
                end for;
                Block := VerticalJoin(<Matrix(Zps, blocki) : blocki in BlockList>);
                M := HorizontalJoin(allzero, VerticalJoin(Random(V), Block));
            end if;
            M[1,1] := 1;

            // Construct a set with the rows of (star)M
            N := M^(-1);
            MStarRows := { N[1] } join { N[1] + N[i] : i in [2..T[1]] } join
                 &join{ {N[1]+p^(k-1)*N[i] : i in [&+T[1..(k-1)]+1..&+T[1..k]]} : k in [2..s]};     

            // Check whether the matrix M has the required property, that is, whether 
            // the rows of (star)M are all different and disjoint from the allStarRows set
            Mfound := (#MStarRows eq sumT) and IsDisjoint(allStarRows, MStarRows);
            if Mfound then
                allStarRows := allStarRows join MStarRows;
                for row in MStarRows do
                    Exclude(~V, V!Remove(Eltseq(row), 1));
                end for;
            end if;
        until Mfound or (Cputime(iniTime) gt MAXTIME);

        if Mfound then
            Append(~PDSetMatrices, M);
        else
            // If not found, return false and an empty sequence
            return false, [];
        end if;     

    until (#PDSetMatrices eq (r+1));
    
    return true, PDSetMatrices;
end function;

/****************************************************************/
/*                                                              */
/* Function name: PDSetHadamardCodeZps_Random                   */
/* Parameters: p, T                                             */
/* Function description: Given a sequence T=[t1,...,ts] and a   */
/*   prime number p, construct a sequence of r+1 invertible     */
/*   matrices which give an r-PD-set for the Hadamard code H    */
/*   over Z/p^s of type T given by the ZpHadamardCode function. */
/*   The function returns two sequences containing the          */
/*   permutations in Sym(n) and Sym(np^(s-1)) associated to the */
/*   sequence of invertible matrices, and the sequence of       */
/*   invertible matrices found.                                 */
/*   The function gives an r-PD-set of size r+1 such that       */
/*   ftilda_p^(t1,...,ts) <= r <= f_p^(t1,...,ts), where        */
/*   ftilda_p^(t1,...,ts)=max(f_p^(t1,0,..,0), f_p^(1,t2,0...0),*/
/*   ..., f_p^(1,0...0,ts))                                     */
/*   The function starts from the maximum value of r and        */
/*   decreases it when the r-PD-set is not found after a time   */
/*   out. In case the parameter r=ftilda_p^(t1,...,ts), the     */
/*   r-PD-set is constructed using the deterministic method.    */
/* Input parameters description:                                */
/*   - T : A sequence with the type t1,...,ts                   */
/*   - p : A primer number                                      */
/* Output parameters description:                               */
/*   - A sequence of permutations in Sym(n), which form an      */
/*     r-PD-set for H                                           */
/*   - A sequence of permutations in Sym(np^(s-1)), which form  */
/*     an r-PD-set for Hp=Phi_s(H)                              */
/*   - A sequence of invertible matrices in GL(t1+...+ts), which*/
/*     are equivalent to the permutations given in output 1     */
/*                                                              */
/* Reference: "Partial permutation decoding and PD-sets for     */
/*   Z/p^s-linear generalized Hadamard codes," by A. Torres and */
/*   M. Villanueva, Finite Fields Their Appl., 102316, 2023.    */ 
/*                                                              */
/* Function initially developed by Roland Barrolleta and        */
/*                    generalized by Adrián Torres              */
/*                                                              */
/****************************************************************/
forward PDSetHadamardCodeZps_DetermFree;
forward MaximumSizePDSetHadamardFree;
forward MaximumSizePDSetHadamard;
PDSetHadamardCodeZps_Random := function(p, T)
    s := #T; 
    Zps := Integers(p^s);
    sumT := &+T;	    	
    allzero := Matrix(Zps, sumT, 1, [0^^sumT]);

    // Construct a sequence with the general linear groups
    seqGLZps := ((T[1]-1) ne 0) select [* GL(T[1]-1, Zps) *] else [* [ZeroMatrix(Zps, 0, 0)] *];
    for i in [2..s] do
        Zpsi := Integers(p^(s-i+1));
        seqGLZps := seqGLZps cat [* (T[i] ne 0) select GL(T[i], Zpsi) 
                else [ZeroMatrix(Zpsi, 0, 0)] *];
    end for;

    // Construct all vectors with t1-1 coordinates over Z/p^s,
    // t2 coordinates in pZ/p^s, t3 coordinates in p^2Z/p^s, etc.
    _, G := ZpHadamardCode(p, T);
    V := Rows(Transpose(RemoveRow(G, 1)));

    minimumR, ti, indexTi := MaximumSizePDSetHadamardFree(p, T, 0);
    maximumR := MaximumSizePDSetHadamard(p, T); 
    maxRfound := 0;

    // from the maximum possible r, a random method is applied NUMITER times 
    // before decreasing the value of r up to the minimum possible r
    r := maximumR;
    while (r gt minimumR) and (r gt maxRfound) do
        iteration := 1;

        ///////////////////////////////////////////////////////////////////
        vprintf PDSetHadamardFlag, 1: "Trying to find an %o-PD-set of size %o randomly... \n", r, r+1;
        ///////////////////////////////////////////////////////////////////

        repeat 
            Mfound, PDSetMatrices := FindPDSetRandomHadamard(T, 
                                                        allzero, V, seqGLZps, p, r);
            if #PDSetMatrices-1 gt maxRfound then
                maxRfound := #PDSetMatrices-1;
                maxPDSetMatrices := PDSetMatrices;
            end if;
            iteration +:= 1;
        until (iteration gt NUMITER) or Mfound;
        if Mfound then
            
            ///////////////////////////////////////////////////////////////////
            vprintf PDSetHadamardFlag, 1: "A %o-PD-set has been found in iteration %o of %o!! \n", 
                    r, iteration, NUMITER;
            ///////////////////////////////////////////////////////////////////

            break;
        else
            r -:= 1;
        end if;
    end while;

    // when the minimum r is achieved, the deterministic method is applied
    if (r eq minimumR) then

        ///////////////////////////////////////////////////////////////////
        freeType := [1] cat [0^^(s-1)]; 
        freeType[indexTi] := ti;
        vprintf PDSetHadamardFlag, 1: 
                "A %o-PD-set of size %o has been obtained by applying a recursive construction\n", r, r+1;
        vprintf PDSetHadamardFlag, 1: "  from a %o-PD-set for a Hadamard code of type %o. \n", r, freeType; 
        ///////////////////////////////////////////////////////////////////

        S, Sp, PDSetMatrices := PDSetHadamardCodeZps_DetermFree(ti, p, s, indexTi);
        timesExtension := TimesExtendPerm(T, indexTi);
        S := [ExtendPerm(perm, p^timesExtension) : perm in S];
        Sp := [ExtendPerm(pperm, p^timesExtension) : pperm in Sp];
    	return S, Sp, PDSetMatrices;
	end if;

    // Construction of the corresponding permutations in Sym(n)
    PDSetPermsZps := [MatrixToPermZps(M, G, T, p) : M in maxPDSetMatrices];

    // Construction of the corresponding permutations in Sym(p^(s-1)*n)
    PDSetPermsZp := [PermZpsToPermZp(permZps, p, s) : permZps in PDSetPermsZps];

    return PDSetPermsZps, PDSetPermsZp, maxPDSetMatrices;
end function;

/****************************************************************/
/*                                                              */
/* Function name: GaloisRingSequence                            */
/* Parameters: p, s, m                                          */
/* Function description: Return a sequence R with all the       */
/*   elements in the Galois ring over Z/p^s of dimension m>1    */
/*   using the p-adic representation and a lexicographical      */
/*   order.                                                     */
/* Input parameters description:                                */
/*   - p : A prime number                                       */
/*   - s : A positive integer                                   */
/*   - m : A positive integer m greater than 1                  */
/* Output parameters description:                               */
/*   - A sequence containing all elements in R,                 */
/*     lexicographically ordered                                */
/*                                                              */
/* Function developed by Adrián Torres                          */
/*                                                              */
/****************************************************************/
function GaloisRingSequence(p, s, m)
    ZX<z> := PolynomialRing(Z);
    ZpX := PolynomialRing(Integers(p));
    Zps := Integers(p^s);
    ZpsX := PolynomialRing(Zps);
    
    polyG := z^(p^m-1) - 1;
    factorsG := [w[1] : w in Factorization(ZpX!polyG)];
    henselfact := HenselLift(ZX!polyG, factorsG, ZpsX);
    for i in [1..#henselfact] do
        polyF := henselfact[i];
        if (Degree(polyF) eq m) and (IsPrimitive(PolynomialRing(GF(p))!polyF)) then
            break;
        end if;
    end for;
    GR<a> := GaloisRing(p^s, ZX!polyF);

    T := [0] cat [a^i : i in [0..p^m-2]];
    R := [&+[r[i]*p^(s-i) : i in Reverse([1..s])] : r in CartesianProduct([T^^s])];
    return R;
end function;

/****************************************************************/
/*                                                              */
/* Function name: PDSetHadamardCodeZps_DetermFree               */
/* Parameters: ti, p, s, i                                      */
/* Function description: Given an integer ti, a prime number p, */
/*   and two integers s and i, construct a sequence of r+1      */
/*   invertible matrices which give an r-PD-set for the Hadamard*/
/*   code over Z/p^s of type (ti,0..0) if i=1 or (1,0..0,ti,0.. */
/*   0) otherwise, given by the ZpHadamardCode function.        */
/*   The function returns two sequences containing the          */
/*   permutations in Sym(n) and Sym(np^(s-1)) associated to the */
/*   sequence of invertible matrices. Finally, it also returns  */
/*   the sequence of r+1 invertible matrices, where             */
/*   r=floor((p^(s*t1-s)-t1)/t1) if i=1 or otherwise            */
/*   r=floor((p^((s-i+1)ti)-t1-1)/(ti+1)). The rows are not     */
/*   always given in their more reduced form, mod p^(s-i+1).    */
/* Input parameters description:                                */
/*   - ti : A positive integer                                  */
/*   - p : A prime number                                       */
/*   - s : A positive integer                                   */
/*   - i : A positive integer                                   */
/* Output parameters description:                               */
/*   - A sequence of permutations in Sym(n) which form an       */
/*     r-PD-set for H                                           */
/*   - A sequence of permutations in Sym(p^(s-1)*n) which form  */
/*     an r-PD-set for Hp=Phi_s(H)                              */
/*   - A sequence of invertible matrices in GL(ti) which are    */
/*     equivalent to the permutations given in output 1         */
/*                                                              */
/* Reference: "Partial permutation decoding and PD-sets for     */
/*   Z/p^s-linear generalized Hadamard codes," by A. Torres and */
/*   M. Villanueva, Finite Fields Their Appl., 102316, 2023.    */ 
/*                                                              */
/* Function initially developed by Roland Barrolleta and        */
/*                    generalized by Adrián Torres              */
/*                                                              */
/****************************************************************/
PDSetHadamardCodeZps_DetermFree := function(ti, p, s, i)
    newT1 := (i eq 1) select ti else ti + 1;
    newS := s-i+1;
    r := Floor((p^(newS*(newT1-1))-newT1)/newT1);
    m := newT1 - 1;
    Zps := Integers(p^s);
    
    // when m=0, that is, newT1=1 and i=1 
    if IsZero(m) then
        PDSetMatrices := [Matrix(Zps, [[1]])];

    // when m=1, that is, newT1=2
    elif IsOne(m) then
        PDSetMatrices := [];
        for j in [0..r] do
            M := Matrix(Zps, [[1,p^(i-1)*(2*j)], [0,1]]); 
            Append(~PDSetMatrices, M^(-1));
        end for;

    // when m>1, that is, ti>1
    else
        R := GaloisRingSequence(p, newS, m);

        zerocolumn := Matrix(Zps, newT1, 1, [0^^newT1]);
        PDSetMatricesInv := [ ];
        PDSetMatrices := [ ];
        for k in [0..r] do
            a := R[newT1*k + 1];
            row := p^(i-1)*Matrix(Zps, [Eltseq(a)]);
            Bi := Matrix(Zps, [Eltseq(x) : x in [R[newT1*k + j] - a : j in [2..newT1]]]);
            Bi := HorizontalJoin(zerocolumn, VerticalJoin(row, Bi));
            Bi[1][1] := 1;
            Append(~PDSetMatricesInv, Bi); 
            Append(~PDSetMatrices, Bi^(-1));
        end for;
    end if;
    
    // Type of the Hadamard code H over Z/p^s
    T := [0^^s];
    T[1] := 1;
    T[i] := ti;
    // Construction of the corresponding permutations in Sym(n)
	_, G := ZpHadamardCode(p, T);
    PDSetPermsZps := [MatrixToPermZps(M, G, T, p) : M in PDSetMatrices];
    // Construction of the corresponding permutations in Sym(p^(s-1)*n)
    PDSetPermsZp := [PermZpsToPermZp(permZps, p, s) : permZps in PDSetPermsZps];

    return PDSetPermsZps, PDSetPermsZp, PDSetMatrices;
end function;

/****************************************************************/
/*                                                              */
/* Function name: obtainAlpha                                   */
/* Parameters: p, T                                             */
/* Function description: Given a prime number p and a sequence  */
/*   T=[t1,...,ts], return the maximum integer that satisfies   */
/*   the conditions given in the main theorem of the reference  */
/*   to construct an r-PD-set for the Hadamard code over Z/p^s  */
/*   of type [n;t1,...,ts].                                     */
/*   It is used in function PDSetHadamardCodeZps_DetermGeneral  */
/* Input parameters description:                                */
/*   - p : A primer number                                      */
/*   - T : A sequence with the type t1,...,ts                   */
/* Output parameters description:                               */
/*   - An integer with the maximum value                        */
/*   - A sequence with the values of d                          */
/*   - A sequence with the values of numClass                   */
/*                                                              */
/* Reference: "Explicit construction of r-PD-sets for non-free  */
/*   Z/p^s-linear generalized Hadamard codes," J. Rifà, A.      */
/*   Torres-Martín and M. Villanueva, submitted to Finite Fields*/
/*   Their Appl., 2023.                                         */ 
/*                                                              */
/* Function developed by Adrián Torres                          */
/*                                                              */
/****************************************************************/
function obtainAlpha(p, T)
    s := #T;
    numClass1 := p^(T[1]-1);
    numClass := [numClass1^i : i in [1..s-1]];

    dp := Floor(numClass[1]/T[1]);
    d := [dp] cat [dp*numClass[i] : i in [1..s-2]];

    alpha := 0;
    repeat
        alpha +:= d[s-1];
        m := [(T[i+1] gt 0) 
                select d[i]*T[1]*Floor((numClass[s-i]-alpha/d[i])/(&+T[i+1..s])) 
                else d[s-1]*numClass[1]: i in [1..s-1]
             ];
    until &or[alpha gt mi : mi in m];

    return alpha-d[s-1], d, numClass;
end function;

/****************************************************************/
/*                                                              */
/* Function name: InverseMatrixPDSetHadamard                    */
/* Parameters: M, invA                                          */
/* Function description: Given a matrix M such that is of the   */
/*   form [[1,v], [0,A]] and the inverse of matrix A, return    */
/*   the inverse of matrix M.                                   */
/*   It is used in function PDSetHadamardCodeZps_DetermGeneral  */
/* Input parameters description:                                */
/*   - M : A matrix                                             */
/*   - Ainv: A matrix, which is the inverse of a submatrix of M */
/* Output parameters description:                               */
/*   - The inverse matrix of M                                  */                      
/*                                                              */
/****************************************************************/
function InverseMatrixPDSetHadamard(M, Ainv)
    Minv := M; 
    InsertBlock(~Minv, Ainv, 2, 2); 
    firstRowM := SubmatrixRange(M, 1, 2, 1, Ncols(M));
    firstRow := -firstRowM * Ainv;
    InsertBlock(~Minv, firstRow, 1, 2); 
    return Minv; 
end function;

/****************************************************************/
/*                                                              */
/* Function name: PDSetHadamardCodeZps_DetermGeneral            */
/* Parameters: p, T, newT                                       */
/* Function description: Given a prime number, a sequence       */
/*   T=[t1,...,ts] and a sequence newT=[r1,...,rs'] such that   */
/*   r1>=2, construct a sequence of r+1 invertible              */
/*   matrices which give an r-PD-set for the Hadamard code H    */
/*   over Z/p^s of type T given by the ZpHadamardCode function. */
/*   The function returns two sequences containing the          */
/*   permutations in Sym(n) and Sym(np^(s-1)) associated to the */
/*   sequence of invertible matrices, and a sequence with the   */
/*   invertible matrices.                                       */
/*   The function gives an r-PD-set of size r+1 with an r that  */
/*   satisfy the conditions given in the main theorem of the    */
/*   reference.                                                 */
/*   If t1>1, then newT=T, otherwise newT=[1+ti,ti+1,...,ts]    */
/*   with T = [1,0,..,0,ti,ti+1,..., ts].                       */               
/* Input parameters description:                                */
/*   - p : A primer number                                      */
/*   - T : A sequence with the type t1,...,ts,                  */
/*   - newT : A sequence [r1,...,rs'] with r1>1                 */
/* Output parameters description:                               */
/*   - A sequence of permutations in Sym(n), which form an      */
/*     s-PD-set for H                                           */
/*   - A sequence of permutations in Sum(p^(s-1)), which form   */
/*     an s-PD-set for Hp=Phi_s(H)                              */
/*   - A sequence of invertible matrices in GL(t1), which are   */
/*     equivalent to the permutations given in output 1         */
/*                                                              */
/* Reference: "Explicit construction of r-PD-sets for non-free  */
/*   Z/p^s-linear generalized Hadamard codes," J. Rifà, A.      */
/*   Torres-Martín and M. Villanueva, submitted to Finite Fields*/
/*   Their Appl., 2023.                                         */ 
/*                                                              */
/* Function developed by Adrián Torres                          */
/*                                                              */
/****************************************************************/
function PDSetHadamardCodeZps_DetermGeneral(p, T, newT)
    difSizeType := #T - #newT;

    alpha, d, numClass := obtainAlpha(p, newT);
    M := [Z!(alpha/di) : di in d];
    r := numClass[1]-d[1]*newT[1];
    newS := #newT;
    s := #T;

    // Construct matrices Ap, Ap^2,..., Ap^(s-1) given in the reference
    Ap := [[lambda*numClass[1]+newT[1]*mu + x : lambda in [M[1]..numClass[newS-1]-1]] : 
                                             x in [1..newT[1]], mu in [0..d[1]-1]];
    A := [Ap];
    for i in [2..newS-1] do
        Api := [[lambda*numClass[i]+kappa*numClass[1]+newT[1]*mu + x : 
                                             lambda in [M[i]..numClass[newS-i]-1]] :
                  x in [1..newT[1]], mu in [0..d[1]-1], kappa in [0..numClass[i-1]-1]];
        Append(~A, Api);
    end for;

    // Select indices for the first t1 rows and last ts rows of each matrix
    PDSetMatricesIndex := [];
    for mu in [0..d[1]-1] do
        indices := (newS gt 2) select [{0..numClass[1]-1} : i in [1..newS-2]] else [{0}];
        carProd := CartesianProduct(indices);
        for kappa in carProd do
            for lambda in [0..M[newS-1]-1] do
                // Construct a block of newT[1] consecutive elements
                index := lambda*numClass[newS-1]+ ((newS gt 2) select &+[kappa[i]*numClass[newS-i-1] : i in [1..newS-2]] else 0)+mu*newT[1];
                matrix := [index + 1..index + newT[1]]; 
                // Determine position of the leader
                position := (lambda+((newS gt 2) select &+[kappa[i]*M[i+1] : i in [1..newS-2]] else 0)+mu*M[1]) mod newT[1]+1; 
                // Move the leader to the first row
                firstRow := matrix[1];
                matrix[1] := matrix[position]; 
                matrix[position] := firstRow;

                // Determine row of Aps-1 corresponding to the leader
                rowA := (matrix[1]-1 - ((matrix[1]-1) div numClass[1])*r) mod (newT[1]*d[newS-1]) + 1; 
                for gamma in [1..newT[newS]] do
                    x := A[newS-1][rowA][1];
                    matrix cat:= [x];
                    for i in [1..newS-1] do
                        if (newT[i+1] gt 0) then
                            Exclude(~A[i][((rowA-1) mod (d[i]*newT[1]))+1], x);
                        end if;
                    end for;
                end for;
                Append(~PDSetMatricesIndex, matrix);
            end for;
        end for;
    end for;

    // Select indices for the remaining t2 rows of each matrix
    for i in [2..newS-1] do
        for m in [1..#PDSetMatricesIndex] do
            leader := PDSetMatricesIndex[m][1];
            rowAi := (leader-1 - ((leader-1) div numClass[1])*r) mod (newT[1]*d[newS-i]) + 1;
            for gamma in [1..newT[newS-i+1]] do
                x := A[newS-i][rowAi][1];
                Append(~PDSetMatricesIndex[m], x);
                Exclude(~A[newS-i][rowAi], x);
                for j in [1..newS-i] do
                    Exclude(~A[j][((rowA-1) mod (d[j]*newT[1]))+1], x);
                end for;
            end for;
        end for;
    end for;

    // when t1=2
    if newT[1] eq 2 then
        R := [0..p^newS-1];

    // when t1>2
    else
        // Create the Galois Ring of dimension t1-1 over Z/p^s
        R := GaloisRingSequence(p, newS, newT[1]-1);
    end if;

    // Construct all possible vectors for the first row that lead to different 
    // variations for each constructed matrix
    _, G := ZpHadamardCode(p, [1] cat newT[2..newS]);
    V := Rows(Transpose(RemoveRow(G, 1)));

    // Construct all matrices in the PDSet, with their corresponding variations
    PDSetMatrices := [];
    sumT := &+newT;
    sumTminusT1 := sumT-newT[1];
    newZps := Integers(p^newS);
    Zps := Integers(p^s);
    for matrix in PDSetMatricesIndex do
        listRows := matrix[1..newT[1]] cat Reverse(matrix[newT[1]..sumT]);
        firstrow := Matrix(newZps, [[1] cat Eltseq(R[listRows[1]])]);
        rowsZpi := <(newT[1] gt 1) select Matrix(newZps, [[0] cat Eltseq(x) : 
                                       x in [(R[listRows[j]]-R[listRows[1]]) : 
                                       j in [2..newT[1]]]])
                                else ZeroMatrix(newZps, newT[1]-1, newT[1])> cat
                   <(newT[i] gt 0) select Matrix(newZps, [[0] cat Eltseq(x) : 
                                       x in [(R[listRows[j]]-R[listRows[1]]) div p^(i-1) : 
                                       j in [&+newT[1..i-1]+1..&+newT[1..i]]]])
                                else ZeroMatrix(newZps, newT[i], newT[1]) : i in [2..newS]>;

        blockZpi := <ZeroMatrix(newZps, newT[1]-1, sumTminusT1)> cat
                    <HorizontalJoin(<ZeroMatrix(newZps, newT[i], &+newT[2..i-1]), 
                                     IdentityMatrix(newZps, newT[i]),
                                     ZeroMatrix(newZps, newT[i], &+newT[i+1..newS])>) : i in [2..newS]>;

        zeroVector := Vector(newZps, [0^^sumTminusT1]);
        M_0 := HorizontalJoin(VerticalJoin(<firstrow> cat rowsZpi),
                              VerticalJoin(<zeroVector> cat blockZpi));
        subM_0_inv := Submatrix(M_0, 2, 2, sumT-1, sumT-1)^(-1);

        if not IsZero(difSizeType) then 
            subM_0_inv := ChangeRing(subM_0_inv, Zps);
            
            M_0 := MultiplyRow(ChangeRing(M_0, Zps), p^difSizeType, 1);
            M_0[1,1] := 1;
        end if;

        // Construct all possible variations  
        for vector in V do 
            if IsZero(difSizeType) then 
                M := InsertBlock(M_0, vector, 1, newT[1]+1);
            else
                M := InsertBlock(M_0, p^difSizeType*ChangeRing(vector, Zps), 1, newT[1]+1);
            end if; 
            // Compute the inverse of M using the inverse of a submatrix 
            Append(~PDSetMatrices, InverseMatrixPDSetHadamard(M, subM_0_inv));
        end for;
    end for;

    // Construction of the corresponding permutations in Sym(n)
    _, G := ZpHadamardCode(p, T);
    PDSetPermsZps := [MatrixToPermZps(M, G, T, p) : M in PDSetMatrices];
    // Construction of the corresponding permutations in Sym(p^(s-1)*n)
    PDSetPermsZp := [PermZpsToPermZp(permZps, p, s) : permZps in PDSetPermsZps];

    return PDSetPermsZps, PDSetPermsZp, PDSetMatrices;
end function;

/****************************************************************/
/*                                                              */
/* Function name: InformationSetHadamard                        */
/* Parameters: p, T                                             */
/* Function description: Given a prime number p and a sequence  */
/*   with the type [t1,..,ts] of a Hadamard code H over Z/p^s,  */
/*   return an information set I for the code H, generated by G,*/
/*   where G is the matrix given by the ZpHadamardCode          */
/*   function. It also returns an information set Ip for the    */
/*   corresponding code Hp=Phi_s(H) over Zp.                    */
/* Input parameters description:                                */
/*   - p : A prime number                                       */
/*   - T : A sequence with the type t1,...,ts                   */
/* Output parameters description:                               */
/*   - A sequence of coordinate positions in [1..Length(C)]     */
/*     containing an information set for H                      */
/*   - A sequence of coordinate positions in [1..               */
/*     p^(s-1)*Length(C)] containing an information set for Hp  */
/*                                                              */
/* Reference: "Partial permutation decoding and PD-sets for     */
/*   Z/p^s-linear generalized Hadamard codes," by A. Torres and */
/*   M. Villanueva, Finite Fields Their Appl., 102316, 2023.    */ 
/*                                                              */
/* Function developed by Adrián Torres                          */
/*                                                              */
/****************************************************************/
function InformationSetHadamard(p, T)
    s := #T;
    T[1] -:= 1;
    n := 1;
    I := [1];
    for i in [1..s] do
        powerP := p^(s-i+1);
        for j in [1..T[i]] do
            Append(~I, n+1);
            n *:= powerP;
        end for;
    end for;

    T[1] +:= 1;
    Ibin := &cat[&cat[[p^(s-1)*(i-1)+1] cat 
                       [p^(s-1)*(i-1)+1+p^(k-1+j) : j in [0..(s-k-1)]]
                       : i in I[(&+T[1..(k-1)]+1)..(&+T[1..k])]]
                 : k in [1..s]]; 
 
    return I, Ibin;
end function;

/****************************************************************/
/*                                                              */
/* Function name: ReduceTypeT1equal1                            */
/* Parameters: T                                                */
/* Function description: Given a sequence T=[1,t2,...,ts],      */
/*   return the sequence [1+ti, ti+1,...,ts], where ti is the   */
/*   first element in T[2..s] different to zero. It also returns*/
/*   the difference between the size of T and the new sequence. */ 
/*   It is used in intrinsic PDSetZpHadamardCode.               */
/* Input parameters description:                                */
/*   - T : A sequence with the type 1,t2,...,ts                 */
/* Output parameters description:                               */
/*   - A sequence with the new type 1+ti,ti+1,...,ts            */
/*   - An integer with the difference between the size of T and */
/*     the size of the new type, #T-#newT                       */
/*                                                              */
/****************************************************************/
function ReduceTypeT1equal1(T)
    s := #T;
    i := 2;
    while IsZero(T[i]) do
        i +:= 1;
    end while; 
    return [1+T[i]] cat [T[j] : j in [i+1..s]], i-1;
end function;

/****************************************************************/
/*                                                              */
/* Function name: MaximumSizePDSetHadamardFree                  */
/* Parameters: p, newT, difSizeType                             */
/* Function description: Given a prime number p, a sequence     */
/*   newT=[t1,t2,...,ts], where t1>=2, and the difference       */
/*   between the size of the original type T and the size of the*/
/*   reduced sequence, return three integers: the maximum value */
/*   r of a r-PD-set for a Hadamard code over Z/p^s of type     */
/*   (t1,0..0) or (1,0..0,ti,0..0), the value of ti; and the    */
/*   index i where this maximum is achieved.                    */
/*   From an r-PD-set of size r+1 for a Hadamard code of type   */
/*   (t1,0..0) or (1,0..0,ti,0..0), we can obtain one for a     */
/*   Hadamard code of type T, recursively.                      */  
/*   It is used in intrinsic PDSetZpHadamardCode.               */
/* Input parameters description:                                */
/*   - p : A prime number                                       */
/*   - newT: A sequence with the reduced type [t1,...,ts], that */
/*           is, a type having t1>=2                            */
/*   - difSizeType: An integer with the difference between      */
/*                 the size of T and the size of newT, #T-#newT */
/* Output parameters description:                               */
/*   - An integer r with the maximum value of the r-PD-set      */
/*   - An integer ti from [t1,...,ts]                           */
/*   - An integer i, such that the type is [t1,0..0] if i=1     */
/*                   and is [1,0..0,ti,0...0] otherwise         */                      
/*                                                              */
/****************************************************************/
function MaximumSizePDSetHadamardFree(p, newT, difSizeType)
    newS := #newT;
    newT[1] -:= 1;
    maxSizePDset, maxIndex := Max([Floor((p^((newS-i+1)*newT[i])-newT[i]-1)/(newT[i]+1)) 
                                                                      : i in [1..newS]]);
    if IsZero(difSizeType) then
        newT[1] +:= 1;
    end if;
    return maxSizePDset, newT[maxIndex], maxIndex + difSizeType;
end function;

/****************************************************************/
/*                                                              */
/* Function name: MaximumSizePDSetHadamard                      */
/* Parameters: p, T                                             */
/* Function description: Given a prime number p, and a sequence */
/*   T=[t1,t2,...,ts] with the type of a Hadamard code over     */
/*   Z/p^s, return the theoretical upper bound for the value of */
/*   r such that an r-PD-set of size r+1 may exist.             */
/*   It is used in intrinsic PDSetZpHadamardCode.               */
/* Input parameters description:                                */
/*   - p : A prime number                                       */
/*   - T : A sequence with the type t1,...,ts                   */
/* Output parameters description:                               */
/*   - An integer r with the maximum value of the r-PD-set      */                      
/*                                                              */
/****************************************************************/
function MaximumSizePDSetHadamard(p, T)
    s := #T;
    sumT := &+T;
    return Floor((p^(&+[(s-i+1)*T[i] : i in [1..s]]-s) - sumT)/sumT); 
end function;

/*******************************************************************************/
/*                                                                             */
/* Function name: PDSetZpHadamardCode                                          */
/* Parameters: p, T                                                            */
/* Function description: Given a prime number p and a sequence of nonnegative  */
/*   integers [t1,...,ts] with t1>0 and s>1, the generalized Hadamard code C   */
/*   over Z/p^s of type (n; t1,...,ts), given by function ZpHadamardCode(p,    */
/*   [t1,...,ts]), is considered. The function returns an information set I in */
/*   [1,...,n] for C together with a subset S of the permutation automorphism  */
/*   group of C such that Phi(S) is an r-PD-set for C_p=Phi_s(C) with respect  */
/*   to Phi(I), where Phi_s is the Carlet's generalized Gray map as defined    */
/*   above, Phi(I) and Phi(S) as defined in the manual. The function also      */
/*   returns the information set Phi(I) and the r-PD-set Phi(S).               */
/* Input parameters description:                                               */
/*   - p : A prime number                                                      */
/*   - T : A sequence with the type t1,...,ts                                  */
/* Output parameters description:                                              */
/*   - An information set from {1..n} for C                                    */ 
/*   - A sequence S such that Phi(S) is a r-PD-set for the Hadamard code C     */
/*   - An information set from {1..np^(s-1)} for Phi_s(C)                      */
/*   - An r-PD-set Phi(S) for the Hadamard code Phi_s(C) given as a sequence   */ 
/*   - A sequence of invertible matrices, associated to the r-PD-set           */  
/*                                                                             */
/* References: [1] "Partial permutation decoding and PD-sets for Z/p^s-linear  */
/*   generalized Hadamard codes," by A. Torres and M. Villanueva, Finite       */
/*   Fields Their Appl., 102316, 2023.                                         */
/*       [2] "Explicit construction of r-PD-sets for non-free Z/p^s-linear     */
/*   generalized Hadamard codes," J. Rifà, A. Torres-Martín and M. Villanueva, */
/*   submitted to Finite Fields Their Appl., 2023.                             */ 
/*                                                                             */
/* Function developed by Adrián Torres                                         */
/*                                                                             */
/* Signature: (<RngIntElt> p, <[ModTupRngElt]> T)                              */   
/*              -> SeqEnum[RngIntElt], SeqEnum[GrpPermElt],                    */
/*                 SeqEnum[RngIntElt], SeqEnum[GrpPermElt], SeqEnum[AlgMatElt] */
/*                                                                             */
/*******************************************************************************/
intrinsic PDSetZpHadamardCode(p::RngIntElt, T::[RngIntElt] : AlgMethod := 
             "Deterministic") -> SeqEnum[RngIntElt], SeqEnum[GrpPermElt],
             SeqEnum[RngIntElt], SeqEnum[GrpPermElt], SeqEnum[AlgMatElt]
{
Given a prime number p and a sequence of nonnegative integers [t1,...,ts] with t1>0 
and s>1, the generalized Hadamard code C over Z/p^s of type (n; t1,...,ts), given 
by function ZpHadamardCode(p, [t1,...,ts]), is considered. The function returns 
an information set I in [1,...,n] for C together with a subset S of the permutation 
automorphism group of C such that Phi(S) is an r-PD-set for C_p=Phi_s(C) with respect
to Phi(I), where Phi_s is the Carlet's generalized Gray map as defined above, Phi(I) 
and Phi(S) as defined in the manual. The function also returns the information set 
Phi(I) and the r-PD-set Phi(S).

Note that for p=2 and [t1,...,ts] = [1,0,...,0,ts] or [t1,...,ts] = [1,0,...,0,1,ts], 
we have that C_p is linear, so it is possible to find an r-PD-set of size r+1 for C_p, 
for any r ≤ floor(2^m/(m+1)), by using function PDSetHadamardCode(m) with m = ts + s 
if [t1,...,ts] = [1,0,...,0,ts] and m = ts + s + 2 if [t1,..., ts] = [1,0,...,0,1,ts].

The information sets I and Phi(I) are returned as sequences of t1+···+ts and 
s(t1-1)+(s-1)t2+···+ts integers, giving the coordinate positions that correspond 
to the information sets for C and C_p, respectively. The sets S and Phi(S) are also 
returned as sequences of elements in the symmetric groups Sym(n) and Sym(p^(s-1)*n) 
of permutations on the set [1,...,n] and [1,...,p^(s-1)*n], respectively. 

A deterministic algorithm is used by default. In this case, when t1>=2, the function
first computes r as the maximum of g_p^(t1,...,ts) and ftilda_p^(t1,...,ts), where 
g_p^(t1,...,ts) is given by the general construction descrided in reference [1] and
ftilda_p^(t1,...,ts)=max(f_p^(t1,0,..,0), f_p^(1,t2,0...0), ..., f_p^(1,0...0,ts))
<= f_p^(t1,...,ts), where f_p^(t1,...,ts) is the is the theoretical upper bound 
for the value of r such that there exists an r-PD-set of size r+1 for a Hadamard 
code over Z/p^s of type (n; t1,...,ts), given in reference [2]. Let i be the first 
index i ≥ 1 such that the maximum ftilda_p^(t1,...,ts) is achieved. Then, if 
ftilda_p^(t1,...,ts) > g_p^(t1,...,ts), it constructs an r-PD-set of size r+1 for 
the generalized Hadamard code over Z/p^s of type (n; t1,0,...,0) if i=1 or 
(n; 1,0,...,0,ti,0,...,0) otherwise, which is transformed into an r-PD-set for C, 
by using the recursive construction defined in the reference. The value of r remains 
unchanged after the recursive construction, so r = ftilda_p^(t1,...,ts). On the 
other hand, if g_p^(t1,...,ts) > ftilda_p^(t1,...,ts), the general construction 
is applied and an r-PD-set of size r+1 with r = g_p^(t1,...,ts) is obtained. When
t1=1, and there exists an index i>=2 such that ti<>0 and t2=...=t_(i-1)=0, it first 
constructs an r-PD-set of size r+1 for the generalized Hadamard code over Z/p^s 
of type (n; 1+ti, t_(i+1),...,ts), which satisfies that 1+ti >= 2, by following the 
same process as described above when t1 ≥ 2. Then, the obtained r-PD-set is 
transformed into an r-PD-set for C as described in references [1,2]. 

If the parameter AlgMethod is assigned the value "Nondeterministic", the function 
tries to improve the previous results by finding an r-PD-set of size r+1 such that
ftilda_p^(t1,...,ts) <= r <= f_p^(t1,...,ts). In this case, the function starts
from the maximum value of r = f_p^(t1,...,ts) and attempts to find an r-PD-set
within a time limit of 5.0 seconds of “user time”. This is performed 5 times, 
each time starting from an empty set and trying to construct the r-PD set 
incrementally, by adding elements randomly. If such an r-PD-set is not found, 
the value of r decreases by one and the same process takes place with the new 
value of r. The value of r keeps decreasing until an r-PD-set is found or 
r = ftilda_p^(t1,...,ts).

The verbose flag PDsetHadamardFlag is set to level 0 by default. If it is set to 
level 1, some information about the process of constructing the r-PD-set of size 
r+1 is shown. Moreover, the value of the theoretical upper bound f_p^(t1,...,ts) 
given in reference [1] is also shown.

If p = 2 and s = 2, then function PDSetHadamardCodeZ4(t1, 2t1 + t2 − 1) can also 
be used. Both function only coincide when t2 = 0. When t2 > 0, the output parameters 
I and Phi(I) coincide as sets and PDSetZpHadamardCode(2, [t1, t2]) may give a 
larger r-PD-set
}
    require IsPrime(p) : "Argument 1 must be a prime number";
    s := #T;
    require s ge 2 : "The number of elements in the sequence must be greater than or equal to 2";
    require T[1] ge 1 : "The first element in the sequence must be greater than or equal to 1";
    require Min(T) ge 0: "Argument 2 must be a sequence of nonnegative integers";
    require Type(AlgMethod) eq MonStgElt: 
                 "The optional parameter AlgMethod must be a string";

    I, Ip := InformationSetHadamard(p, T);

    ///////////////////////////////////////////////////////////////////
    maximumR := MaximumSizePDSetHadamard(p, T);
    vprintf PDSetHadamardFlag, 1: "The upper bound for an r-PD-set of size r+1 is %o. \n", 
            maximumR;
    ///////////////////////////////////////////////////////////////////

    if AlgMethod eq "Nondeterministic" then
        S, Sp, SMat := PDSetHadamardCodeZps_Random(p, T);
        return I, S, Ip, Sp, SMat;

    // if AlgMethod "Deterministic" (by default)
    else 
        max, maxI := Max(T[2..s]);
        T2 := T[2..s];
        T2[maxI] := 0;

        // when the type is (n; t1,0..0)
        if (&+T[2..s] eq 0) then

            /////////////////////////////////////////////////////////////////// 
            vprintf PDSetHadamardFlag, 1: "A %o-PD-set of size %o has been obtained.\n", 
                    maximumR, maximumR+1;
            ///////////////////////////////////////////////////////////////////
        
            S, Sp, SMat := PDSetHadamardCodeZps_DetermFree(T[1], p, s, 1);
            return I, S, Ip, Sp, SMat;

        // when the type is (n; 1,0..0,ti,0..0)
        elif (T[1] eq 1) and IsZero(T2) then

            /////////////////////////////////////////////////////////////////// 
            vprintf PDSetHadamardFlag, 1: "A %o-PD-set of size %o has been obtained.\n", 
                    maximumR, maximumR+1;
            ///////////////////////////////////////////////////////////////////

            S, Sp, SMat := PDSetHadamardCodeZps_DetermFree(max, p, s, maxI+1);
            return I, S, Ip, Sp, SMat;

        // for any other type use the same as when type is (n;t1,0..0) extending 
        // the permutations in the PD-set or use the general construction 
        else 
            newT := (T[1] eq 1) select ReduceTypeT1equal1(T) else T;
            newS := #newT;
            difSizeType := s - newS;

            // r estimation by using the general construction for non-free Z/p^s-linear 
            // Hadamard codes with t1>1
            r_determGeneral := obtainAlpha(p, newT) * p^(&+[newT[i]*(newS-i+1) : i in [2..newS]])-1;
            // r estimation by using a recursive construction for free Z/p^s-linear 
            // Hadamard codes of type [t1,0...0] or of type [1,0..0,ti,0..0]
            // the parameters ti and indexTi are with respect to the initial type T
            r_determRecursive, ti, indexTi := MaximumSizePDSetHadamardFree(p, newT, difSizeType);

            // when it is better to use a recursive construction
            if (r_determRecursive gt r_determGeneral) then  
                
                ///////////////////////////////////////////////////////////////////
                freeType := [1] cat [0^^(s-1)]; 
                freeType[indexTi] := ti;  
                vprintf PDSetHadamardFlag, 1: 
                        "A %o-PD-set of size %o has been obtained by applying a recursive construction\n", 
                        r_determRecursive, r_determRecursive+1;
                vprintf PDSetHadamardFlag, 1: "  from a %o-PD-set for a Hadamard code of type %o. \n", 
                        r_determRecursive, freeType; 
                ///////////////////////////////////////////////////////////////////

                S, Sp, SMat := PDSetHadamardCodeZps_DetermFree(ti, p, s, indexTi); 

                timesExtension := TimesExtendPerm(T, indexTi); 
                S := [ExtendPerm(perm, p^timesExtension) : perm in S];
                Sp := [ExtendPerm(pperm, p^timesExtension) : pperm in Sp];
                MatId := IdentityMatrix(BaseRing(SMat[1]), Nrows(SMat[1]) + timesExtension);
                SMat := [ InsertBlock(MatId, M, 1, 1) : M in SMat ];

                return I, S, Ip, Sp, SMat;

            // when it is better to use the general construction 
            else

                ///////////////////////////////////////////////////////////////////
                vprintf PDSetHadamardFlag, 1: 
                        "A %o-PD-set of size %o has been obtained by applying the general construction\n",
                        r_determGeneral, r_determGeneral+1;
                vprintf PDSetHadamardFlag, 1:
                        "  for a Hadamard code of type %o and then adapted for a Hadamard code of type %o. \n", 
                         newT, T;
                ///////////////////////////////////////////////////////////////////

                S, Sp, SMat := PDSetHadamardCodeZps_DetermGeneral(p, T, newT);
                return I, S, Ip, Sp, SMat;
            end if;
        end if;    
    end if;
    
end intrinsic;

/****************************************************************/
/*                                                              */  
/* Function name: PermutationDecodeAlg                          */
/* Parameters: C, Sp, r, u, Ip, grayMapC                        */
/* Function description: Given a vector u over GF(p), perform   */
/*   the permutation decoding algorithm for the code Cp=Phi_s(C)*/
/*   with information set Ip and r-PD-set Sp. If the algorithm  */
/*   succeeds and the decoded codeword is cp, then the  function*/
/*   returns true, (Phi_s)^(-1)(cp) and cp. If the algorithm    */
/*   does not suceed, then the function returns false, 0, 0.    */
/* Input parameters description:                                */
/*   - C : A linear code over Z/p^s                             */
/*   - Sp : A set of permutations in PAut(Phi_s(C))             */
/*   - r : A positive integer                                   */
/*   - u : A vector over GF(p)                                  */
/*   - Ip : A sequence of positive integers                     */
/*   - grayMapC : a map from C to Phi_s(C)                      */
/* Output parameters description:                               */
/*   - A boolean value indicating whether the algorithm has been*/
/*     successful                                               */ 
/*   - A codeword in C                                          */
/*   - A codeword in Phi_s(C)                                   */
/*                                                              */
/* Reference: "Partial permutation decoding and PD-sets for     */
/*   Z/p^s-linear generalized Hadamard codes," by A. Torres and */
/*   M. Villanueva, Finite Fields Their Appl., 102316, 2023.    */ 
/*                                                              */
/* Function initially developed by Roland Barrolleta and        */
/*                    generalized by Adrián Torres              */
/*                                                              */
/****************************************************************/
PermutationDecodeAlg := function(C, Sp, r, u, Ip, grayMapC)
	x, xp := SystematicEncoding(C, Eltseq([u[i] : i in Ip]));
	if Weight(u - xp[1]) le r then
		return true, x[1], xp[1];

	else 
		for sigma in Sp do
            u_sigma := u^sigma;
			x, xp := SystematicEncoding(C, Eltseq([(u_sigma)[i]: i in Ip]));
            if Weight(u_sigma - xp[1]) le r then
				x_invSigma := (xp[1]^(sigma^-1));
				return true, x_invSigma@@grayMapC, x_invSigma;
			end if;
		end for;
	    return false, C!0, grayMapC(C!0);

	end if;
end function;

/*******************************************************************************/
/*                                                                             */
/* Function name: ZpPermutationDecode                                          */
/* Parameters: C, Ip, S, r, u                                                  */
/* Function description: Given a linear code C over Z/p^s of type (n; t1,...,  */
/*   ts), an information set Ip in {1..np^(s-1)} for Cp=Phi_s(C) as a sequence */
/*   of coordinate positions, a sequence S such that either S or Phi(S) is an  */
/*   r-PD-set for Cp with respect to the information set Ip, an integer t in   */
/*   {1,...,t}, where t is the error-correcting capability of Cp, and a vector */
/*   u which can be defined from the ambient space U=(Z/p^s)^n or              */
/*   Up=Zp^(np^(s-1)), the function attempts to decode u in Up (or Phi_s(u) in */
/*   Up if u in U) with respect to the code Cp (or C if u in U), assuming a    */
/*   systematic encoding with respect to the information set Ip. If the        */
/*   decoding algorithm succeeds in computing a codeword u' in Cp as the       */
/*   decoded version of u in Up (or Phi_s(u) in Up if u in U), then the        */
/*   function returns true, the preimage of u' by Carlet's generalized Gray    */
/*   map Phi_s and finally u'. If the decoding algorithm does not succeed in   */
/*   decoding u, then the function returns false, the zero codeword in C and   */
/*   the zero codeword in Cp.                                                  */
/* Input parameters description:                                               */
/*   - C : A linear code over Z/p^s                                            */
/*   - I : A subset of integers in {1..n}                                      */
/*   - S : A set of permutations acting on {1..n} or {1..np^(s-1)}             */
/*   - r : An integer in {1..t}, where t is the error-correcting capability    */
/*   - u : A received vector to be decoded                                     */
/* Output parameters description:                                              */
/*   - Boolean, true if u is decoded and false otherwise                       */
/*   - Decoded codeword in C or the zero codeword                              */
/*   - Decoded codeword in Cp or the zero codeword                             */
/*                                                                             */
/* Function initially developed by Roland Barrolleta and                       */
/*                    generalized by Adrián Torres                             */
/*                                                                             */
/* case u over Z/p^s                                                           */
/* Signature: (<CodeLinRng> C, <[RngIntElt]> I,                                */
/*           <[GrpPermElt]> S, <RngIntElt> s, <ModTupFldElt> u)                */
/*           -> BoolElt, ModTupRngElt, ModTupFldElt                            */
/* case u over Zp                                                              */
/* Signature: (<CodeLinRng> C, <[RngIntElt]> I,                                */
/*           <[GrpPermElt]> S, <RngIntElt> s, <ModTupRngElt> u)                */
/*           -> BoolElt, ModTupRngElt, ModTupFldElt                            */
/*                                                                             */
/*******************************************************************************/
intrinsic ZpPermutationDecode(C::CodeLinRng, Ip::[RngIntElt], S::[GrpPermElt], 
          r::RngIntElt, u::ModTupFldElt) -> BoolElt, ModTupRngElt, ModTupFldElt
{
Given a linear code C over Z/p^s of type (n; t1,...,ts), an information set Ip in 
[1..np^(s-1)] for Cp=Phi_s(C) as a sequence of coordinate positions, a sequence S 
such that either S or Phi(S) is an r-PD-set for Cp with respect to the information 
set Ip, an integer t in [1,...,t], where t is the error-correcting capability of Cp, 
and a vector u from the ambient space Up=Zp^(np^(s-1)), the function attempts to 
decode u in Up with respect to the code Cp, assuming a systematic encoding with 
respect to the information set Ip. If the decoding algorithm succeeds in computing 
a codeword u' in Cp as the decoded version of u in Up, then the function returns true, 
the preimage of u' by Carlet's generalized Gray map Phi_s and finally u'. If the 
decoding algorithm does not succeed in decoding u, then the function returns false, 
the zero codeword in C and the zero codeword in Cp. 

The permutation decoding algorithm consists in moving all errors in a received vector
u = c + e, where u in Up, c in Cp and e in Up is an error vector with at most t 
errors, out of the information positions, that is, moving the nonzero coordinates 
of e out of the information set Ip for Cp, by using a permutation in S subset PAut(Cp) 
(or Phi(S) subset PAut(Cp) if S subset PAut(C)). If S subset PAut(C) subset Sym(n), 
then Phi(S) subset PAut(Cp) subset Sym(np^(s−1)) is computed by using the map Phi 
defined in the manual. The function does not check whether Ip is an information set
for Cp, whether S or Phi(S) is an r-PD-set for Cp with respect to Ip, or whether r ≤ t.

If C is over Z4, function ZpPermutationDecode(C, Ip, S, r, u) coincides with function 
PermutationDecode(C, I, S, r, u), which works only for linear codes over Z4. However, 
the former only accepts an information set Ip subset [1,...,2n] for C_2, while the 
latter also accepts an information set I subset [1,..,n] for C as long as the sequence 
S subset PAut(C).
}
    isOverZps, p, s := IsLinearCodeOverZps(C);
    require isOverZps and (s ge 2): "The code C must be over Z/p^s with s>1";
    require PseudoDimension(C) gt 0: "The code C cannot be the zero code";
    requirege r, 1;
    
    k := #Ip;
    n := Length(C);
    np := p^(s-1)*n; 
    require (k ge 1): "Argument 2 cannot be an empty sequence";
    require #C eq p^k: "Argument 2 is not an information set for Argument 1";
    require (Min(Ip) ge 1) and (Max(Ip) le np): 
                                    "Argument 2 should be a subset of", [1..np]; 
                                                                         
    require not IsEmpty(S): "Argument 3 cannot be an empty sequence";   
    numSymG := Degree(Parent(S[1]));  
    require ((numSymG eq n) or (numSymG eq np)): 
          "Argument 3 should contain permutations acting on a set of cardinality", 
                                                                    n, "or", np;

    require (BaseRing(u) eq GF(p)): "Argument 4 must be a vector over GF(%o)", p;
    require (OverDimension(u) eq np): "Argument 4 must be of length", np;	
    
    if numSymG eq n then
        S := [ PermZpsToPermZp(sigma, p, s) : sigma in S ];
    end if;
 
    grayMapC := CarletGrayMap(C);
    Exclude(~S, Sym(np)!1);

    return PermutationDecodeAlg(C, S, r, u, Ip, grayMapC);

end intrinsic;

/****************************************************************/
intrinsic ZpPermutationDecode(C::CodeLinRng, Ip::[RngIntElt], S::[GrpPermElt], 
          r::RngIntElt, u::ModTupRngElt) -> BoolElt, ModTupRngElt, ModTupFldElt
{
Given a linear code C over Z/p^s of type (n; t1,...,ts), an information set Ip in 
[1..np^(s-1)] for Cp=Phi_s(C) as a sequence of coordinate positions, a sequence S 
such that either S or Phi(S) is an r-PD-set for Cp with respect to the information 
set Ip, an integer t in [1,...,t], where t is the error-correcting capability of Cp, 
and a vector u from the ambient space U=(Z/p^s)^n, the function attempts to decode 
Phi_s(u) with respect to the code Cp, assuming a systematic encoding with respect 
to the information set Ip. If the decoding algorithm succeeds in computing a codeword 
u' in Cp as the decoded version of Phi_s(u) in Up=Zp^(np^(s-1)), then the function 
returns true, the preimage of u' by Carlet's generalized Gray map Phi_s and 
finally u'. If the decoding algorithm does not succeed in decoding u, then the 
function returns false, the zero codeword in C and the zero codeword in Cp. 

The permutation decoding algorithm consists in moving all errors in a received vector
u = c + e, where u in Up, c in Cp and e in Up is an error vector with at most t 
errors, out of the information positions, that is, moving the nonzero coordinates 
of e out of the information set Ip for Cp, by using a permutation in S subset PAut(Cp) 
(or Phi(S) subset PAut(Cp) if S subset PAut(C)). If S subset PAut(C) subset Sym(n), 
then Phi(S) subset PAut(Cp) subset Sym(np^(s−1)) is computed by using the map Phi 
defined in the manual. The function does not check whether Ip is an information set
for Cp, whether S or Phi(S) is an r-PD-set for Cp with respect to Ip, or whether r ≤ t.

If C is over Z4, function ZpPermutationDecode(C, Ip, S, r, u) coincides with function 
PermutationDecode(C, I, S, r, u), which works only for linear codes over Z4. However, 
the former only accepts an information set Ip subset [1,...,2n] for C_2, while the 
latter also accepts an information set I subset [1,..,n] for C as long as the sequence 
S subset PAut(C).
}
    isOverZps, p, s := IsLinearCodeOverZps(C);
    require isOverZps and (s ge 2): "The code C must be over Z/p^s with s>1";
    require PseudoDimension(C) gt 0: "The code C cannot be the zero code";
    requirege r, 1;
    
    k := #Ip;
    n := Length(C);
    np := p^(s-1)*n; 
    require (k ge 1): "Argument 2 cannot be an empty sequence";
    require #C eq p^k: "Argument 2 is not an information set for Argument 1";
    require (Min(Ip) ge 1) and (Max(Ip) le np): 
                                    "Argument 2 should be a subset of", [1..np]; 
                                                                         
    require not IsEmpty(S): "Argument 3 cannot be an empty sequence";   
    numSymG := Degree(Parent(S[1]));  
    require ((numSymG eq n) or (numSymG eq np)): 
          "Argument 3 should contain permutations acting on a set of cardinality", 
                                                                    n, "or", np;
    Zps := Integers(p^s);                                                                
    require (BaseRing(u) eq Zps): "Argument 4 must be a vector over Z/%o", p^s;
    require (OverDimension(u) eq n): "Argument 4 must be of length", n;	
    
    if numSymG eq n then
        S := [ PermZpsToPermZp(sigma, p, s) : sigma in S ];
    end if;
 
    grayMapC := CarletGrayMap(C);
    Exclude(~S, Sym(np)!1);
    up := CarletGrayMap(UniverseCode(Zps, n))(u);

    return PermutationDecodeAlg(C, S, r, up, Ip, grayMapC);

end intrinsic;

/*******************************************************************************/
/*                                                                             */
/* Function name: ZpPermutationDecode                                          */
/* Parameters: C, Ip, S, r, Q                                                  */
/* Function description: Given a linear code C over Z/p^s of type (n; t1,...,  */
/*   ts), an information set Ip in {1..np^(s-1)} for Cp=Phi_s(C) as a sequence */
/*   of coordinate positions, a sequence S such that either S or Phi(S) is an  */
/*   r-PD-set for Cp with respect to the information set Ip, an integer t in   */
/*   {1,...,t}, where t is the error-correcting capability of Cp, and a        */
/*   sequence Q of vectors which can be defined from the ambient space         */
/*   U=(Z/p^s)^n or Up=Zp^(np^(s-1)), the function attempts to decode all      */
/*   vectors u in Q if Q subset Up (or Phi_s(u) in Q if Q subset U), with      */
/*   respect to the code Cp (or C if Q subset U), assuming a systematic        */
/*   encoding with respect to the information set Ip. The function returns     */
/*   three values: a sequence of booleans representing whether the decoding    */
/*   process have been successful for each u in Q, a sequence of codewords of  */
/*   C and a sequence of codewords of Cp. For each u in Q, if the decoding     */
/*   algorithm succeeds in computing a codeword u′ in Cp as the decoded version*/
/*   of u in Up (or Phi_s(u) in Up if u in U), then it returns the preimage of */
/*   u′ by Carlet’s generalized Gray map Phi_s in the second sequence, and u′  */
/*   in the third sequence. If the decoding algorithm does not succeed in      */
/*   decoding u, then the function returns the zero codeword in C and the zero */
/*   codeword in Cp in the corresponding positions of the second and third     */
/*   sequences, respectively.                                                  */
/* Input parameters description:                                               */
/*   - C : A linear code over Z/p^s                                            */
/*   - I : A subset of integers in {1..n}                                      */
/*   - S : A set of permutations acting on {1..n} or {1..np^(s-1)}             */
/*   - r : An integer in {1..t}, where t is the error-correcting capability    */
/*   - Q : S sequence of received vectors to be decoded                        */
/* Output parameters description:                                              */
/*   - Boolean, true if u is decoded and false otherwise                       */
/*   - Sequence of decoded codewords in C                                      */
/*   - Sequence of decoded codewords in Cp                                     */
/*                                                                             */
/* Function initially developed by Roland Barrolleta and                       */
/*                    generalized by Adrián Torres                             */
/*                                                                             */
/* case sequence Q over Z/p^s                                                  */
/* Signature: (<CodeLinRng> C, <[RngIntElt]> I,                                */
/*          <[GrpPermElt]> S, <RngIntElt> s, <[ModTupFldElt]> Q)               */
/*           -> [BoolElt], [ModTupRngElt], [ModTupFldElt]                      */
/* case sequence Q over Zp                                                     */
/* Signature: (<CodeLinRng> C, <[RngIntElt]> I,                                */
/*          <[GrpPermElt]> S, <RngIntElt> s, <[ModTupRngElt]> Q)               */
/*           -> [BoolElt], [ModTupRngElt], [ModTupFldElt]                      */
/*                                                                             */
/*******************************************************************************/
intrinsic ZpPermutationDecode(C::CodeLinRng, Ip::[RngIntElt], S::[GrpPermElt], 
          r::RngIntElt, Q::[ModTupFldElt]) -> SeqEnum[BoolElt], SeqEnum[ModTupRngElt],
          SeqEnum[ModTupFldElt]
{
Given a linear code C over Z/p^s of type (n; t1,...,ts), an information set Ip in 
[1..np^(s-1)] for Cp=Phi_s(C) as a sequence of coordinate positions, a sequence S 
such that either S or Phi(S) is an r-PD-set for Cp with respect to the information 
set Ip, an integer t in [1,...,t], where t is the error-correcting capability of Cp, 
and a sequence Q of vectors from the ambient space Up=Zp^(np^(s-1)), the function 
attempts to decode all vectors u in Q subset Up, with respect to the code Cp, assuming 
a systematic encoding with respect to the information set Ip. The function returns     
three values: a sequence of booleans representing whether the decoding process 
have been successful for each u in Q, a sequence of codewords of C and a sequence 
of codewords of Cp. For each u in Q, if the decoding algorithm succeeds in computing 
a codeword u′ in Cp as the decoded version of u in Up, then it returns the preimage
of u′ by Carlet’s generalized Gray map Phi_s in the second sequence, and u′ in the 
third sequence. If the decoding algorithm does not succeed in decoding u, then the 
function returns the zero codeword in C and the zero codeword in Cp in the 
corresponding positions of the second and third sequences, respectively.

The permutation decoding algorithm consists in moving all errors in a received vector
u = c + e, where u in Up, c in Cp and e in Up is an error vector with at most t 
errors, out of the information positions, that is, moving the nonzero coordinates 
of e out of the information set Ip for Cp, by using a permutation in S subset PAut(Cp) 
(or Phi(S) subset PAut(Cp) if S subset PAut(C)). If S subset PAut(C) subset Sym(n), 
then Phi(S) subset PAut(Cp) subset Sym(np^(s−1)) is computed by using the map Phi 
defined in the manual. The function does not check whether Ip is an information set
for Cp, whether S or Phi(S) is an r-PD-set for Cp with respect to Ip, or whether r ≤ t.

If C is over Z4, function ZpPermutationDecode(C, Ip, S, r, Q) coincides with function 
PermutationDecode(C, I, S, r, Q), which works only for linear codes over Z4. However, 
the former only accepts an information set Ip subset [1,...,2n] for C_2, while the 
latter also accepts an information set I subset [1,...,n] for C as long as the sequence 
S subset PAut(C).
}
    isOverZps, p, s := IsLinearCodeOverZps(C);
    require isOverZps and (s ge 2): "The code C must be over Z/p^s with s>1";
    require PseudoDimension(C) gt 0: "The code C cannot be the zero code";
    requirege r, 1;
    
    k := #Ip;
    n := Length(C);
    np := p^(s-1)*n; 
    require (k ge 1): "Argument 2 cannot be an empty sequence";
    require #C eq p^k: "Argument 2 is not an information set for Argument 1";
    require (Min(Ip) ge 1) and (Max(Ip) le np): 
                                    "Argument 2 should be a subset of", [1..np]; 
                                                                         
    require not IsEmpty(S): "Argument 3 cannot be an empty sequence";   
    numSymG := Degree(Parent(S[1]));  
    require ((numSymG eq n) or (numSymG eq np)): 
          "Argument 3 should contain permutations acting on a set of cardinality", 
                                                                    n, "or", np;
                                                                    
    require not IsEmpty(Q): "Argument 4 cannot be an empty sequence";
    require (BaseRing(Q[1]) eq GF(p)): "Argument 4 must be a sequence of vectors over GF(%o)", p;
    require (OverDimension(Q[1]) eq np): "Argument 4 must contain vectors of length", np;	
    
    if numSymG eq n then
        S := [ PermZpsToPermZp(sigma, p, s) : sigma in S ];
    end if;
 
    grayMapC := CarletGrayMap(C);
    Exclude(~S, Sym(np)!1);

    isDecodedSeq := [];
    uDecodedSeq := [];
    upDecodedSeq := [];
    for u in Q do
        isDecoded, uDecoded, upDecoded := PermutationDecodeAlg(C, S, r, u, Ip, grayMapC);
        Append(~isDecodedSeq, isDecoded);
        Append(~uDecodedSeq, uDecoded);
        Append(~upDecodedSeq, upDecoded);
    end for;

    return isDecodedSeq, uDecodedSeq, upDecodedSeq;

end intrinsic;

/*************************************************************/
intrinsic ZpPermutationDecode(C::CodeLinRng, Ip::[RngIntElt], S::[GrpPermElt], 
          r::RngIntElt, Q::[ModTupRngElt]) -> SeqEnum[BoolElt], SeqEnum[ModTupRngElt],
          SeqEnum[ModTupFldElt]
{
Given a linear code C over Z/p^s of type (n; t1,...,ts), an information set Ip in 
[1..np^(s-1)] for Cp=Phi_s(C) as a sequence of coordinate positions, a sequence S 
such that either S or Phi(S) is an r-PD-set for Cp with respect to the information 
set Ip, an integer t in [1,...,t], where t is the error-correcting capability of Cp, 
and a sequence Q of vectors from the ambient space U=(Z/p^s)^n, the function attempts 
to decode all vectors Phi_s(u), with respect to the code Cp, assuming a systematic        
encoding with respect to the information set Ip. The function returns three values: 
a sequence of booleans representing whether the decoding process have been successful 
for each u in Q, a sequence of codewords of C and a sequence of codewords of Cp. 
For each u in Q, if the decoding algorithm succeeds in computing a codeword u′ in 
Cp as the decoded version of Phi_s(u) in Up, then it returns the preimage of u′ by 
Carlet’s generalized Gray map Phi_s in the second sequence, and u′ in the third 
sequence. If the decoding algorithm does not succeed in decoding u, then the function 
returns the zero codeword in C and the zero codeword in Cp in the corresponding 
positions of the second and third sequences, respectively. 

The permutation decoding algorithm consists in moving all errors in a received vector
u = c + e, where u in Up, c in Cp and e in Up is an error vector with at most t 
errors, out of the information positions, that is, moving the nonzero coordinates 
of e out of the information set Ip for Cp, by using a permutation in S subset PAut(Cp) 
(or Phi(S) subset PAut(Cp) if S subset PAut(C)). If S subset PAut(C) subset Sym(n), 
then Phi(S) subset PAut(Cp) subset Sym(np^(s−1)) is computed by using the map Phi 
defined in the manual. The function does not check whether Ip is an information set
for Cp, whether S or Phi(S) is an r-PD-set for Cp with respect to Ip, or whether r ≤ t.
}
    isOverZps, p, s := IsLinearCodeOverZps(C);
    require isOverZps and (s ge 2): "The code C must be over Z/p^s with s>1";
    require PseudoDimension(C) gt 0: "The code C cannot be the zero code";
    requirege r, 1;
    
    k := #Ip;
    n := Length(C);
    np := p^(s-1)*n; 
    require (k ge 1): "Argument 2 cannot be an empty sequence";
    require #C eq p^k: "Argument 2 is not an information set for Argument 1";
    require (Min(Ip) ge 1) and (Max(Ip) le np): 
                                    "Argument 2 should be a subset of", [1..np]; 
                                                                         
    require not IsEmpty(S): "Argument 3 cannot be an empty sequence";   
    numSymG := Degree(Parent(S[1]));  
    require ((numSymG eq n) or (numSymG eq np)): 
          "Argument 3 should contain permutations acting on a set of cardinality", 
                                                                    n, "or", np;

    Zps := Integers(p^s);
    require not IsEmpty(Q): "Argument 4 cannot be an empty sequence";
    require (BaseRing(Q[1]) eq Zps): "Argument 4 must be a sequence of vectors over Z/%o", p^s;
    require (OverDimension(Q[1]) eq n): "Argument 4 must contain vectors of length", n;	
    
    if numSymG eq n then
        S := [ PermZpsToPermZp(sigma, p, s) : sigma in S ];
    end if;
 
    grayMapC := CarletGrayMap(C);
    Exclude(~S, Sym(np)!1);

    isDecodedSeq := [];
    uDecodedSeq := [];
    upDecodedSeq := [];
    grayMapV := CarletGrayMap(UniverseCode(Zps, n));
    for u in Q do
        isDecoded, uDecoded, upDecoded := PermutationDecodeAlg(C, S, r, grayMapV(u), Ip, grayMapC);
        Append(~isDecodedSeq, isDecoded);
        Append(~uDecodedSeq, uDecoded);
        Append(~upDecodedSeq, upDecoded);
    end for;

    return isDecodedSeq, uDecodedSeq, upDecodedSeq;

end intrinsic;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////				SYNDROME DECODING FUNCTIONS                         ///////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/****************************************************************/
/*                                                              */
/* Function name: MinimumHomogeneousWeightCoset                 */
/* Parameters:  C, rep                                          */
/* Function description: Given a linear code C over Z/p^s and   */
/*   a vector of the same ambient space, rep, return a vector   */
/*   of minimum homogeneous weight in the coset C+rep. It also  */
/*   return the homogeneous weight of the vector.               */
/* Input parameters description:                                */
/*   - C : A code over Zps of length n                          */
/*   - rep : A vetor in the ambient space of the code           */
/* Output parameters description:                               */
/*   - A vector of minim homogeneous weight in the coset C+rep  */ 
/*   - The homogeneous weight of the vector                     */
/*                                                              */
/* Function developed by Abdullah Irfan Basheer                 */
/*                                                              */
/* Signature: (<CodeLinRng> C, <ModTupRngElt> v) ->             */
/*                                   ModTupRngElt, RngIntElt    */
/*                                                              */
/****************************************************************/ 
MinimumHomogeneousWeightCoset := function(C, rep)
    //brute force members of Coset and minimize according to Homogeneous Weight
    minWord := rep;
    minWeight := HomogeneousWeight(rep);
    for x in C do
        if HomogeneousWeight(rep + x) lt minWeight then
            minWord := rep + x;
            minWeight := HomogeneousWeight(rep + x);
        end if;
    end for;

    return minWord, minWeight;
end function;

/****************************************************************/
/*                                                              */
/* Function name: SyndromeSpace                                 */
/* Parameters: C                                                */
/* Function description: Given a code C over Z4 of length n and */
/*   type 2^gamma 4^delta, return the Z4-submodule of Z4^(n-    */
/*   delta) isomorphic to Z2^gamma x Z4^(n-gamma-delta) such    */
/*   that the first gamma coordinates are of order two, that is,*/
/*   the space of syndrome vectors for C. The function also     */
/*   returns the (2n-2delta-gamma)-dimensional binary vector    */
/*   space,  which is the space of syndrome vectors for the     */
/*   corresponding binary code Cbin=Phi(C), where Phi is the    */
/*   Gray map. Note that these spaces are computed by using the */
/*   function ZpInformationSpace(C) applied to the dual code of */
/*   C, given by function DualZ4(C).                            */
/* Input parameters description:                                */
/*   - C : A linear code over Z/p^s                             */
/* Output parameters description:                               */
/*   - Z/p^s-submodule of length n-delta                        */
/*   - Vector space over GF(p) of dimension 2n-2delta-gamma     */
/*                                                              */
/* Signature: (<CodeLinRng> C) -> ModTupRng, ModTupFld          */
/*                                                              */
/****************************************************************/
intrinsic ZpSyndromeSpace(C::CodeLinRng) -> ModTupRng, ModTupFld
{
}  
    isOverZps, p, s := IsLinearCodeOverZps(C);
    require isOverZps and (s ge 2): "The code C must be over Z/p^s with s>1";
    
    InfRSpace, InfVSpace := ZpInformationSpace(ZpDual(C));
    return InfRSpace, InfVSpace;

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: SyndromeDecode                                */
/* Parameters:  C, v                                            */
/* Function description: Given a code C over Zps of lenth n,    */
/*    and a vector v from the ambient space of C, attempt to    */
/*    decode v with respect to C.                               */
/* Input parameters description:                                */
/*   - C : A code over Zps of length n                          */
/*   - v : received word to decode.                             */
/* Output parameters description:                               */
/*   - isDecoded : Boolean that is true if decoding succcesful  */
/*   - w : output decoded word in C                             */
/*                                                              */
/* Function developed by Abdullah Irfan Basheer                 */
/*                                                              */
/* Signature: (<CodeLinRng> C, <RngIntElt> i) -> CodeLinFld     */
/*                                                              */
/****************************************************************/ 
intrinsic ZpSyndromeDecode(C::CodeLinRng, v::ModTupRngElt) -> Bool, ModTupRngElt
{
Given a code C over Zps of lenth n, and a vector v from the ambient space of C, attempt to
decode v with respect to C. 
}
    isOverZps, p, s := IsLinearCodeOverZps(C); 
    require isOverZps and (s ge 2): "The code C must be over Z/p^s with s>1";
    
    errorCapability := Floor((MinimumHomogeneousWeight(C) - 1)/2);
    minWord, minWeight := MinimumHomogeneousWeightCoset(C, v);

    if minWeight gt errorCapability then
        return false, v-minWord;
    else
        return true, v-minWord;
    end if;

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: SyndromeDecode                                */
/* Parameters:  C, v                                            */
/* Function description: Given a code C over Zps of lenth n,    */
/*    and a vector v from the ambient space of C, attempt to    */
/*    decode v with respect to C.                               */
/* Input parameters description:                                */
/*   - C : A code over Zps of length n                          */
/*   - v : received word to decode.                             */
/* Output parameters description:                               */
/*   - isDecoded : Boolean that is true if decoding succcesful  */
/*   - w : output decoded word in C                             */
/*                                                              */
/* Function developed by Abdullah Irfan Basheer                 */
/*                                                              */
/* Signature: (<CodeLinRng> C, <RngIntElt> i) -> CodeLinFld     */
/*                                                              */
/****************************************************************/ 
intrinsic ZpSyndromeDecode(C::CodeLinRng, vs::[ModTupRngElt]) -> Bool, ModTupRngElt
{
Given a code C over Zps of lenth n, and a sequence of vectors v from the ambient space of C, attempt to
decode the vectors with respect to C. Returns a sequence of booleans and decoded vectors.
}
    isOverZps, p, s := IsLinearCodeOverZps(C); 
    require isOverZps and (s ge 2): "The code C must be over Z/p^s with s>1";
    
    bools := [];
    sols := [];
    for x in vs do
        decoded, sol := ZpSyndromeDecode(C, x);
        Append(~bools, decoded);
        Append(~sols, sol);
    end for;

    return bools, sols; 

end intrinsic;
