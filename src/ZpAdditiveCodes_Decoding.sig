175,0
S,ZpAdditiveCodes_Decoding_version,Return the current version of this package,0,0,0,0,0,0,0,82,-38,-38,-38,-38,-38
V,IsPDSetFlag,3
V,PDSetHadamardFlag,2
S,IsZpPermutationDecodeSet,"Given a linear code C over Z/p^s of type (n; t1,...,ts), a sequence I in [1,...,np^(s-1)], a sequence S of elements in the symmetric group Sym(np^(s-1)) of permutations on the set [1,...,np^(s-1)], and an integer r>=1, return true if and only if S is a r-PD-set for Cp=Phi_s(C), where Phi_s is Carlet's generalized Gray map, with respect to the information set I. The parameters I and S can also be given as a sequence I in [1,...,n] and a sequence S of elements in the symmetric group Sym(n) of permutations on the set [1,...,n], respectively. In this case, the function returns true if and only if Phi(S) is a r-PD-set for Cp=Phi_s(C) with respect to the information set Phi(I), where Phi(I) and Phi(S) are the sequences defined as in the manual. Depending on the length of the code C, its type, and the integer r, this function could take some time to compute whether S or Phi(S) is a r-PD-set for Cp with respect to I or Phi(I), respectively. Specifically, if the function returns true, it is necessary to check sum_(i=1)^r (|I| choose i)*((N-|I|) choose (r-i)) r-sets, where N=n and |I|=t1+···+ts when I is given as an information set for C, or N=np^(s-1) and |I|=s*t1+···+ts when I is given as an information set for Cp=Phi_s(C). The verbose flag IsPDsetFlag is set to level 0 by default. If it is set to level 1, the total time used to check the condition is shown. Moreover, the reason why the function returns false is also shown, that is, whether I is not an information set, S is not a subset of the permutation automorphism group or S is not an r-PD-set. If it is set to level 2, the percentage of the computation process performed is also printed. If C is over Z4, function IsZpPermutationDecodeSet(C, I, S, r) coincides with function IsPermutationDecodeSet(C, I, S, r), which works only for linear codes over Z4, but the former may perform less efficiently when I subset [1,...,2n] because it calls function IsZpInformationSet(C, I) instead of IsInformationSet(C, I)",2,1,1,82,0,148,2,1,82,0,222,4,0,0,0,0,0,0,0,148,,0,0,82,,0,0,82,,0,0,202,,36,-38,-38,-38,-38,-38
S,PDSetZpHadamardCode,"Given a prime number p and a sequence of nonnegative integers [t1,...,ts] with t1>0 and s>1, the generalized Hadamard code C over Z/p^s of type (n; t1,...,ts), given by function ZpHadamardCode(p, [t1,...,ts]), is considered. The function returns an information set I in [1,...,n] for C together with a subset S of the permutation automorphism group of C such that Phi(S) is an r-PD-set for C_p=Phi_s(C) with respect to Phi(I), where Phi_s is the Carlet's generalized Gray map as defined above, Phi(I) and Phi(S) as defined in the manual. The function also returns the information set Phi(I) and the r-PD-set Phi(S). Note that for p=2 and [t1,...,ts] = [1,0,...,0,ts] or [t1,...,ts] = [1,0,...,0,1,ts], we have that C_p is linear, so it is possible to find an r-PD-set of size r+1 for C_p, for any r ≤ floor(2^m/(m+1)), by using function PDSetHadamardCode(m) with m = ts + s if [t1,...,ts] = [1,0,...,0,ts] and m = ts + s + 2 if [t1,..., ts] = [1,0,...,0,1,ts]. The information sets I and Phi(I) are returned as sequences of t1+···+ts and s(t1-1)+(s-1)t2+···+ts integers, giving the coordinate positions that correspond to the information sets for C and C_p, respectively. The sets S and Phi(S) are also returned as sequences of elements in the symmetric groups Sym(n) and Sym(p^(s-1)*n) of permutations on the set [1,...,n] and [1,...,p^(s-1)*n], respectively. A deterministic algorithm is used by default. In this case, when t1>=2, the function first computes r as the maximum of g_p^(t1,...,ts) and ftilda_p^(t1,...,ts), where g_p^(t1,...,ts) is given by the general construction descrided in reference [1] and ftilda_p^(t1,...,ts)=max(f_p^(t1,0,..,0), f_p^(1,t2,0...0), ..., f_p^(1,0...0,ts)) <= f_p^(t1,...,ts), where f_p^(t1,...,ts) is the is the theoretical upper bound for the value of r such that there exists an r-PD-set of size r+1 for a Hadamard code over Z/p^s of type (n; t1,...,ts), given in reference [2]. Let i be the first index i ≥ 1 such that the maximum ftilda_p^(t1,...,ts) is achieved. Then, if ftilda_p^(t1,...,ts) > g_p^(t1,...,ts), it constructs an r-PD-set of size r+1 for the generalized Hadamard code over Z/p^s of type (n; t1,0,...,0) if i=1 or (n; 1,0,...,0,ti,0,...,0) otherwise, which is transformed into an r-PD-set for C, by using the recursive construction defined in the reference. The value of r remains unchanged after the recursive construction, so r = ftilda_p^(t1,...,ts). On the other hand, if g_p^(t1,...,ts) > ftilda_p^(t1,...,ts), the general construction is applied and an r-PD-set of size r+1 with r = g_p^(t1,...,ts) is obtained. When t1=1, and there exists an index i>=2 such that ti<>0 and t2=...=t_(i-1)=0, it first constructs an r-PD-set of size r+1 for the generalized Hadamard code over Z/p^s of type (n; 1+ti, t_(i+1),...,ts), which satisfies that 1+ti >= 2, by following the same process as described above when t1 ≥ 2. Then, the obtained r-PD-set is transformed into an r-PD-set for C as described in references [1,2]. If the parameter AlgMethod is assigned the value ""Nondeterministic"", the function tries to improve the previous results by finding an r-PD-set of size r+1 such that ftilda_p^(t1,...,ts) <= r <= f_p^(t1,...,ts). In this case, the function starts from the maximum value of r = f_p^(t1,...,ts) and attempts to find an r-PD-set within a time limit of 5.0 seconds of “user time”. This is performed 5 times, each time starting from an empty set and trying to construct the r-PD set incrementally, by adding elements randomly. If such an r-PD-set is not found, the value of r decreases by one and the same process takes place with the new value of r. The value of r keeps decreasing until an r-PD-set is found or r = ftilda_p^(t1,...,ts). The verbose flag PDsetHadamardFlag is set to level 0 by default. If it is set to level 1, some information about the process of constructing the r-PD-set of size r+1 is shown. Moreover, the value of the theoretical upper bound f_p^(t1,...,ts) given in reference [1] is also shown. If p = 2 and s = 2, then function PDSetHadamardCodeZ4(t1, 2t1 + t2 − 1) can also be used. Both function only coincide when t2 = 0. When t2 > 0, the output parameters I and Phi(I) coincide as sets and PDSetZpHadamardCode(2, [t1, t2]) may give a larger r-PD-set",1,1,1,82,0,148,2,0,0,0,0,0,0,0,82,,0,0,148,,82,82,82,82,82,-38
S,ZpPermutationDecode,"Given a linear code C over Z/p^s of type (n; t1,...,ts), an information set Ip in [1..np^(s-1)] for Cp=Phi_s(C) as a sequence of coordinate positions, a sequence S such that either S or Phi(S) is an r-PD-set for Cp with respect to the information set Ip, an integer t in [1,...,t], where t is the error-correcting capability of Cp, and a vector u from the ambient space Up=Zp^(np^(s-1)), the function attempts to decode u in Up with respect to the code Cp, assuming a systematic encoding with respect to the information set Ip. If the decoding algorithm succeeds in computing a codeword u' in Cp as the decoded version of u in Up, then the function returns true, the preimage of u' by Carlet's generalized Gray map Phi_s and finally u'. If the decoding algorithm does not succeed in decoding u, then the function returns false, the zero codeword in C and the zero codeword in Cp. The permutation decoding algorithm consists in moving all errors in a received vector u = c + e, where u in Up, c in Cp and e in Up is an error vector with at most t errors, out of the information positions, that is, moving the nonzero coordinates of e out of the information set Ip for Cp, by using a permutation in S subset PAut(Cp) (or Phi(S) subset PAut(Cp) if S subset PAut(C)). If S subset PAut(C) subset Sym(n), then Phi(S) subset PAut(Cp) subset Sym(np^(s−1)) is computed by using the map Phi defined in the manual. The function does not check whether Ip is an information set for Cp, whether S or Phi(S) is an r-PD-set for Cp with respect to Ip, or whether r ≤ t. If C is over Z4, function ZpPermutationDecode(C, Ip, S, r, u) coincides with function PermutationDecode(C, I, S, r, u), which works only for linear codes over Z4. However, the former only accepts an information set Ip subset [1,...,2n] for C_2, while the latter also accepts an information set I subset [1,..,n] for C as long as the sequence S subset PAut(C)",2,1,1,82,0,148,2,1,82,0,222,5,0,0,0,0,0,0,0,160,,0,0,148,,0,0,82,,0,0,82,,0,0,202,,36,192,160,-38,-38,-38
S,ZpPermutationDecode,"Given a linear code C over Z/p^s of type (n; t1,...,ts), an information set Ip in [1..np^(s-1)] for Cp=Phi_s(C) as a sequence of coordinate positions, a sequence S such that either S or Phi(S) is an r-PD-set for Cp with respect to the information set Ip, an integer t in [1,...,t], where t is the error-correcting capability of Cp, and a vector u from the ambient space U=(Z/p^s)^n, the function attempts to decode Phi_s(u) with respect to the code Cp, assuming a systematic encoding with respect to the information set Ip. If the decoding algorithm succeeds in computing a codeword u' in Cp as the decoded version of Phi_s(u) in Up=Zp^(np^(s-1)), then the function returns true, the preimage of u' by Carlet's generalized Gray map Phi_s and finally u'. If the decoding algorithm does not succeed in decoding u, then the function returns false, the zero codeword in C and the zero codeword in Cp. The permutation decoding algorithm consists in moving all errors in a received vector u = c + e, where u in Up, c in Cp and e in Up is an error vector with at most t errors, out of the information positions, that is, moving the nonzero coordinates of e out of the information set Ip for Cp, by using a permutation in S subset PAut(Cp) (or Phi(S) subset PAut(Cp) if S subset PAut(C)). If S subset PAut(C) subset Sym(n), then Phi(S) subset PAut(Cp) subset Sym(np^(s−1)) is computed by using the map Phi defined in the manual. The function does not check whether Ip is an information set for Cp, whether S or Phi(S) is an r-PD-set for Cp with respect to Ip, or whether r ≤ t. If C is over Z4, function ZpPermutationDecode(C, Ip, S, r, u) coincides with function PermutationDecode(C, I, S, r, u), which works only for linear codes over Z4. However, the former only accepts an information set Ip subset [1,...,2n] for C_2, while the latter also accepts an information set I subset [1,..,n] for C as long as the sequence S subset PAut(C)",2,1,1,82,0,148,2,1,82,0,222,5,0,0,0,0,0,0,0,192,,0,0,148,,0,0,82,,0,0,82,,0,0,202,,36,192,160,-38,-38,-38
S,ZpPermutationDecode,"Given a linear code C over Z/p^s of type (n; t1,...,ts), an information set Ip in [1..np^(s-1)] for Cp=Phi_s(C) as a sequence of coordinate positions, a sequence S such that either S or Phi(S) is an r-PD-set for Cp with respect to the information set Ip, an integer t in [1,...,t], where t is the error-correcting capability of Cp, and a sequence Q of vectors from the ambient space Up=Zp^(np^(s-1)), the function attempts to decode all vectors u in Q subset Up, with respect to the code Cp, assuming a systematic encoding with respect to the information set Ip. The function returns three values: a sequence of booleans representing whether the decoding process have been successful for each u in Q, a sequence of codewords of C and a sequence of codewords of Cp. For each u in Q, if the decoding algorithm succeeds in computing a codeword u′ in Cp as the decoded version of u in Up, then it returns the preimage of u′ by Carlet’s generalized Gray map Phi_s in the second sequence, and u′ in the third sequence. If the decoding algorithm does not succeed in decoding u, then the function returns the zero codeword in C and the zero codeword in Cp in the corresponding positions of the second and third sequences, respectively. The permutation decoding algorithm consists in moving all errors in a received vector u = c + e, where u in Up, c in Cp and e in Up is an error vector with at most t errors, out of the information positions, that is, moving the nonzero coordinates of e out of the information set Ip for Cp, by using a permutation in S subset PAut(Cp) (or Phi(S) subset PAut(Cp) if S subset PAut(C)). If S subset PAut(C) subset Sym(n), then Phi(S) subset PAut(Cp) subset Sym(np^(s−1)) is computed by using the map Phi defined in the manual. The function does not check whether Ip is an information set for Cp, whether S or Phi(S) is an r-PD-set for Cp with respect to Ip, or whether r ≤ t. If C is over Z4, function ZpPermutationDecode(C, Ip, S, r, Q) coincides with function PermutationDecode(C, I, S, r, Q), which works only for linear codes over Z4. However, the former only accepts an information set Ip subset [1,...,2n] for C_2, while the latter also accepts an information set I subset [1,...,n] for C as long as the sequence S subset PAut(C)",3,1,1,82,0,148,2,1,82,0,222,4,1,82,0,160,5,0,0,0,0,0,0,0,82,,0,0,148,,0,0,82,,0,0,82,,0,0,202,,82,82,82,-38,-38,-38
S,ZpPermutationDecode,"Given a linear code C over Z/p^s of type (n; t1,...,ts), an information set Ip in [1..np^(s-1)] for Cp=Phi_s(C) as a sequence of coordinate positions, a sequence S such that either S or Phi(S) is an r-PD-set for Cp with respect to the information set Ip, an integer t in [1,...,t], where t is the error-correcting capability of Cp, and a sequence Q of vectors from the ambient space U=(Z/p^s)^n, the function attempts to decode all vectors Phi_s(u), with respect to the code Cp, assuming a systematic encoding with respect to the information set Ip. The function returns three values: a sequence of booleans representing whether the decoding process have been successful for each u in Q, a sequence of codewords of C and a sequence of codewords of Cp. For each u in Q, if the decoding algorithm succeeds in computing a codeword u′ in Cp as the decoded version of Phi_s(u) in Up, then it returns the preimage of u′ by Carlet’s generalized Gray map Phi_s in the second sequence, and u′ in the third sequence. If the decoding algorithm does not succeed in decoding u, then the function returns the zero codeword in C and the zero codeword in Cp in the corresponding positions of the second and third sequences, respectively. The permutation decoding algorithm consists in moving all errors in a received vector u = c + e, where u in Up, c in Cp and e in Up is an error vector with at most t errors, out of the information positions, that is, moving the nonzero coordinates of e out of the information set Ip for Cp, by using a permutation in S subset PAut(Cp) (or Phi(S) subset PAut(Cp) if S subset PAut(C)). If S subset PAut(C) subset Sym(n), then Phi(S) subset PAut(Cp) subset Sym(np^(s−1)) is computed by using the map Phi defined in the manual. The function does not check whether Ip is an information set for Cp, whether S or Phi(S) is an r-PD-set for Cp with respect to Ip, or whether r ≤ t",3,1,1,82,0,148,2,1,82,0,222,4,1,82,0,192,5,0,0,0,0,0,0,0,82,,0,0,148,,0,0,82,,0,0,82,,0,0,202,,82,82,82,-38,-38,-38
S,ZpSyndromeSpace,,0,1,0,0,0,0,0,0,0,202,,191,159,-38,-38,-38,-38
S,ZpSyndromeDecode,"Given a code C over Zps of lenth n, and a vector v from the ambient space of C, attempt to decode v with respect to C",0,2,0,0,0,0,0,0,0,192,,0,0,202,,37,192,-38,-38,-38,-38
S,ZpSyndromeDecode,"Given a code C over Zps of lenth n, and a sequence of vectors v from the ambient space of C, attempt to decode the vectors with respect to C. Returns a sequence of booleans and decoded vectors",1,1,1,82,0,192,2,0,0,0,0,0,0,0,82,,0,0,202,,37,192,-38,-38,-38,-38