print "Example: H1E21";
ei := GetEchoInput();
SetEchoInput(true);

C := LinearCode<Integers(27), 5 | [[1,2,5,8,6],
                                   [0,3,6,12,21],
                                   [0,0,9,9,18]]>;
V, Vp, f, fp := ZpInformationSpace(C);

I, Ip := ZpInformationSet(C);
I;
Ip;
mapGrayC := CarletGrayMap(C);
mapGrayV := CarletGrayMap(3, ZpType(C));
[fp(vp) : vp in Vp] eq [mapGrayC(f(vp @@ mapGrayV)) : vp in Vp ];
[vp : vp in Vp] eq [ Vector([fp(vp)[i] : i in Ip ]) : vp in Vp ];

i := V![21,12,18];
ip := Vp![2,2,0,1,2,2];
ip eq mapGrayV(i);
SystematicEncoding(C, ip) eq SystematicEncoding(C, i);
c, cp := SystematicEncoding(C, i);
(c eq f(i)) and (cp eq fp(ip));
cp eq mapGrayC(c);
ip eq Vp![cp[j] : j in Ip];

L := Eltseq(i) cat Eltseq(V![13,3,9]);
SystematicEncoding(C, L);
SystematicEncoding(C, L)[1] eq SystematicEncoding(C, i);

M := Matrix(Integers(27), 3, L);
SystematicEncoding(C, M) eq SystematicEncoding(C, L);

SetEchoInput(ei); 