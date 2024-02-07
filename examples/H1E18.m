print "Example: H1E18";
ei := GetEchoInput();
SetEchoInput(true);

C := LinearCode<Integers(27), 5 | [[1,2,5,8,6],
                                   [0,3,6,12,21],  
                                   [0,0,9,9,18]]>;
V, Vp, f, fp := ZpInformationSpace(C);

(#V eq #C) and (#Vp eq #C);
Set([f(i) : i in V]) eq Set(C);
Set([fp(i) : i in Vp]) eq Set(CarletGrayMapImage(C));

i := V![17,21,18];
c := f(i);
c;
c in C;

ip := Vp![1,0,0,2,0,2];
cp := fp(ip);
cp;
cp in CarletGrayMapImage(C);

mapGrayC := CarletGrayMap(C);
mapGrayV := CarletGrayMap(3, ZpType(C));
Set(Vp) eq {mapGrayV(v) : v in V};
ip eq mapGrayV(i);
cp eq mapGrayC(c);
[mapGrayC(f(v)) : v in V] eq [fp(mapGrayV(v)) : v in V];

SetEchoInput(ei); 