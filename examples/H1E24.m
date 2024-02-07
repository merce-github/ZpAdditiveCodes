print "Example: H1E24";
ei := GetEchoInput();
SetEchoInput(true);

p := 3;
type := [3,0,0];
s := #type;
C := ZpHadamardCode(p, type);
n := Length(C);
np := n*p^(s-1);

U := RSpace(Integers(p^s), n);
Up := VectorSpace(GF(p), np);
grayMap := CarletGrayMap(UniverseCode(Integers(p^s), n));

I, S, Ip, Sp := PDSetZpHadamardCode(p, type);
r := #Sp-1; r;
t1 := type[1];
r eq Floor((p^(s*(t1-1))-t1)/t1);

c := C![6^^n];
cp := grayMap(c);
ep := Up![1^^r,0^^(np-r)];
up := cp + ep;

isDecoded, uDecoded, upDecoded := ZpPermutationDecode(C, Ip, Sp, r, up);
isDecoded;
uDecoded eq c;
upDecoded eq cp;


p := 3;
type := [1,0,1,1];
s := #type;
C := ZpHadamardCode(p, type);
n := Length(C);
np := n*p^(s-1);

SetVerbose("PDSetHadamardFlag", 1);
I, S, Ip, Sp := PDSetZpHadamardCode(p, type);
r := #Sp-1; r;
IsZpPermutationDecodeSet(C, I, S, r);

I, S, Ip, Sp := PDSetZpHadamardCode(p, type : AlgMethod := "Nondeterministic");
r := #Sp-1; r;
IsZpPermutationDecodeSet(C, I, S, r); 


I, S, Ibin, Sbin := PDSetZpHadamardCode(2, [3,0]);
I2, S2, I2bin, S2bin := PDSetHadamardCodeZ4(3, 5);

I eq I2; I;
S eq S2; S;
Ibin eq I2bin; Ibin;
Sbin eq S2bin; Sbin;

SetEchoInput(ei); 