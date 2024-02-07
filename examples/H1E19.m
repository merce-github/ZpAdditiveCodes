print "Example: H1E19";
ei := GetEchoInput();
SetEchoInput(true);

C := LinearCode<Integers(27), 5 | [[1,2,5,8,6],
                                   [0,3,6,12,21],  
                                   [0,0,9,9,18]]>;
V, Vp, f, fp := ZpInformationSpace(C);

I, Ip := ZpInformationSet(C);
I;
Ip;

#PunctureCode(C, {1..5} diff Set(I)) eq #C;
Cp := CarletGrayMapImage(C);
#{Eltseq(cp)[Ip] : cp in Cp} eq #Cp;
 
IsZpInformationSet(C, I);
IsZpInformationSet(C, Ip);
IsZpInformationSet(C, [1,2,3,4,5,6]);

ip := Vp![1,0,0,2,0,2];
cp := fp(ip);
Vp![cp[i] : i in Ip] eq ip;
#[ip : ip in Vp | Vp![fp(ip)[i] : i in Ip] eq ip ] eq #Vp; 

SetEchoInput(ei); 