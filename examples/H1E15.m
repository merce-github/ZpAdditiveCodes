print "Example: H1E15";
ei := GetEchoInput();
SetEchoInput(true);

C := LinearCode<Integers(27), 4 | [[1,1,2,20],
                                   [0,3,0,12],
                                   [0,0,9,18]]>;
L := CosetRepresentatives(C);
Set(RSpace(Integers(27),4)) eq {v+ci : v in Set(C), ci in L};

 
S := LinearCode<Integers(27), 4 | [[1,1,2,20],
                                   [0,3,0,12]]>; 
S subset C;
L := CosetRepresentatives(C, S);
Set(C) eq {v+ci : v in Set(S), ci in L};

SetEchoInput(ei); 