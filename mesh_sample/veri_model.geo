// -------------------------------------
lc=0.02; // element size 

rc=0.5;  // radius of circle

xc=0.0;   
yc=0.0;  // center of circle (xc,yc)

d=1.5;
// -------------------------------------

Point(1)={xc   ,yc   ,0.0,lc}; // center of circle
Point(2)={xc+rc,yc   ,0.0,lc};
Point(3)={xc   ,yc+rc,0.0,lc};
Point(4)={xc-rc,yc   ,0.0,lc};
Point(5)={xc   ,yc-rc,0.0,lc};

Circle(1)={2,1,3};
Circle(2)={3,1,4};
Circle(3)={4,1,5};
Circle(4)={5,1,2};

o1[]=Translate { 1.0*d, 0, 0} { Duplicata{ Line{1,2,3,4}; } };
o2[]=Translate {-1.0*d, 0, 0} { Duplicata{ Line{1,2,3,4}; } };
o3[]=Translate { 2.0*d, 0, 0} { Duplicata{ Line{1,2,3,4}; } };
o4[]=Translate {-2.0*d, 0, 0} { Duplicata{ Line{1,2,3,4}; } };
o5[]=Translate { 3.0*d, 0, 0} { Duplicata{ Line{1,2,3,4}; } };
o6[]=Translate {-3.0*d, 0, 0} { Duplicata{ Line{1,2,3,4}; } };
o7[]=Translate { 4.0*d, 0, 0} { Duplicata{ Line{1,2,3,4}; } };
o8[]=Translate {-4.0*d, 0, 0} { Duplicata{ Line{1,2,3,4}; } };
o9[]=Translate { 5.0*d, 0, 0} { Duplicata{ Line{1,2,3,4}; } };
o0[]=Translate {-5.0*d, 0, 0} { Duplicata{ Line{1,2,3,4}; } };
p1[]=Translate { 6.0*d, 0, 0} { Duplicata{ Line{1,2,3,4}; } };
p2[]=Translate {-6.0*d, 0, 0} { Duplicata{ Line{1,2,3,4}; } };
p3[]=Translate { 7.0*d, 0, 0} { Duplicata{ Line{1,2,3,4}; } };
p4[]=Translate {-7.0*d, 0, 0} { Duplicata{ Line{1,2,3,4}; } };
p5[]=Translate { 8.0*d, 0, 0} { Duplicata{ Line{1,2,3,4}; } };
p6[]=Translate {-8.0*d, 0, 0} { Duplicata{ Line{1,2,3,4}; } };

Physical Line(99)={-1,-2,-3,-4,
                   -o1[0],-o1[1],-o1[2],-o1[3],-o2[0],-o2[1],-o2[2],-o2[3],
                   -o3[0],-o3[1],-o3[2],-o3[3],-o4[0],-o4[1],-o4[2],-o4[3],
                   -o5[0],-o5[1],-o5[2],-o5[3],-o6[0],-o6[1],-o6[2],-o6[3],
                   -o7[0],-o7[1],-o7[2],-o7[3],-o8[0],-o8[1],-o8[2],-o8[3],
                   -o9[0],-o9[1],-o9[2],-o9[3],-o0[0],-o0[1],-o0[2],-o0[3],
                   -p1[0],-p1[1],-p1[2],-p1[3],-p2[0],-p2[1],-p2[2],-p2[3],
                   -p3[0],-p3[1],-p3[2],-p3[3],-p4[0],-p4[1],-p4[2],-p4[3],
                   -p5[0],-p5[1],-p5[2],-p5[3],-p6[0],-p6[1],-p6[2],-p6[3]}; // open region
Physical Line( 1)={ 1, 2, 3, 4,
                    o1[0], o1[1], o1[2], o1[3], o2[0], o2[1], o2[2], o2[3],
                    o3[0], o3[1], o3[2], o3[3], o4[0], o4[1], o4[2], o4[3],
                    o5[0], o5[1], o5[2], o5[3], o6[0], o6[1], o6[2], o6[3],
                    o7[0], o7[1], o7[2], o7[3], o8[0], o8[1], o8[2], o8[3],
                    o9[0], o9[1], o9[2], o9[3], o0[0], o0[1], o0[2], o0[3],
                    p1[0], p1[1], p1[2], p1[3], p2[0], p2[1], p2[2], p2[3],
                    p3[0], p3[1], p3[2], p3[3], p4[0], p4[1], p4[2], p4[3],
                    p5[0], p5[1], p5[2], p5[3], p6[0], p6[1], p6[2], p6[3]}; // domain 1

