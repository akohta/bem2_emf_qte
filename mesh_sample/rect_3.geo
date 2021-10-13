//--------------------
lc0=0.005; // element size 
lc1=0.020;  // element size

lx=0.5-lc0*0.1;  // width of x-direction
ly=1.0;  // width of y-direction
lt=0.95;// thickness 

xc=0.0;  //
yc=0.5;  // center of rectangle
// -------------------

Point(1)={xc+0.5*lx,yc+0.5*ly,0.0,lc1}; 
Point(2)={xc-0.5*lx,yc+0.5*ly,0.0,lc1};
Point(3)={xc-0.5*lx,yc-0.5*ly,0.0,lc0};
Point(4)={xc+0.5*lx,yc-0.5*ly,0.0,lc0};
Point(5)={xc+0.5*lx,yc+0.5*ly-lt,0.0,lc0}; 
Point(6)={xc-0.5*lx,yc+0.5*ly-lt,0.0,lc0};

Line(1)={1,2};
Line(2)={2,6};
Line(3)={6,5};
Line(4)={5,1};
Line(5)={6,3};
Line(6)={3,4};
Line(7)={4,5};

Physical Line(99)={-1,-2,-5,-6,-7,-4}; // open region
Physical Line( 1)={ 1, 2, 3, 4}; // domain 1
Physical Line( 2)={-3, 5, 6, 7}; // domain 2

