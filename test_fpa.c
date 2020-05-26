
#include <stdlib.h>
#include <stdio.h>
#include<math.h>
#include<float.h>
#define C 134217729 

double rand_double() {                      // 0.0 to 1.0, uniformly distributed random double values 
   return rand()/(double)RAND_MAX;
}

double rand_double2(double x, double y) {   // a to b 
   return (y - x)*rand_double() + x;
}

void VeltkampSplit(double x, double *xh, double *xl) {                                        //unsigned long C = B^(t/2) + 1;  //s=27 & B=2;
double p = C * x;
double q = x - p;
*xh = p + q; 
*xl = x - *xh;
}

void F2Sum(double a, double b,double *s, double *t) {
double z;
*s= a+b;
z = *s - a;
*t= b-z;
return ;
}
void TwoSum(double a, double b,double *s, double *t) {
double z ,a1 ,b1 ,del_a ,del_b;
*s= a+b;
a1= *s-b;
b1= *s-a1;
del_a = a-a1;
del_b = b-b1;
*t = del_a+del_b;
return ;
}

//PRIST doubly compensated summation algorithm
void Summation(double xh1, double xl1, double xh2, double xl2,double *xh3, double *xl3 ){   
  double t,s,u,v,z,k;
  F2Sum(xh1, xh2, &t, &s);
  F2Sum(xl1, xl2, &u, &v);
  z = s+v+u;
  F2Sum(t, z, xh3, xl3);
return ;
}

//Pichat and Neumaier’s summation algorithm
void Summation2(double xh1, double xl1, double xh2, double xl2,double *xh3, double *xl3 ){  
  double t,s,u,v,z,k;
  TwoSum(xh1, xh2, &t, &s);
  TwoSum(xl1, xl2, &u, &v);
  z=u+v+s; 
  F2Sum(t, z,xh3,xl3);
 return ;
}   

// Dekker's & Polynomial Multiplication
void Mult(double xh1, double xl1, double xh2, double xl2,double *xh3, double *xl3) {
double xh_1, xl_1, xh_2, xl_2;
VeltkampSplit(xh1,&xh_1, &xl_1); 
VeltkampSplit(xh2,&xh_2, &xl_2); 
*xh3 = xh1*xh2;
*xl3  =  xh_1*xh_2 - *xh3;
*xl3 +=  xh_1*xl_2;
*xl3 +=  xl_1*xh_2;
*xl3 +=  xl_1*xl_2+ xl1*xh2 + xh1*xl2;
return;
}

void Division(double xh1, double xl1, double xh2, double xl2,double *xh3, double *xl3){
double xh,x_h,xl,x_l,u1,u0,d0,d1;
*xh3 =xh1/xh2;
VeltkampSplit(*xh3, &xh, &xl) ;
VeltkampSplit(xh2, &x_h, &x_l) ;
Mult(xh, xl, x_h, x_l,&u0,&u1); 
*xl3 = xh1- u0;
*xl3 -= u1;
*xl3 += xl1;
*xl3 -= *xh3 * xl2;
*xl3 /= xh2;  
return;
} 


void SquareRoot(double xh1, double xl1,double *xh3, double *xl3) {   //sqrt(x) just 1 input 
double u1,u0,xh,xl,y0,y1;
if (xh1<0)
printf("error \n");
else
y0 = sqrt(xh1);
VeltkampSplit(y0,&xh, &xl) ;
Mult(xh, xl, xh, xl,&u0,&u1); 
y1 = xh1 - u0;
y1 -= u1;
y1 += xl1 ;
if(y0<0) { 
printf("error \n"); }
y1 = (y1 / y0)*0.5;  
F2Sum(y0,y1,xh3,xl3);
} 

 void Dot(double x1[3], double x2[3], double x3[3], double x4[3],double *xh3, double *xl3) {  //as pairs
      double u0, u1, v0, v1, t0, t1 ,t2, t3;
      Mult(x1[0],x2[0],x3[0],x4[0],&u0,&u1); 
      Mult(x1[1],x2[1],x3[1],x4[1],&v0,&v1); 
      Mult(x1[2],x2[2],x3[2],x4[2],&t0,&t1);  
      Summation2(u0,u1,v0,v1,&t2,&t3);
      Summation2(t0,t1,t2,t3,xh3,xl3);
  } 

  void Inverse_SquareRoot(double x1[3], double x2[3], double x3[3], double x4[3],double *xh3, double *xl3) {  //as pairs
      double xh, xl, x_h, x_l, y0, y1, y2, y3, d0, d1, R3_0, R3_1;
      double u0, u1, uu0, uu1, v0, v1, p0, p1, t0, t1 ,t2, t3, t4, t5, t6, t7;

      Mult(x1[0],x2[0],x3[0],x4[0],&u0,&u1); // starting to scalar product 
      Mult(x1[1],x2[1],x3[1],x4[1],&v0,&v1); 
      Mult(x1[2],x2[2],x3[2],x4[2],&t0,&t1);  
      Summation2(u0,u1,v0,v1,&t2,&t3);
      Summation2(t0,t1,t2,t3,&t4,&t5);

     y0 = sqrt(t4);                          //starting to the square root 
     VeltkampSplit(y0,&xh, &xl) ;
     Mult(xh, xl, xh, xl,&p0,&p1); 
     y1 = t4 - p0;
     y1 -= p1;
     y1 += t5 ;
     if(y0<0) { 
     printf("error \n"); }
     y1 = (y1 / y0)*0.5;         // y0 and y1 are pairs of the squareroot operation

     Mult(y0, y1, y0, y1, &uu0,&uu1); 
     Mult(y0, y1, uu0, uu1, &R3_0,&R3_1);    //pairs of R^3

    // divide with 1's pairs with R^3 's pairs
     y2=1;
     *xh3 = y2/R3_0;
     VeltkampSplit(*xh3, &xh, &xl) ;
     VeltkampSplit(y2, &x_h, &x_l) ;
     Mult(xh, xl, x_h, x_l,&t6,&t7); 
     *xl3 = y2 - t6;
     *xl3 -= t7;
     *xl3 += y3;
     *xl3 -= *xh3 * R3_1;
     *xl3 /= R3_0;      
     }

int main(int argc, char **argv)
{
  int i;
   int i,j,a,b,N;
double xh1=0.1777,xl1=0.000293,xh2=0.335593,xl2=0.0001233, x, y, X, Y, X2, Y2 ; 
double x_pair1[] ={0.3,0.5,0.7};         // P1={ah,bh,ch}
double x_pair2[] ={0.003,0.005,0.007};   // P2={al,bl,cl}
double xh3, xl3, xh4, xl4, xh5, xl5, xh8, xl8, xh88, xl88, xh99, xl99, xhh, xll, xhhh, xlll;
double xh6, xl6, xh66, xl66, xh7, xl7, xh77, xl77, xh9, xl9, xh10, xl10, xh, xl, hxh, lxl;
double relative_errorr[i],relative_error_calc[i], relative_error_ord[i], sum, sum2, mult, mult2;

   FILE *kk,*ff;
    kk=fopen("relative_errorSUM.txt","w");
    ff=fopen("relative_errorMULT.txt","w");
    srand((long)time(NULL));                   // initialize rand() - seed generator for pseudo-random numbers 
    printf("RAND_MAX: %i\n",  RAND_MAX);


   x = 2.0; y = 5.0; N = 10;
    double *relative_error = (double*)malloc(512);
    double *relative_error2 = (double*)malloc(512);
    if (relative_error != NULL) {
       if (relative_error2 != NULL) {
    for (i = 0; i < N ; i++) { 
    xl6 = rand_double();
    xl7 = rand_double();
    xh6 =rand_double2(x, y);
    xh7 =rand_double2(x, y);
    Summation2(xh6, xl6, xh7, xl7, &xh66, &xl66);
    Mult(xh6, xl6, xh7, xl7, &xh77, &xl77);
    X= xh6+xl6;
    Y= xh7+xl7;
    X2= xh6*xh7;
    Y2= xh6*xl7+xh7*xl6; 
    relative_error[i]= (double)fabs((xh66+xl66)-(X+Y))/(xh66+xl66);
    relative_error2[i]= (double)fabs((xh77+xl77)-(X2+Y2))/(xh77+xl77);
    fprintf(kk,"%d %e \n",i,relative_error[i]);
    fprintf(ff,"%d %e \n",i,relative_error2[i]);
    }  


/*  double x = 123.4674567; 
    char buf1[MAX]; 
    double y = 122.4324567; 
    char buf2[MAX]; 
  
    gcvt(x, 16, buf1); 
    gcvt(y, 16, buf2); */


        FILE *kk,*ff;
        kk=fopen("test1.txt","w");
        kk=fopen("test2.txt","w");
       
  /*    for (i=0;i<1000;i++) {
      xh1 += i;
      xl1 += i;
      xh2 += i;
      xl2 += i;
      Summation(xh1,xl1,xh2,xl2,&xh3,&xl3);
      fprintf(kk,"%.10f %.10f \n",xh3,xl3);
      i += 0.00002
      } */

    printf(" Summation 1: \n");
        Summation(xh1,xl1,xh2,xl2,&xh3,&xl3);
    printf("xh = %1.16e\n", xh3); 
    printf("xl = %1.16e \n\n", xl3);

    printf(" Summation 2: \n");
        Summation2(xh1,xl1,xh2,xl2,&xh4,&xl4);
    printf("xh = %1.16e\n", xh4); 
    printf("xl = %1.16e \n\n", xl4);

    printf(" without summation algorithm: \n");
    printf("xh = %1.16e \n", xh1+xh2);
    printf("xl = %1.16e \n\n", xl1+xl2);

        printf(" Multiplication : \n");
        Mult(xh1,xl1,xh2,xl2,&xh5,&xl5);
    printf("xh = %1.16e\n", xh5); 
    printf("xl = %1.16e \n\n", xl5); 
        printf(" without multiplication algorithm: \n");
    printf("xh = %1.16e \n", xh1*xh2);
    printf("xl = %1.16e \n\n", xl1*xh2+xl2*xh1);

        printf("Inverse_SquareRoot : \n\n");
        Inverse_SquareRoot(x_pair1,x_pair2,x_pair1,x_pair2,&xh8,&xl8);
    printf("xh = %1.16e\n", xh8); 
    printf("xl = %1.16e\n", xl8); 

return (0) ;
}

