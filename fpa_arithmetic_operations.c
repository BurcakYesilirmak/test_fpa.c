
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<float.h>
#include<time.h>
#include<limits.h>
#define C  134217729 

//define C2 4294967297
//unsigned long C = B^(t/2) + 1;  
//s=27 & B=2;

void VeltkampSplit(double x, double *xh, double *xl) {        
double p = C * x;
double q = x - p;
*xh = p + q; 
*xl = x - *xh;
}

void F2Sum(double a, double b,double *xh, double *xl) {
double z;
*xh= a+b;
z = *xh-a;
*xl= b-z;
}

void TwoSum(double a, double b,double *xh, double *xl) {
double a1 ,b1 ,del_a ,del_b;
*xh= a+b;
a1= *xh-b;
b1= *xh-a1;
del_a = a-a1;
del_b = b-b1;
*xl = del_a+del_b;
}

//PRIST doubly compensated summation algorithm
void Summation(double xh1, double xl1, double xh2, double xl2,double *xh3, double *xl3 ){   
  double t,s,u,v,z;
  F2Sum(xh1, xh2, &t, &s);
  F2Sum(xl1, xl2, &u, &v);
  z = v+s+u; 
  F2Sum(t, z,xh3,xl3);
  /* *xh3=t;
     *xl3=z; */
}

//Pichat and Neumaier’s summation algorithm
void Summation2(double xh1, double xl1, double xh2, double xl2,double *xh3, double *xl3 ){  
  double t,s,u,v,z;
  TwoSum(xh1, xh2, &t, &s);
  TwoSum(xl1, xl2, &u, &v);
  z= v+s+u; 
  //F2Sum(t, z, xh3, xl3);
  *xh3=t;
  *xl3=z; 
}   

// Dekker's & Polynomial Multiplication
/* void Mult(double xh1, double xl1, double xh2, double xl2,double *xh3, double *xl3) {
double xh_1, xl_1, xh_2, xl_2, t, z;
VeltkampSplit(xh1, &xh_1, &xl_1); 
VeltkampSplit(xh2, &xh_2, &xl_2); 
*xh3 =  xh1*xh2;
*xl3 =  xh_1*xh_2 - *xh3;
*xl3 +=  xh_1*xl_2;
*xl3 +=  xl_1*xh_2;
*xl3 +=  xl_1*xl_2+ xl1*xh2 + xh1*xl2;
return;
} */

    // Dekker's & Polynomial Multiplication
  void Mult(double xh1, double xl1, double xh2, double xl2,double *t, double *z) {
    double xh_1, xl_1, xh_2, xl_2, xh3, xl3;
    VeltkampSplit(xh1, &xh_1, &xl_1); 
    VeltkampSplit(xh2, &xh_2, &xl_2); 
    xh3 =  xh1*xh2;
    xl3 =  xh_1*xh_2 - xh3;
    xl3 +=  xh_1*xl_2;
    xl3 +=  xl_1*xh_2;
    xl3 +=  xl_1*xl_2+ xl1*xh2 + xh1*xl2;
    F2Sum(xh3, xl3, t, z);
    return;
    }

// inverse multiplication
void Division2(double xh1, double xl1, double xh_inv, double xl_inv,double *xh3, double *xl3) {
  double xh_1, xl_1, xh_2, xl_2, t, z;
  double xh2,xl2;
  if ( (xh_inv != 0) || (xl_inv != 0) ) {  //check order of operations seperately
  xh2= 1/xh_inv;
  xl2= 1/xl_inv; }
  else {  
  xh2=xh_inv; 
  xl2=xl_inv; } 
//  xh2= 1/xh_inv;
//  xl2= xl_inv;
  Mult(xh1, xl1, xh2, xl2, xh3, xl3);
/* VeltkampSplit(xh1, &xh_1, &xl_1); 
VeltkampSplit(xh2, &xh_2, &xl_2); 
*xh3 =  xh1*xh2;
*xl3 =  xh_1*xh_2 - *xh3;
*xl3 +=  xh_1*xl_2;
*xl3 +=  xl_1*xh_2;
*xl3 +=  xl_1*xl_2+ xl1*xh2 + xh1*xl2; */
}

// with iteration
void Division(double xh1, double xl1, double xh2, double xl2,double *xh3, double *xl3){
double u1, u0, t, z;
t =xh1/xh2;
Mult(t, 0, xh2, 0,&u0,&u1); 
z = xh1- u0;
z -= u1;
z += xl1;
z -= t * xl2;
z /= xh2;  
F2Sum(t, z, xh3, xl3);
return;
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
     y1 = (y1 / y0)*0.5;         // y0 and y1 are pairs of the square root operation

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

   void unit_vector_longdouble(long double xh1, long double xh2, long double xh3, long double *X_unit, long double *Y_unit, long double *Z_unit) {
   long double alpha, F, dF, eps , scalar_new, max ,Z_new, Y_new, X_new ;
   int count,i;
   long double k[]={ xh1, xh2, xh3 }; 
    max=fabsl(k[1]);
    for(i=0;i<3;i++) {
    if(fabsl(k[i])>max) 
    max=fabsl(k[i]);
    }
    //initilial value of alpha long double
    alpha = 1/max;
    eps=1e-15;
    count =0;
    do {
        count ++ ;
        //for double 
       X_new= xh1*alpha;
       Y_new= xh2*alpha;
       Z_new= xh3*alpha;
       scalar_new = X_new*X_new + Y_new*Y_new + Z_new*Z_new ;
       F = scalar_new-1 ;
       dF =  2/alpha*scalar_new ;
       alpha = (F/dF)*-1 +alpha; } 
     while (scalar_new - 1 >= eps);
    *X_unit = xh1*alpha;
    *Y_unit = xh2*alpha;
    *Z_unit = xh3*alpha;   
        printf("scalar product long double= %Le \n",scalar_new);
    }


   void unit_vector(double xh1, double xh2, double xh3, double *X_unit, double *Y_unit, double *Z_unit) {
        
        double alpha, F, dF, eps , scalar_new, max ,Z_new, Y_new, X_new ;
        int count,i;
        long double X2, Y2, Z2, I_long, J_long, K_long;
        X2= xh1;
        Y2= xh2;
        Z2= xh3; 

    //finding maximum absolute value of the vector
    double k[]={ xh1, xh2, xh3 };  
    max=fabs(k[1]);
    for(i=0;i<3;i++) {
    if(fabs(k[i])>max) 
    max=fabs(k[i]);
    }

    //initilial value of alpha with double
    alpha = 1/max;
    eps=1e-15;
    do {
        count ++ ;
        //for double 
       X_new= xh1*alpha;
       Y_new= xh2*alpha;
       Z_new= xh3*alpha;
       scalar_new = X_new*X_new + Y_new*Y_new + Z_new*Z_new ;
       F = scalar_new-1 ;
       dF =  2/alpha*scalar_new ;
       alpha = (F/dF)*-1 +alpha; } 
     while (scalar_new - 1 >= eps);
    *X_unit = xh1*alpha;
    *Y_unit = xh2*alpha;
    *Z_unit = xh3*alpha;  
        printf("scalar product for double= %le \n",scalar_new);
    }  

   void scalar_product(double xh[3], double xl[3], double *xh3, double *xl3) {  // the square of the vector itself
      double u0, u1, v0, v1, t0, t1 ,t2, t3;
      Mult(xh[0], xl[0], xh[0], xl[0], &u0, &u1); 
      Mult(xh[1], xl[1], xh[1], xl[1], &v0, &v1); 
      Mult(xh[2], xl[2], xh[2], xl[2], &t0, &t1);  
      Summation2(u0,u1,v0,v1,&t2,&t3);
      Summation2(t0,t1,t2,t3,xh3,xl3);
  }

   void unit_vector_pairwise(double xh[3], double xl[3], double *x_h, double *x_l) {
        double alpha, F, dF, eps , scalar_new, max ,Z_new, Y_new, X_new ;
        double xh11, xl11, xh22, xl22, xh33, xl33, yhh, yll, zh33, zl33, fhh, fll, dfhh, dfll;
        double xh44, xl44, xh55, xl55, xh66, xl66, xh555, xl555, xh666, xl666, xh333, xl333, xh444, xl444, xh222, xl222;
        int i,count;

    //finding maximum absolute value of the vector
    double k[]={ xh[0], xh[1], xh[2] };  
    max=fabs(k[1]);
    for(i=0;i<3;i++) {
    if(fabs(k[i])>max) 
    max=fabs(k[i]);
    }
    //initilial value of alpha with double
    alpha = 1/max;
    eps=1e-15;
    do {
        count ++ ;
        //for double 
       X_new= xh[0]*alpha;
       Y_new= xh[1]*alpha;
       Z_new= xh[2]*alpha;
       scalar_new = X_new*X_new + Y_new*Y_new + Z_new*Z_new ;
       F = scalar_new-1 ;
       dF =  2/alpha*scalar_new ;
       alpha -= F/dF;} 
     while (scalar_new - 1 >= eps);    // get approximate good alpha 

         //  calculating alpha biggest part as pairs
        Mult(alpha, 0.0, xh[0], xl[0], &xh44, &xl44);
        Mult(alpha, 0.0, xh[1], xl[1], &xh55, &xl55);
        Mult(alpha, 0.0, xh[2], xl[2], &xh66, &xl66); 
       double xx[3]= {xh44, xh55, xh66};
       double yy[3]= {xl44, xl55, xl66};
       scalar_product(xx, yy, &xh33, &xl33); 

       // alpha = alpha - (alpha*0.5) * (sc-1)/sc  
       Mult(alpha, 0.0, 0.5, 0.0, &yhh, &yll);      // alpha*0.5                     
       Division(xh33-1, xl33, xh33, xl33, &dfhh, &dfll);        // (sc-1)/sc
       Mult(dfhh, dfll, yhh, yll, &xh555, &xl555);       
       Mult(xh555, xl555, -1, 0.0, &xh22, &xl22);
       Summation2(xh22, xl22, alpha, 0.0, &zh33, &zl33); 

     // xh[]*alpha
    Mult(xh[0], 0.0, zh33, zl33, &xh444, &xl444);
    Mult(xh[1], 0.0, zh33, zl33, &xh222, &xl222); 
    Mult(xh[2], 0.0, zh33, zl33, &xh333, &xl333);   
    x_h[0] =  xh444;   // *(x_h+0)=x_h[0]
    x_h[1] =  xh222;
    x_h[2] =  xh333; 
    x_l[0] =  xl444;
    x_l[1] =  xl222;
    x_l[2] =  xl333; 
   }


int main(int argc, char **argv)
{
int i;

//double xh1=0.37797,xl1=0.000293,xh2=0.335593,xl2=0.0001233, xh3=0.3355674, I, J, K;
double x1= -10.0, x2= 6.0, x3= 4.0 ;
long double x1long= -10.0, x2long= 6.0, x3long= 4.0 ;

double xh_pair[]={-10.0, 6.0, 4.0}, xh_p[3];
double xl_pair[]={0, 0, 0}, xl_p[3];

double I, J, K;
long double T, P, L;

double x_pair1[] ={0.3,0.5,0.7};         // P1={ah,bh,ch}
double x_pair2[] ={0.003,0.005,0.007};   // P2={al,bl,cl}

double sum, sum_h, sum_l, sum_h2, sum_l2, prod, prod_h, prod_l, div_h, div_l; 
double xh6, xl6, xh, xl, ch, cl, ch3, cl3, c, c3, c4, t_h, t_l ;

long double summ, cc, prod_long, c_long, c_4;
long int c_44, N;

// SUMMATION ACCURACY TEST

    FILE *kk,*kkkk;
    kk=fopen("SUM2.txt","w");
    //usual double
    sum=0.0;
    c=3.0;

    //usual long double 
    summ=0.0;
    cc=3.0;

    // with algoritm
    sum_h=0.0;
    sum_l=0.0;
    ch =3.0;
    cl = 0.0;
    xh6=3.0;
    xl6=0.0;

    for (i=0; i<45; i++){
          //ORDINARY DOUBLE SUMMATION
    c=c/3.0;
    sum=sum+c;
         //ORDINARY LONG-DOUBLE SUMMATION
    cc=cc/3.0;
    summ=summ+cc;
         //SUMMATION WITH ARITHMETIC OPERATION
    Division(ch, cl, xh6, xl6, &ch, &cl);
    Division(ch, cl, xh6, xl6, &ch3, &cl3);
    Summation2(sum_h, sum_l, ch, cl, &sum_h, &sum_l);
    Summation2(sum_h, sum_l, ch3, cl3, &sum_h2, &sum_l2);
    fprintf(kk,"%d %.20Le %.20Le %.20Le %.20Le\n",i, summ, -(sum-summ)/summ, -(sum_h-summ+sum_l)/summ, -(sum_h2-summ+sum_l2)/summ);  
    } 
   
    FILE *kkk,*fff;
    kkk=fopen("MULT4.txt","w");
    //fff=fopen("MULT2.txt","w");

    prod=1.0;
    c3= 1/1000.0+1.0;

    prod_long=1.0;
    c_long= 1/1000.0+1.0;
  
    prod_h=1.0;
    prod_l=0.0;    
    t_h=1.0;
    t_l=0.0; 
    Division(1.0, 0.0, 1000.0, 0.0, &xh, &xl); // c3= 1/100000.0+1.0;
    printf("xh= %.16le xl=%.16le\n",xh,xl);
    Summation2(xh, xl, 1, 0, &ch, &cl);    
    printf("ch= %.16le cl=%.16le\n",ch,cl);
    printf("c3=%.16le\n",c3);
    printf("c_long=%.16Le\n",c_long); 

    for(i=0; i<100000; i++) { 
    // ORDINARY DOUBLE MULTIPLICATION
    prod= prod*c3;

    // ORDINARY LONG DOUBLE MULTIPLICATION
    prod_long = prod_long*c_long;

    // MULTIPLICATION ACCURACY TEST       
    Mult(prod_h, prod_l, c3, 0, &prod_h, &prod_l);
    Mult(t_h, t_l, ch, cl, &t_h, &t_l);
    fprintf(kkk,"%d %.22Le %.22Le %.22Le  \n",i, prod_long, -(prod-prod_long)/prod_long, -(prod_h-prod_long+prod_l)/prod_long);
    //fprintf(fff,"%d %.22Le %.22Le %.22Le  \n",i, prod_long, -(prod-prod_long)/prod_long, -(t_h-prod_long+t_l)/prod_long);
    } 
    //}
            printf("prod_h= %.16le prod_l=%.16le\n",prod_h,prod_l);
            printf("th= %.16le tl=%.16le\n",t_h,t_l); 

  // DIVISION ACCURACY TEST  

    FILE *tt;
    tt=fopen("DIV.txt","w"); 
    
    N= 100000; //number of division operation
    c4  = 1.0;    //double
    c_4 = 1.0;    //long double  
    
    div_h = 1.0; //pairs
    div_l = 0;

    for (i=0; i<=600; i++){
          //ORDINARY DOUBLE DIVISION
    c4 =c4/3.0;

         //ORDINARY LONG-DOUBLE DIVISION
    c_4 =c_4/3.0;

         //SUMMATION WITH ARITHMETIC OPERATION
    Division(div_h, div_l, 3.0, 0.0, &div_h, &div_l);

    fprintf(tt,"%d %.20Le %.20Le \n",i, -(c4-c_4)/c_4, -(div_h-c_4+div_l)/c_4 );  
    } 


    // Unit vector with pair / double / longdouble ACCURACY TEST
   unit_vector(x1, x2, x3, &I, &J, &K);
   unit_vector_longdouble(x1long, x2long, x3long, &T, &P, &L);  
   unit_vector_pairwise(xh_pair, xl_pair, xh_p, xl_p);

    printf("relative errors DOUBLE / LONG-DOUBLE\n");
    printf("  for X component= %.16Le \n",-((I-T)/T));
    printf("  for Y component= %.16Le \n",-((J-P)/P));
    printf("  for Z component= %.16Le \n",-((K-L)/L)); 

    printf("relative errors PAIRWISE / LONG-DOUBLE\n");
    printf("  for X component= %.16Le \n",-((xh_p[0]-T+xl_p[0])/T));
    printf("  for Y component= %.16Le \n",-((xh_p[1]-P+xl_p[1])/P));
    printf("  for Z component= %.16Le \n",-((xh_p[2]-L+xl_p[2])/L)); 
    
 fclose(kk);
 fclose(kkk); 
fclose(fff); 
 fclose(tt); 
return (0);  
}
