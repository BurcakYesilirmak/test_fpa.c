
#include <stdlib.h>
#include <stdio.h>
#include<math.h>
#include<float.h>
#include <time.h>
#include <limits.h>
#define C 134217729 
#define e 2.71828182845904523536028747  // --> 353602874713527

void VeltkampSplit(double x, double *xh, double *xl) {        //unsigned long C = B^(t/2) + 1;  //s=27 & B=2;
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
  double t,s,u,v,z;
  F2Sum(xh1, xh2, &t, &s);
  F2Sum(xl1, xl2, &u, &v);
  z = s+v+u;
  F2Sum(t, z,xh3,xl3);
  *xh3=t;
  *xl3=z;
return ;
}

//Pichat and Neumaier’s summation algorithm
void Summation2(double xh1, double xl1, double xh2, double xl2,double *xh3, double *xl3 ){  
  double t,s,u,v,z;
  TwoSum(xh1, xh2, &t, &s);
  TwoSum(xl1, xl2, &u, &v);
  z= s+v+u; 
  //F2Sum(t, z,xh3,xl3);
 *xh3=t;
  *xl3=z; 
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
double u1,u0;
*xh3 =xh1/xh2;
Mult(*xh3, 0, xh2, 0,&u0,&u1); 
*xl3 = xh1- u0;
*xl3 -= u1;
*xl3 += xl1;
*xl3 -= *xh3 * xl2;
*xl3 /= xh2;  
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


    double f(double R_newh, double unit)
    {
    return R_newh/unit-1 ;
    }

    void regula (double *x, double x0, double x1, double fx0, double fx1, int *itr)
    {
    *x = x0 - ((x1 - x0) / (fx1 - fx0))*fx0;
     ++(*itr);
    }

    void squareroot_unit(double X, double Y, double Z, double x1, double x2, double x3) {
    double i_new, j_new, k_new ;
    double Xnew_h, Xnew_l, Ynew_h, Ynew_l, Znew_h, Znew_l;
    double max_h, max_l, old_alp_h, old_alp_l, new_alp_h, new_alp_l, eps, max ;
    double X_reg, Y_reg, Z_reg, a , b;
    int i, itr =0, max_itr= 30;

    //finding maximum absolute value of the vector
    double k[]={ X, Y, Z }; 
    max=fabs(k[1]);
    for(i=0;i<3;i++) {
    if(fabs(k[i])>max) 
    max=fabs(k[i]);
    }
    printf("max=%f\n",max);
    Mult(max, 0, max, 0, &max_h, &max_l);
    Division(1, 0, max_h, max_l, &old_alp_h, &old_alp_l); // alpha = 1/r^2 's maximum absolute value
    printf("alpha_old_h = %6.4f\n", old_alp_h);
    printf("alpha_old_l = %6.4f\n", old_alp_l); 

    eps = 1e-13;
    a= 0.5;
    b= 1.5;

                      // multiply alpha with components of r vector 
    Mult(old_alp_h, old_alp_l, X, 0, &Xnew_h, &Xnew_l);
    Mult(old_alp_h, old_alp_l, Y, 0, &Ynew_h, &Ynew_l);
    Mult(old_alp_h, old_alp_l, Z, 0, &Znew_h, &Znew_l);
    printf("X_newh = %6.4f\n", Xnew_h);
    printf("Y_newh = %6.4f\n", Ynew_h);
    printf("Z_newh = %6.4f\n", Znew_h);

        // Starting Regular Falsi Method 
    regula (&x1, a, b, f(Xnew_h,a), f(Xnew_h,b), &itr);
    regula (&x2, a, b, f(Ynew_h,a), f(Ynew_h,b), &itr);
    regula (&x3, a, b, f(Znew_h,a), f(Znew_h,b), &itr);
    printf("x1 = %6.4f\n", x1);  
    printf("x2 = %6.4f\n", x2);  
    printf("x3 = %6.4f\n", x3);    

        do
         {     
            // FOR i 
        if (f(Xnew_h,a)*f(Xnew_h,x1) < 0)
            a= x1;
        else
            b= x1;
        regula(&i_new, a, b, f(Xnew_h,a), f(Xnew_h,b), &itr);
        if (fabs(i_new-x1) < eps)
        {
            printf("After %d iterations, root = %6.4f\n", itr, i_new);
        }
        x1=i_new;
             // FOR j 
       /* if (f(Ynew_h,a)*f(Ynew_h,x2) < 0)
            a= x2;
        else
            b= x2;
        regula(&j_new, a, b, f(Ynew_h,a), f(Ynew_h,b), &itr);
        if (fabs(j_new-x2) < eps)
        {
            printf("After %d iterations, root = %6.4f\n", itr, j_new);
        }
        x2=j_new; 
            // FOR k 
        if (f(Znew_h,a)*f(Znew_h,x3) < 0)
            a= x3;
        else
            b= x3;
        regula(&k_new, a, b, f(Znew_h,a), f(Znew_h,b), &itr);
        if (fabs(k_new-x3) < eps)
        {
            printf("After %d iterations, root = %6.4f\n", itr, k_new);
        }
        x3=k_new; */
          }
    while (itr<max_itr);
    printf("Solution does not converge or iterations not sufficient:\n");    
 } 

int main(int argc, char **argv)
{
   int i;
double xh1=0.37797,xl1=0.000293,xh2=0.335593,xl2=0.0001233, xh3=0.3355674, I, J, K;
double x_pair1[] ={0.3,0.5,0.7};         // P1={ah,bh,ch}
double x_pair2[] ={0.003,0.005,0.007};   // P2={al,bl,cl}
double c, c3, sum, sum_h, sum_l, chh, cll, sum_hh, sum_ll, prod, prod_h, prod_l, ch, cl;
double xh6, xl6, xh, xl;
double relative_errorr[i],relative_error_calc[i];
long double summ, cc, prod_long, c_long, z, t;

    FILE *kk;
    kk=fopen("SUM.txt","w");
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
    cl= 0.0;

    for (i=0; i<33; i++){
          //ORDINARY Double SUMMATION
    c=c/3.0;
    sum=sum+c;
         // ORDINARY Long-Double SUMMATION
    cc=cc/3.0;
    summ=summ+cc;
       //SUMMATION WITH ARITHMETIC OPERATION
    Division(ch, cl, xh6, xl6, &ch, &cl);
    Summation2(sum_h, sum_l, ch, cl, &sum_h, &sum_l);
    /* ch = chh;
    cl = cll;
    sum_h = sum_hh;
    sum_l = sum_ll; */
    // relative_errorr[i]= fabs(sum_h-1.5)/1.5;
    // log(fabs(sum-1.5)),log(fabsl(summ-1.5)),log(relative_errorr[i])
    fprintf(kk,"%d %.16Le %.16Le %.16Le\n",i, summ, sum-summ, sum_h+sum_l-summ);
    }

    FILE *kkk;
    kkk=fopen("MULT.txt","w");
    prod=1.0;
    c3= 1.0+1/100000.0;
    prod_long=1.0;
    c_long= 1.0+1/100000.0;
    prod_h=1.0;
    prod_l=0.0;   

    Division(1, 0, 100000.0, 0, &xh, &xl); // c3= 1.0+1/100000.0;
    Summation2(1, 0, xh, xl, &ch, &cl); 

    double *relative_error_ord = (double*)malloc(500);
    if (relative_error_ord != NULL) {  

    for(i=0; i<100000; i++) { 
    // ORDINARY DOUBLE MULTIPLICATION
    prod= prod*c3;

    // ORDINARY LONG DOUBLE MULTIPLICATION
    prod_long = prod_long*c_long;

    // MULTIPLICATION WITH ARITHMETIC OPERATION           
     Mult(prod_h, prod_l, ch, cl, &prod_h, &prod_l);
    // relative_error_ord[i]= -(prod_h-e+prod_l)/e;
    // log(fabs(prod-e)), log(fabsl(prod_long-e)), log(relative_error_ord[i])
    fprintf(kkk,"%d %.16Le %.16Le %.16Le  \n",i, prod_long, prod-prod_long, prod_h+prod_l-prod_long);
    }

    /* printf("ordinary summation = %.16le \n",sum); //  with geometric sqeuence or progression series 
    printf("relative error = %.16le \n\n",-(sum-1.5)/1.5); 
    printf("ordinary summation with Long double = %.16Le \n",summ); //  with geometric sqeuence or progression series 
    printf("relative error Long double = %.16Le \n\n",-(summ-1.5)/1.5); 
    printf(" sum_h = %.16le \n", sum_h); 
    printf(" sum_l = %le \n", sum_l);  
    printf(" relative error = %.16le \n", -(sum_h-1.5)/1.5);
    printf(" relative error2 = %.16le \n\n", -(sum_h-1.5+sum_l)/1.5);
    printf(" ordinary mult = %.16le \n",prod);
    printf(" relative error ordinary = %.16le \n\n",fabs(prod-e)/e);
    printf(" ordinary mult long = %.16Le \n",prod_long);
    printf(" relative error ordinary with long = %.16Le \n\n",fabsl(prod_long-e)/e);
    printf(" prod_h = %.16le \n",prod_h);
    printf(" prod_l = %.16le \n",prod_l);
    printf(" calculated error  = %.16le \n",fabs(prod_h-e)/e);
    printf(" calculated error2 = %.16le \n\n",fabs(prod_h-e+prod_l)/e);  */

    squareroot_unit(xh1, xh2, xh3,I, J, K);

fclose(kk);
fclose(kkk);
return (0); 
} 
}



