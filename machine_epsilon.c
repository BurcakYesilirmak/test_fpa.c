// program can discover the avaible precision for the machine and produce machine epsilon
// it is executing on at execution time
// Burçak Yeşilırmak / 11 Sep 2021 / Akdeniz University

#include<stdio.h>
#include<float.h>
int main()
{
    int count=0;
   //long double eps=1.0;
    double eps=1.0;
    while(1.0+eps>1.0){
        eps=eps*0.5;
        count=count+1;
        printf(" %.40le %d\n",eps,count);
    }
    printf("     eps= %.60le\n", eps*2);
    printf("constant= %.60le\n", DBL_EPSILON);
    return 0;
}

