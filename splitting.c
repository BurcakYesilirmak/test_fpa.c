  #include <stdlib.h>
  #include <stdio.h>
  #include <math.h>
  #include <string.h>
/* void split(char *buffer, double *xh, double *xl){
     FILE *fp; 
     char ch[24], res[30], rest[10] , rest2[10], ress[20];
     char *ptr;
     double number1, number4;
     long int number2, number3, x1;
     int m, i, k, xx, digit_num_after_dec, zeros;
     fp = fopen("myfile.txt", "r");
     if(!fp)
    {
        printf("Could not open file\n");
    } 
        fgets(buffer, 50, fp);
        sscanf(buffer, "%s" , ch);
        sscanf(buffer, "%lf" , &number1); //original number taken from the buffer
        sprintf(res,"%.23lf", number1);  // assign number1 to a string

     xx = sizeof(ch) ;   //length of the array
     gcvt(number1,23,ress);	// convert number1 to a string array
     // compare original and converted strings and find the starting different point 
	printf("buffer converted is: %s\n",ress);
    printf("buffer2 original is: %s\n",ch);
     for(m=0; m<23; m++) {                 //find the similar finish point (m)
             if(ch[m]!=ress[m])
                           break;
     }                                // m=17 
          printf("point :\n%d\n", m); //different point with integer and comma 
   // convert different character parts to long integer
                for(int i=m, k=0; i<=xx && k<xx-m ;i++, k++) {  
     rest[k]= ch[i];
     rest2[k]= ress[i];}     
     // converts both strings to the long integer 
    sscanf(rest, "%ld" , &number2);
    sscanf(rest2, "%ld" , &number3);  
          printf("residual string for converted number:\n%ld\n", number3);  
          printf("residual string for original number:\n%ld\n", number2); 
          x1=(number3-number2);                   //find subtraction of long integers -- should it be labs?
          printf("difference:\n%ld\n", x1); 
          number4=x1*pow(10,-xx+2); // 2-> take into account integer and comma 
        printf("xl:\n%le\n", number4);
        printf("xh:\n%.25lf\n", number1); 
        *xh=number1;
        *xl=number4;
             fclose(fp);
}  */
  int main(int argc, char* argv[])
{
    char buffer[50];
FILE *fp; 
     char ch[24], res[30], rest[10] , rest2[10], ress[20];
     char *ptr;
     double number1, number4, number5;
     float number44, number11;
     long int number2, number3, x1;
     int m, i, k, xx, xxx;
     fp = fopen("myfile.txt", "r");
     if(!fp)
    {
        printf("Could not open file\n");
    } 
        fgets(buffer, 50, fp);
        sscanf(buffer, "%s" , ch);
        sscanf(buffer, "%lf" , &number1); //original number taken from the buffer
        sprintf(res,"%.23lf", number1);  // assign number1 to a string

     xx = sizeof(ch) ;   //length of the array
     gcvt(number1,24,ress);	// convert number1 to a string array
     // compare original and converted strings and find the starting different point 
	//printf("buffer converted is: %s\n",ress);
    //printf("buffer2 original is: %s\n",ch);
     for(m=0; m<24; m++) {                 //find the similar finish point (m)
             if(ch[m]!=ress[m])
                           break;
     }                                // m=17 
          printf("point :\n%d\n", m); 
   // convert different character parts to long integer
                for(int i=m, k=0; i<=xx && k<xx-m ;i++, k++) {  
     rest[k]= ch[i];
     rest2[k]= ress[i];}     
     // converts both strings to the long integer 
    sscanf(rest, "%ld" , &number2);
    sscanf(rest2, "%ld" , &number3);  
         printf("residual string for converted number:\n%ld\n", number3);  
         printf("residual string for original number:\n%ld\n", number2); 
          x1=(number2-number3); //find subtraction of long integers -- should it be labs?
          printf("difference:\n%ld\n",x1);                  
          m=m-2; //different point without integer and comma 
          number4=x1*pow(10,-xx+2); // 2-> take into account integer and comma (lowest part)
        number11=number1;
        number44=number4;
        number5=3.1000058789845856868798;
        printf("relatve error double %le\n", -(number1-number5+number4)/number5);
        printf("relative error float %le\n", -(number11-number5+number44)/number5);  
        printf("xh double:\n %.23lf\n",number1);
        printf("xl double:\n%.23lf\n",number4); 
        //printf("original:\n3.1000058789456845856868798 \n");  
        printf("original:\n 3.1000058789845856868798\n");  
        printf("xh float:\n %.22le\n", number11); 
        printf("xl float:\n%.22le\n", number44);
        printf("original:\n 3.1000058789456845856868798 \n");  
             fclose(fp);
     return (0);
}
   







