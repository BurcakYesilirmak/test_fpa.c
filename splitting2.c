  #include <stdlib.h>
  #include <stdio.h>
  #include <ctype.h>
  #include <math.h>
  #include <string.h>
typedef struct{
  char str[23], res[30], res1[30], rest[10], rest1[10];
  double value, xl;
  int xx, m;
  long int num, num1, x1;
} numbers;

  int main(int argc, char* argv[])
{
    numbers inp[6];
    char buffer[200];
     FILE *fp; 
     long int num1, num11, x1, num2, num22, x2, num3, num33, x3;
     double numa, numb;
     int m, i, k, l, o, p, t, xx, count, countt;
     fp = fopen("denemeee.txt", "r");
     if(!fp)
    {
        printf("Could not open file\n");
    } 
    // original numbers taken from the file
        count=1;
        
        while (fgets(buffer,200, fp)!= NULL) { 
        //  printf("buffer=%s\n", buffer);
          sscanf(buffer, "%s", inp[0].str);
          sscanf(buffer, "%lf", &inp[0].value);
         // printf("string=%s value=%lf \n", inp[0].str, inp[0].value) ;  
          for(p=0;p<200; p++) {
            if(isspace(buffer[p])) { 
              sscanf(buffer+p+1, "%s", inp[count].str);
              sscanf(buffer+p+1, "%lf", &inp[count].value);
           // printf("string=%s value=%lf \n", inp[count].str, inp[count].value) ; 
              count++; }
           } }
          //printf("%d\n",count);
        for(i=0; i<count; i++) {
        inp[i].xx = strlen(inp[i].str) ;         //size
        if(inp[i].value > 0) {
        gcvt(inp[i].value,inp[i].xx,inp[i].res1); }
        else ;
         gcvt(inp[i].value,inp[i].xx+1,inp[i].res1);
         // convert numbers to a string array     -ise 1 eksik
       // printf("string1= %s string12= %s \n", inp[i].str, inp[i].res1);
        }

     // compare original and converted strings and find the starting different point 
      for(i=0; i<count; i++) {
      for(inp[i].m=0; inp[i].m<inp[i].xx ; inp[i].m++) {                 //find the similar finish point (m)
             if(inp[i].str[inp[i].m]!= inp[i].res1[inp[i].m])
                break;     
                }   // printf("point :\n%d\n", inp[i].m); 
      //printf("string1= %s string12= %.22s x= %d value=%c  value2=%c\n", inp[i].str, inp[i].res1,  inp[i].xx, inp[i].str[inp[i].m], inp[i].res1[inp[i].m]); 
       }
 

     for(i=0; i<count; i++) {
     for( t=inp[i].m, k=0; t<=inp[0].xx && k<inp[0].xx-inp[i].m ; t++, k++) {  
     inp[i].rest[k]= inp[i].str[t];
     inp[i].rest1[k]= inp[i].res1[t];}  
      printf("string1= %s string2= %s \n", inp[i].rest, inp[i].rest1); }


    
     // converts both strings to the long integer 
    for(i=0; i<count; i++) {
    sscanf(inp[i].rest, "%ld" , &inp[i].num);
    sscanf(inp[i].rest1, "%ld" , &inp[i].num1);    
    inp[i].x1 =(inp[i].num-inp[i].num1);          

    //printf("string=%s value=%s \n", inp[i].rest, inp[i].rest1) ;  
    printf("num1=%ld num2=%ld \n", inp[i].num, inp[i].num1) ; 
    printf("x1=%ld \n\n", inp[i].x1) ;
    }  

    countt=0;
        for(i=0; i<count; i++) {
         // m=m-3; //different point without integer and comma 
       // if (inp[i].x1<0 && inp[i].value <0) { 
         inp[i].xl =-inp[i].x1*pow(10,-(inp[0].xx)+3); } // 2-> take into account integer and comma (lowest part)
      /*  else if (inp[i].x1<0 && inp[i].value >0)
        { inp[i].xl = -inp[i].x1*pow(10,-(inp[0].xx)+2); }
        else if (inp[i].x1>0 && inp[i].value <0)
        { inp[i].xl = inp[i].x1*pow(10,-(inp[0].xx)+3); }
        else ; 
        inp[i].xl = inp[i].x1*pow(10,-(inp[0].xx)+2); } */
     countt=0;

        for(i=0; i<count; i++) {
        printf("xh double:\n%.21lf\n",inp[i].value);
        printf("xl double:\n%.21lf\n",inp[i].xl);
        printf("%s\n\n\n",inp[i].str);  }  

             fclose(fp);
     return (0);
}
   



// 2 and 5 are problematic



