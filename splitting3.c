  #include <stdlib.h>
  #include <stdio.h>
  #include <ctype.h>
  #include <math.h>
  #include <string.h>
  
  typedef struct{
  char str[22], res[30], res1[30], rest[10], rest1[10];
  double value, xl, xh;
  int xx, m, countt, counttt;
  long int num, num1, x1, num2, num3;
} numbers;

  void pairs(double *x_h, double *x_l) {
      numbers inp[7];
      char buffer[150];
      FILE *fp; 
      long int num1, num11, x1, num2, num22, x2, num3, num33, x3;
      double numa, numb;
      int m, i, k, l, o, p, t, xx, count, pp;
      long int D=8100000;

      fp = fopen("SPLIT3.txt", "r");
      if(!fp)
       {
        printf("Could not open file\n");
       } 
      // original numbers taken from the file
        count=0; 
      while (fgets(buffer,150, fp)!= NULL) { 
         //   printf("buffer=%s\n", buffer);
          for(pp=0; pp<150 ; pp++)  //looking for end of line 
              if(buffer[pp]=='\0')
                  break;
          sscanf(buffer, "%s", inp[count].str);
          sscanf(buffer, "%lf", &inp[count].value);
        //printf("pp=%d count=%d string=%s value=%.20lf \n", pp, count, inp[count].str, inp[count].value) ; 
          count++;

      for( p=0; p<pp-1; p++) {            
          if(isspace(buffer[p])){           //looking for space    
              sscanf(buffer+p+1, "%s", inp[count].str);
              sscanf(buffer+p+1, "%lf", &inp[count].value);
        //      printf("count=%d string=%s value=%lf \n", count, inp[count].str, inp[count].value) ; 
              count++;
              } }
        } 
        inp[3].value= inp[3].value/100;
        inp[4].value= inp[4].value/1000;
        inp[5].value= inp[5].value/1000;
     // inp[6].value= inp[6].value/100; 
      // printf("count = %d \n", count); 
      // exit(1);  

      for(i=0; i<count; i++) {
        inp[i].xx = sizeof(inp[i].str) ;         //size of number
        gcvt(inp[i].value,22,inp[i].res1);       // convert numbers to a string array     -ise 1 eksik
        //printf("size : %d \n\n", inp[i].xx); 
       }

       // compare original and converted strings and find the starting different point 
      for(i=0; i<count; i++) {
        for(inp[i].m=0; inp[i].m<inp[i].xx ; inp[i].m++) {                 //find the similar finish point (m)
          if(inp[i].str[inp[i].m]!= inp[i].res1[inp[i].m])
            break;     
                }    printf("point :\n%d\n", inp[i].m);  
      //   printf("string1= %s string2= %.22s x= %d value=%c  value2=%c\n", inp[i].str, inp[i].res1,  inp[i].xx, inp[i].str[inp[i].m], inp[i].res1[inp[i].m]);
      /* for(=0; i<6; i++) {
        inp[i].rest[k]= inp[i].str[t];
        inp[i].rest1[k]= inp[i].res1[t];  }   */
        } 

      for(i=0; i<count; i++) {
        for( t=inp[i].m, k=0; t<=inp[i].xx && k<inp[i].xx-inp[i].m ; t++, k++) {  
          inp[i].rest[k]= inp[i].str[t];
          inp[i].rest1[k]= inp[i].res1[t];  }   
   //  printf("string1= %s string2= %s \n\n", inp[i].rest, inp[i].rest1); 
      }
                                                                
     // converts both strings to the long integer 
      for(i=0; i<count; i++) {
        sscanf(inp[i].rest, "%ld" , &inp[i].num1);
        sscanf(inp[i].rest1, "%ld" , &inp[i].num2);  
     //    printf("string=%s value=%s \n", inp[i].rest, inp[i].rest1) ;  
     //    printf("num1=%ld num2=%ld \n\n", inp[i].num1, inp[i].num2) ; 
      }

      for(i=0; i<count; i++) {
         inp[i].num  = inp[i].num1 ;
         inp[i].num3  = inp[i].num2 ;
        inp[i].countt=0;
        inp[i].counttt=0;
        while (inp[i].num != 0) {
        inp[i].num /= 10; 
        ++inp[i].countt;
        }

        while (inp[i].num3 != 0) {
        inp[i].num3 /= 10; 
        ++inp[i].counttt; }         
       //   printf("Number of digits: 1 : %d  , 2: %d\n", inp[i].countt, inp[i].counttt);     
       
        if(inp[i].countt != inp[i].counttt){
            inp[i].num1 *= 10;  } 
          else ; //if(inp[i].countt == inp[i].counttt){
            inp[i].num1 *= 1;     
            inp[i].x1 = inp[i].num1-inp[i].num2;  
        //printf("num1=%ld num2=%ld num1-num2 = %ld \n\n", inp[i].num1, inp[i].num2,inp[i].x1) ;
                      }
                            
        //for(i=0; i<count; i++) {
           // --6.685074817026049E-01, 1.706558853671706E+00 -1.123043465687632E+00 -1.207950275095047E-02 -3.384553498022686E-03 3.350958555890135E-03
         // m=m-3; //different point without integer and comma 
         /* inp[0].x1 = inp[0].x1+D;
            printf("xX1 double:\n%lD\n", inp[0].x1); */
         inp[0].xl =-inp[0].x1*pow(10,-(inp[0].xx)+2); 
         inp[1].xl =inp[1].x1*pow(10,-(inp[1].xx)+2); 
         inp[2].xl =-inp[2].x1*pow(10,-(inp[2].xx)+3); 
         inp[3].xl =-inp[3].x1*pow(10,-(inp[3].xx)+1); 
         inp[4].xl =-inp[4].x1*pow(10,-(inp[4].xx)); 
         inp[5].xl = inp[5].x1*pow(10,-(inp[5].xx)-1); //} // 2-> take into account integer and comma (lowest part)
        // inp[6].xl = inp[6].x1*pow(10,-(inp[6].xx)); //} // 2-> take into account integer and comma (lowest part)
            x_h[0] = inp[0].value/10;
            x_h[1] = inp[1].value;
            x_h[2] = inp[2].value;
            x_h[3] = inp[3].value;
            x_h[4] = inp[4].value;
            x_h[5] = inp[5].value;
            //x_h[6] = inp[6].value; 
            for(i=0; i<count; i++) {
            x_l[i] = inp[i].xl;
            printf("xh double:\n%.21lf\n", x_h[i]);
            printf("xl double:\n%.21lf\n",x_l[i]);
            printf("xh+xl=\n%.21lf\n",x_h[i]+x_l[i]);
            printf(" %s\n\n",inp[i].str);  }
             fclose(fp);
      } 
   
    char* a[6]= {"-6.6850748170260490000", "1.7065588536717060000", "-1.1230434656876320000", 
                            "-1.2079502750950470000", "-3.3845534980226860000", "3.3509585558901350000"};
    int main(int argc, char* argv[]) {
      int i, count;
      double xh[6], xl[6];
      //printf("\n a[0] = %s", a[0]);
      //printf("\n a[1] = %s", a[1]);
      FILE *file_01;
      if((file_01   = fopen("SPLIT3.txt","w"))==NULL) {
        printf("CANNOT OPEN FILE OUT_a.txt 'w'\n");
       // exit(1);
       }
     //printf("\n strlen(a[0]) = %lu", strlen(a[0]));
     //printf("\n strlen(a[1]) = %lu", strlen(a[1]));

      fprintf(file_01, "%s %s %s %s %s %s\n", a[0], a[1], a[2], a[3], a[4], a[5]);
      fclose(file_01);

        pairs(xh, xl) ;
       //printf("\n\nFINISH!\n");
        return 0;
      }




















