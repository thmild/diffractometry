/* p7fit.c */
/* Approximation of Thin film Data using Pearson VII Curves */

/* p7 evaluates Pearson VII - Curve at t  
   m = exponent, mu = location, a = scale parameter
   not normalized to give a density
*/

#include <math.h>

double p7 (double t, double m, double mu, double a) {
       double t2, t5 ,t9, t10;
       t2 = pow(t - mu, 0.2e1);
       t5 = a * a;
       t9 = pow(0.1e1 + t2 / m / t5, m);
       t10 = 0.1e1 / t9;
       return(t10);
       }

/* p7fit passt P7-Kurve an:
Reihenfolge des Parametervektors: 
            Int, Steig, masse_1,shape_1,lok_1,breite_1,... */


void p7fit (double *daten,
            double *fit,
            double *resid,
            double *param,
            double *param2,
            double *rss,
            double *versch,
            int *laenge,
            int *numker,
            int *gew,                                      /* 0 = Nein, 1 = Ja */
            double *weights                                /* Vektor der Gewichte */
           ){ int i,j ; 
           double temp1, temp2, temp3, temp4, temp5;
           double *parpointer;

           /* Zunächst Basislinie */

           temp1 = exp(*param++);
           temp2 = *param++;
           temp1 = *versch*(temp1-1)/(temp1+1);
           *param2++ = temp1;
           *param2++ = temp2;
                      
           for (i=0; i<=*laenge-1; i++) {
               fit[i]=temp1+(i+1)*temp2; 
               }
     
               /* Rücktransformation der Parameter */

               parpointer=param2;

               temp2=1;          /* Summanden */
               temp3=1;          /* Z */

           for (i=1; i<=*numker; i++) {

               *param2++=exp(*param++);              /* Masse */
               temp5=exp(*param++);
               *param2++=(1000*temp5+1)/(temp5+1);  /* Shape */
                        
               temp1=exp(*param++);      /* p_i*/      
               temp2=temp2*temp1;                         
               temp3=temp3+temp2;        /* Nenner */                  
               *param2++=temp1;          /* p_i*/        
                                  
               *param2++=exp(*param++);       /* Breite */
                   }
     
              /* Weitere Rücktrafo der Lokationen */
 
              temp2=parpointer[2];       /* p_i */
              temp1=(*laenge-1)/temp3;   /* h_i */
              temp4=1+temp1;             /* Speicher letzte Lok. */
              parpointer[2]=1+temp1;     /* 1. Lokation */ 
     
           for (i=2; i<=*numker; i++) {
              
              temp1=temp1*temp2;
              temp4=temp4+temp1;
              temp2=parpointer[2+(i-1)*4];
              parpointer[2+(i-1)*4]=temp4;
              }     

              /* Ende der Rücktransformation */

              /* Eigentliche Auswertung der P7-Kerne*/

           for (i=1; i<=*numker; i++) {
               temp1=parpointer[0+(i-1)*4];        /* Masse */
               temp2=parpointer[1+(i-1)*4];        /* Shape */
               temp3=parpointer[2+(i-1)*4];        /* Lok */
               temp4=parpointer[3+(i-1)*4];        /* Breite */
               for(j=1; j<=*laenge; j++) {
                        fit[j-1]=fit[j-1]+temp1*p7(j,temp2,temp3,temp4);
                        }
           }

           /* Berechnung der Residuen und RSS */
           *rss=0;
               for(j=1; j<=*laenge; j++) {
               temp1=fit[j-1]-daten[j-1];
               if(*gew==1) temp1=temp1*weights[j-1];   /* Falls Gewichtung erwünscht */
               resid[j-1]=temp1;
               *rss=*rss+temp1*temp1;
               }
       }
