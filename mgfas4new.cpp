#include <stdio.h>
#include <math.h>
#define NPRE 30 
#define NPOST /*1300*//**1000*//*1000*//*5000*/1700/*1400 for 33,65; 1700 for 129*/
#define ALPHA 0.33
#define NGMAX 15
#define Pi 3.141592653589793
#include "Header1new.h"
#include "Header2new.h"
void mgfas4new( double **u, int n, int maxcyc, double m, int a,int b, double *q,  double *x,  double *y, double R0, double K, double H)
{
	int j,i,jpost,lsw,jsw,ipass; 
	double h,res,v1,v2,v3,v4,v5,v6,v7,v8,v9,h2i,foh2,rjac,omega=1.0;

	h=1.0/(n-1); 
	h2i=1.0/(h*h);
	foh2=-4.0*h2i;
	
	rjac=/*0.9999*//*0.999*//*0.5*/0.999/*0.9 for 33,65; 0.99for 129*/;
	for (jpost=1;jpost<=NPOST;jpost++)
	{ 
             for(i=2;i<(int)(n*R0-R0+1);i++)
		{
			for (j=b+2;j<(int)(n*H-H+1)-b;j++)
			{
				        v1=u[i][j+1]-u[i][j];
						v2=u[i][j]-u[i][j-1];
						v4=u[i+1][j]-u[i][j];
						v5=u[i][j]-u[i-1][j];
						v6=u[i+1][j+1]-u[i-1][j+1];
						v7=u[i-1][j-1]-u[i+1][j-1];
		        			res=((v4-v5)+(v1-v2));
				
		        
						u[i][j] -= omega*res/(-4.0);
							if(u[i][j]>=1.5*Pi) u[i][j]=1.5*Pi;
	if(u[i][j]<=0.5*Pi) u[i][j]=0.5*Pi;
	
			
	         	}
		}


	     	for(i=b+2;i<a-b;i++)
		{
			for (j=2;j<=1+b;j++)
			{
			 v1=u[i][j+1]-u[i][j];
						v2=u[i][j]-u[i][j-1];
						v4=u[i+1][j]-u[i][j];
						v5=u[i][j]-u[i-1][j];
						v6=u[i+1][j+1]-u[i-1][j+1];
						v7=u[i-1][j-1]-u[i+1][j-1];
		        			res=((v4-v5)+(v1-v2));
				
		        
						u[i][j] -= omega*res/(-4.0);
							if(u[i][j]>=1.5*Pi) u[i][j]=1.5*Pi;
	if(u[i][j]<=0.5*Pi) u[i][j]=0.5*Pi;	      
			}
			}


		for(i=a+b+1;i<(int)(n*R0-R0+1);i++)
		{
			for (j=2;j<=1+b;j++)
			{
			 v1=u[i][j+1]-u[i][j];
						v2=u[i][j]-u[i][j-1];
						v4=u[i+1][j]-u[i][j];
						v5=u[i][j]-u[i-1][j];
						v6=u[i+1][j+1]-u[i-1][j+1];
						v7=u[i-1][j-1]-u[i+1][j-1];
		        			res=((v4-v5)+(v1-v2));
				
		        
						u[i][j] -= omega*res/(-4.0);
							if(u[i][j]>=1.5*Pi) u[i][j]=1.5*Pi;
	if(u[i][j]<=0.5*Pi) u[i][j]=0.5*Pi;	       
			}
		}


         for(i=2;i<(int)(n*R0-R0+1)-b;i++)
		{
			for (j=(int)(n*H-H+1)-b;j<(int)(n*H-H+1);j++)
			{
			 v1=u[i][j+1]-u[i][j];
						v2=u[i][j]-u[i][j-1];
						v4=u[i+1][j]-u[i][j];
						v5=u[i][j]-u[i-1][j];
						v6=u[i+1][j+1]-u[i-1][j+1];
						v7=u[i-1][j-1]-u[i+1][j-1];
		        			res=((v4-v5)+(v1-v2));
				
		        
						u[i][j] -= omega*res/(-4.0);
							if(u[i][j]>=1.5*Pi) u[i][j]=1.5*Pi;
	if(u[i][j]<=0.5*Pi) u[i][j]=0.5*Pi;	       
			}
		}
		omega=(jpost==1?1.0/(1.0-0.5*rjac*rjac):1.0/(1.0-0.25*rjac*rjac*omega));
	
	}

	
}
