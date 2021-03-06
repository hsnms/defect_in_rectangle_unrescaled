#include <stdio.h>
#include <errno.h>
#include <math.h>
#define Up 1.570796326794897
#define Down1 1.570796326794897
#define Down2 1.570796326794897
#define Right 1.570796326794897
#define Left 1.570796326794897
#define Pi 3.141592653589793
#define S 0.00000000000001
#include "Header1new.h"
#include "Header2new.h"
/*#include "Header3new.h"*/
/*#include "Header4new.h"*/
#include "Header5new.h"
#include "Header6new.h"
#include "Header7new.h"
int main()
{
  int i,j,k,n,ncycle=5,sos,kk,radius,n1,radius2;
    double **u,*x,*y,*z,*z2,*q,**v,**f;
    double result,zero,rev2,rev3,v1,v2,v4,K,t,sum,min;
	double m,H,h,R0; int a,b; 
	FILE *fp,*fp2,*fp3;
	n=17;
	R0=6.0;
	H=1.0;
	K=0.74;
	printf("%d\n",(int)(n*R0-R0+1));
	fp = fopen ("Rectangle.txt","w");
	  if (fp == NULL) {
        printf ("File not created okay, errno = %d\n", errno);
        return 1;
    }
	  b=1;   //radius of the point or small radius of the ring
	 
	m=-1.0; 
    u=dmatrix2new(1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
	v=dmatrix2new(1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
	q=dvectornew(1,(int)(n*H-H+1));
	
	x=dvectornew(1,(int)(n*R0-R0+1));
	y=dvectornew(1,(int)(n*H-H+1));
	z=dvectornew(1,(int)(n*R0-R0+1));
	z2=dvectornew(1,(int)(n*H-H+1));
		for (i=1;i<=(int)(n*R0-R0+1);i++) //set the length
	{
        
		x[i]=R0*(i-1)*(1.0/((int)(n*R0-R0+1)-1));
	}
	for (i=1;i<=(int)(n*H-H+1);i++) //set the width
	{
        
		y[i]=H*(i-1)*(1.0/((int)(n*H-H+1)-1));
	}

	/*	for (j=1;j<=(int)(n*H-H+1);j++)//the side boundary
		{
			
			q[j]=R0;
		
		}
	
	for(i=1;i<=(int)(n*R0-R0+1);i++)//initial condition and boundary condition
		{
			    u[i][(int)(n*H-H+1)]=Up-0.5*m*Pi;
	    }
	for(j=2;j<=(int)(n*H-H+1)-1;j++)
		{
			for(i=2;i<=(int)(n*R0-R0+1)-1;i++)
			    u[i][j]=Up-0.5*m*Pi;
	    }

	for (j=1;j<=(int)(n*H-H+1);j++)
		{
			   
		  u[(int)(n*R0-R0+1)][j]=Right;   
		
			u[1][j]=Left-0.5*m*Pi;
			}*/


	if (m==-1.0) //hyperbolic type
	  { /*//start of the first major part
	
	    for (j=2;j<=1+b;j++) //inner boundary condition
	{
	  u[1+b][j]=-atan((y[j]-y[1])/((x[1+b]-x[1])))+1.5*Pi;
		
	}
	
	for (i=2;i<=1+b;i++)
	{
	  u[i][1+b]=-atan((y[1+b]-y[1])/((x[i]-x[1])))+1.5*Pi;
		
	}
	
		for (i=(int)(n*R0-R0+1)-b;i<=(int)(n*R0-R0+1)-1;i++)
	{
		u[i][(int)(n*H-H+1)-b]=Pi-atan((y[(int)(n*H-H+1)-b]-y[(int)(n*H-H+1)])/((x[i]-x[(int)(n*R0-R0+1)])));
		
	}
	for (j=(int)(n*H-H+1)-b;j<=(int)(n*H-H+1)-1;j++)
    {
		u[(int)(n*R0-R0+1)-b][j]=Pi-atan((y[j]-y[(int)(n*H-H+1)])/((x[(int)(n*R0-R0+1)-b]-x[(int)(n*R0-R0+1)])));
		
		}

	a=2*b+2;

	for (j=2;j<=1+b;j++)
	{
	  u[a+b][j]=atan((y[j]-y[1])/((x[a+b]-x[a])))+0.5*Pi;
	}
	for (i=a+1;i<=a+b;i++)
	{
	 u[i][1+b]=atan((y[1+b]-y[1])/((x[i]-x[a])))+0.5*Pi;
	}
	u[a][1+b]=Pi;
	for (i=a-b;i<=a-1;i++)
    {
     
      u[i][1+b]=-atan((y[1+b]-y[1])/((x[a]-x[i])))+1.5*Pi;
	}
 	for (j=2;j<=1+b;j++)
    {
      u[a-b][j]=-atan((y[j]-y[1])/((x[a]-x[a-b])))+1.5*Pi;
    }
    	for (i=1;i<=a;i++)
		u[i][1]=1.5*Pi;
	for(i=a+1;i<=(int)(n*R0-R0+1);i++)
			    u[i][1]=0.5*Pi;
    

     
	
	mgfas5new(u,n,ncycle,m,b,q,x,y, R0,K, H); //solve the PDE
		f=dmatrix2new(1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
		
		nrfunc4new(f,u,n,q,b,R0,K,H); //get the free energy density
	
n1=1;
 z[1]=integratnew(f, n1,(int)(n*R0-R0+1),n1,(int)(n*H-H+1),n, H, R0)-integratnew(f, n1,n1+b,n1,b+1,n, H, R0)-integratnew(f, (int)(n*R0-R0+1)-b,(int)(n*R0-R0+1),(int)(n*H-H+1)-b,(int)(n*H-H+1),n, H, R0); //get the free energy, but this one is slightly wrong
	zero=Stat3new(u,n,b,R0,K, H);
	result=z[1];
	for(i=1;i<=(int)(n*R0-R0+1);i++)
		for(j=1;j<=(int)(n*H-H+1);j++)
			v[i][j]=u[i][j];
	radius=1;
	
	printf("z[1]=%f\t\nzero=%f\n",z[1],zero);
	fprintf(fp,"z[1]=%f\t\n{0, %f},\nzero=%f\n",z[1],z[1],zero);*/


	

		kk=4+2*b;   //location of the smallest ring
a=kk;
	


	for(i=1;i<=(int)(n*R0-R0+1);i++)
		{
			    u[i][(int)(n*H-H+1)]=Up-0.5*m*Pi;
	    }
	for(i=2;i<=(int)(n*R0-R0+1)-1;i++)
		{
			for(k=2;k<=(int)(n*H-H+1)-1;k++)
			    u[i][k]=Up-0.5*m*Pi;
	    }


	for (j=1;j<=(int)(n*H-H+1);j++)
		{
			u[(int)(n*R0-R0+1)][j]=Right;
			u[1][j]=Left-0.5*m*Pi;
	    }


        for (j=2;j<=1+b;j++) //inner boundary condition
	{
	  u[1+b][j]=-atan((y[j]-y[1])/((x[1+b]-x[1])))+1.5*Pi;
		
	}

	for (i=2;i<=1+b;i++)
	{
	  u[i][1+b]=-atan((y[1+b]-y[1])/((x[i]-x[1])))+1.5*Pi;
		
	}

	for (i=(int)(n*R0-R0+1)-b;i<=(int)(n*R0-R0+1)-1;i++)
	{
		u[i][(int)(n*H-H+1)-b]=Pi-atan((y[(int)(n*H-H+1)-b]-y[(int)(n*H-H+1)])/((x[i]-x[(int)(n*R0-R0+1)])));
		
	}
	
	for (j=(int)(n*H-H+1)-b;j<=(int)(n*H-H+1)-1;j++)
        {
		u[(int)(n*R0-R0+1)-b][j]=Pi-atan((y[j]-y[(int)(n*H-H+1)])/((x[(int)(n*R0-R0+1)-b]-x[(int)(n*R0-R0+1)])));
		
	}

       
	for (j=2;j<=1+b;j++)
	{
         u[a+b][j]=atan((y[j]-y[1])/((x[a+b]-x[a])))+0.5*Pi;
	}
	for (i=a+1;i<=a+b;i++)
	{
           u[i][1+b]=atan((y[1+b]-y[1])/((x[i]-x[a])))+0.5*Pi;
	}
	u[a][1+b]=Pi;
	for (i=a-b;i<=a-1;i++)
       {
      u[i][1+b]=-atan((y[1+b]-y[1])/((x[a]-x[i])))+1.5*Pi;
	        
	}
	for (j=2;j<=1+b;j++)
        {
         u[a-b][j]=-atan((y[j]-y[1])/((x[a]-x[a-b])))+1.5*Pi;
	}
	for (i=1;i<=a;i++)
		u[i][1]=1.5*Pi;
	for(i=a+1;i<=(int)(n*R0-R0+1);i++)
	        u[i][1]=0.5*Pi;
	
	for (j=1;j<=(int)(n*H-H+1);j++)
		q[j]=R0;
       mgfas4new(u,n,ncycle,m,a,b,q,x,y, R0,K, H);
		f=dmatrix2new(1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
		
	nrfunc3new(f,u,n,q,kk,b,R0,K,H);
		
n1=1;
 z[kk]=integratnew(f, n1,(int)(n*R0-R0+1),n1,(int)(n*H-H+1),n, H, R0)-integratnew(f, a-b,a+b,n1,b+1,n, H, R0)-integratnew(f, (int)(n*R0-R0+1)-b,(int)(n*R0-R0+1),(int)(n*H-H+1)-b,(int)(n*H-H+1),n, H, R0)-integratnew(f, n1,1+b,n1,1+b,n, H, R0);//the areas of lattices on the boundary are cut halves
 

 //z[kk]=integratnew(f, n1,(int)(n*R0-R0+1),1+b,(int)(n*H-H+1)-b,n, H, R0)+integratnew(f, n1,(int)(n*R0-R0+1)-b,(int)(n*H-H+1)-b,(int)(n*H-H+1),n, H, R0)+integratnew(f, a+b,(int)(n*R0-R0+1),1,1+b,n, H, R0)+integratnew(f, 1+b,a-b,1,1+b,n, H, R0);
   
	zero=Stat2new(u,n,a,b,R0,K,H);


	result=z[kk];
	radius=kk;
	for(i=1;i<=(int)(n*R0-R0+1);i++)
		for(j=1;j<=(int)(n*H-H+1);j++)
			v[i][j]=u[i][j];
	        

	printf("z[%d]=%f\t\nzero=%f\n",kk,z[kk],zero);
	fprintf(fp,"z[%d]=%f\t\n{%f, %f},\nzero=%f\n",kk,z[kk],1.0*(kk-1)/(n-1),z[kk],zero);

	for(kk=4+2*b;kk<=(int)(n*R0-R0+1)-b-2;kk++)
	{
	 a=kk;
	

	u=dmatrix2new(1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
	q=dvectornew(1,(int)(n*H-H+1));
	
	
	/*for (i=1;i<=(int)(n*R0-R0+1);i++)
	{
        
		x[i]=R0*(i-1)*(1.0/((int)(n*R0-R0+1)-1));
	}
	for (i=1;i<=(int)(n*H-H+1);i++)
	{
        
		y[i]=H*(i-1)*(1.0/((int)(n*H-H+1)-1));
		}*/




        for(i=1;i<=(int)(n*R0-R0+1);i++)
		{
			    u[i][(int)(n*H-H+1)]=Up-0.5*m*Pi;
	    }
	for(i=2;i<=(int)(n*R0-R0+1)-1;i++)
		{
			for(k=2;k<=(int)(n*H-H+1)-1;k++)
			    u[i][k]=Up-0.5*m*Pi;
	    }

	for (j=1;j<=(int)(n*H-H+1);j++)
		{
			u[(int)(n*R0-R0+1)][j]=Right;
			u[1][j]=Left-0.5*m*Pi;
	    }

	for (j=2;j<=1+b;j++)
	{
         u[1+b][j]=-atan((y[j]-y[1])/((x[1+b]-x[1])))+1.5*Pi;
	}

	for (i=2;i<=1+b;i++)
	{
	  u[i][1+b]=-atan((y[1+b]-y[1])/((x[i]-x[1])))+1.5*Pi;
		
	}


	for (i=(int)(n*R0-R0+1)-b;i<=(int)(n*R0-R0+1)-1;i++)
	{
		u[i][(int)(n*H-H+1)-b]=Pi-atan((y[(int)(n*H-H+1)-b]-y[(int)(n*H-H+1)])/((x[i]-x[(int)(n*R0-R0+1)])));
		
	}

	for (j=(int)(n*H-H+1)-b;j<=(int)(n*H-H+1)-1;j++)
        {
		u[(int)(n*R0-R0+1)-b][j]=Pi-atan((y[j]-y[(int)(n*H-H+1)])/((x[(int)(n*R0-R0+1)-b]-x[(int)(n*R0-R0+1)])));
		
	}

	for (j=2;j<=1+b;j++)
	{
         u[a+b][j]=atan((y[j]-y[1])/((x[a+b]-x[a])))+0.5*Pi;
	}
	
	for (i=a+1;i<=a+b;i++)
	{
         u[i][1+b]=atan((y[1+b]-y[1])/((x[i]-x[a])))+0.5*Pi;
	}
	u[a][1+b]=Pi;
	for (i=a-b;i<=a-1;i++)
        {
       u[i][1+b]=-atan((y[1+b]-y[1])/((x[a]-x[i])))+1.5*Pi;
	}
	for (j=2;j<=1+b;j++)
        {
	 u[a-b][j]=-atan((y[j]-y[1])/((x[a]-x[a-b])))+1.5*Pi;
	}
	for (i=1;i<=a;i++)
		u[i][1]=1.5*Pi;
	for(i=a+1;i<=(int)(n*R0-R0+1);i++)
	        u[i][1]=0.5*Pi;

	for (j=1;j<=(int)(n*H-H+1);j++)
		q[j]=R0;
       mgfas4new(u,n,ncycle,m,a,b,q,x,y, R0,K, H);
		f=dmatrix2new(1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
		
nrfunc3new(f,u,n,q,kk,b,R0,K,H);
		
	n1=1;
z[kk]=integratnew(f, n1,(int)(n*R0-R0+1),n1,(int)(n*H-H+1),n, H, R0)-integratnew(f, a-b,a+b,n1,b+1,n, H, R0)-integratnew(f, (int)(n*R0-R0+1)-b,(int)(n*R0-R0+1),(int)(n*H-H+1)-b,(int)(n*H-H+1),n, H, R0)-integratnew(f, n1,1+b,n1,1+b,n, H, R0);//the areas of lattices on the boundary are cut halves

	
	zero=Stat2new(u,n,a,b,R0,K,H);

	if (z[kk]<result)
	{
	result=z[kk];
	radius=kk;
	for(i=1;i<=(int)(n*R0-R0+1);i++)
		for(j=1;j<=(int)(n*H-H+1);j++)
			v[i][j]=u[i][j];
	}
	printf("z[%d]=%f\t\nzero=%f\n",kk,z[kk],zero);
	
	fprintf(fp,"{%f, %f},",1.0*(kk-1)/(n-1),z[kk]);
	}


	} /*end of the first major*/



	/***/

	


			
	printf("result%f\t\nradius=%d\n",result,radius);	
	fprintf(fp,"m=%f\tR0=%f\nresult=%f\tradius=%d\n",m,R0,result,radius);

	min=z[2*b+5]-z[2*b+4]; radius2=2*b+5;
	printf("z[%d]-z[%d]=%f\n",2*b+5,2*b+4,z[2*b+5]-z[2*b+4]);
	fprintf(fp,"{%f, %f},", 1.0*(2*b+4)/(n-1),z[2*b+5]-z[2*b+4]);
for (kk=6+2*b;kk<=(int)(n*R0-R0+1)-b-2;kk++)
  {
    printf("z[%d]-z[%d]=%f\n",kk,kk-1,z[kk]-z[kk-1]);
    fprintf(fp,"{%f, %f},", 1.0*(kk-1)/(n-1),z[kk]-z[kk-1]);
    if(z[kk]-z[kk-1]<min)
      {min=z[kk]-z[kk-1];radius2=kk;}
  }
 printf("min=%f\tradius=%f\n",min,1.0*(radius2-1)/(n-1));
 fprintf(fp, "min=%f\tradius=%f\n",min,1.0*(radius2-1)/(n-1));
	
	fclose (fp);
                   

	
	fp2 = fopen ("RectangleVectorField.txt","w");

  
	if (fp2 == NULL) {
        printf ("File 2 not created okay, errno = %d\n", errno);
        return 1;
    }
	
		
		

		fprintf(fp2,"{");
		for(i=1;i<=(int)(n*R0-R0+1)-1;i++)
		{
			fprintf(fp2,"{");
			for(k=1;k<=(int)(n*H-H+1)-1;k++)
				fprintf(fp2,"{%f,\t%f},\t",sin(v[i][k]),cos(v[i][k]));
			fprintf(fp2,"{%f,\t%f}\t",sin(v[i][(int)(n*H-H+1)]),cos(v[i][(int)(n*H-H+1)]));
			fprintf(fp2,"},");
		}
			fprintf(fp2,"{");
			for(k=1;k<=(int)(n*H-H+1)-1;k++)
			  fprintf(fp2,"{%f,\t%f},\t",sin(v[(int)(n*R0-R0+1)][k]),cos(v[(int)(n*R0-R0+1)][k]));
			fprintf(fp2,"{%f,\t%f}\t",sin(v[(int)(n*R0-R0+1)][(int)(n*H-H+1)]),cos(v[(int)(n*R0-R0+1)][(int)(n*H-H+1)]));
			fprintf(fp2,"}");
		fprintf(fp2,"}");
		fclose (fp2);
		
    printf ("File created okay\n");


	
	free_dmatrix2new(u,1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
	free_dmatrix2new(v,1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
	free_dvectornew(q,1,(int)(n*H-H+1));
	
	free_dvectornew(x,1,(int)(n*R0-R0+1));
	free_dvectornew(y,1,(int)(n*H-H+1));
	free_dvectornew(z,1,(int)(n*R0-R0+1));
	free_dvectornew(z2,1,(int)(n*H-H+1));
    return 0;
}
