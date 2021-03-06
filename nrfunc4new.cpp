#include <math.h>
#include <stdio.h>
#include "Header2new.h"
/*#define EPI 1.0e-7*/ 
#define Pi 3.141592653589793
void **nrfunc4new(double **f, double **u, int n, double *q,int b,double R0, double K, double H)
{
	double h,v1,v2,v4,v5,v6,v7,v8,v9;
	int i,j;
	h=1.0/(n-1); 	

	for(j=1;j<=(int)(n*H-H+1);j++)
		for(i=1;i<=(int)(n*R0-R0+1);i++)
			f[i][j]=0.0;

	

		for(i=2;i<(int)(n*R0-R0+1);i++) //part 1
	  {
	    for (j=b+2;j<(int)(n*H-H+1)-b;j++)
	    {
	    
			  
			v1=u[i][j+1]-u[i][j];
		    v2=u[i][j]-u[i][j-1];
			v4=u[i+1][j]-u[i][j];
			v5=u[i][j]-u[i-1][j];
		        
	        	f[i][j]=K*pow((v4+v5)/(2.0*h),2)+K*pow((v1+v2)/(2.0*h),2);
				
		}
	   }

	
			for(i=3*b+3;i<(int)(n*R0-R0+1);i++) //part 2
		{
			for (j=2;j<=1+b;j++)
			{
			v1=u[i][j+1]-u[i][j];
		    v2=u[i][j]-u[i][j-1];
			v4=u[i+1][j]-u[i][j];
			v5=u[i][j]-u[i-1][j];
			        
			
			f[i][j]=K*pow((v4+v5)/(2.0*h),2)+K*pow((v1+v2)/(2.0*h),2);
				
		        }
		}	

			// ignore part 2.0 (i=b+2;i<a-b;i++) (j=2;j<=1+b;j++)
	     


			
				for(i=2;i<(int)(n*R0-R0+1)-b;i++) //part 3
		{
			for (j=(int)(n*H-H+1)-b;j<(int)(n*H-H+1);j++)
			{
			v1=u[i][j+1]-u[i][j];
		    v2=u[i][j]-u[i][j-1];
			v4=u[i+1][j]-u[i][j];
			v5=u[i][j]-u[i-1][j];
			       
        	f[i][j]=K*pow((v4+v5)/(2.0*h),2)+K*pow((v1+v2)/(2.0*h),2);
				
		        }
		}

		
		



					for (i=3*b+3;i<=(int)(n*R0-R0+1)-1;i++) // part 4
	{
			v1=u[i][2]-u[i][1];
		    v2=u[i][2]-u[i][1];
			v4=u[i+1][1]-u[i][1];
			v5=u[i][1]-u[i-1][1];
				
			j=1;
		
        	f[i][j]=K*pow((v4+v5)/(2.0*h),2)+K*pow((v1+v2)/(2.0*h),2);
			

		}


	// ignore part 4.0 (i=a+b+1;i<=(int)(n*R0-R0+1)-1;i++) j=1
					

						for (i=2;i<=(int)(n*R0-R0+1)-b-1;i++) //part 5
	{
			v1=u[i][(int)(n*H-H+1)]-u[i][(int)(n*H-H+1)-1];
		    v2=u[i][(int)(n*H-H+1)]-u[i][(int)(n*H-H+1)-1];
			v4=u[i+1][(int)(n*H-H+1)]-u[i][(int)(n*H-H+1)];
			v5=u[i][(int)(n*H-H+1)]-u[i-1][(int)(n*H-H+1)];
		  
			j=(int)(n*H-H+1);
		
	        	f[i][j]=K*pow((v4+v5)/(2.0*h),2)+K*pow((v1+v2)/(2.0*h),2);
			

			}

							for (j=b+2;j<=(int)(n*H-H+1)-1;j++) //part 6
	{
		    v1=u[1][j+1]-u[1][j];
		    v2=u[1][j]-u[1][j-1];
			v4=u[2][j]-u[1][j];
			v5=u[2][j]-u[1][j];
	                i=1;
		
	        	f[i][j]=K*pow((v4+v5)/(2.0*h),2)+K*pow((v1+v2)/(2.0*h),2);
			
			}

								for (j=2;j<=(int)(n*H-H+1)-b-1;j++) //part 7
	{
		    v1=u[(int)(n*R0-R0+1)][j+1]-u[(int)(n*R0-R0+1)][j];
		    v2=u[(int)(n*R0-R0+1)][j]-u[(int)(n*R0-R0+1)][j-1];
			v4=u[(int)(n*R0-R0+1)][j]-u[(int)(n*R0-R0+1)-1][j];
			v5=u[(int)(n*R0-R0+1)][j]-u[(int)(n*R0-R0+1)-1][j];
	        		        
			i=(int)(n*R0-R0+1);
		
	        	f[i][j]=K*pow((v4+v5)/(2.0*h),2)+K*pow((v1+v2)/(2.0*h),2);
			
	}


		
								for (j=(int)(n*H-H+1)-b+1;j<=(int)(n*H-H+1)-1;j++) //part 8
	{
		    v1=u[(int)(n*R0-R0+1)-b][j+1]-u[(int)(n*R0-R0+1)-b][j];
		    v2=u[(int)(n*R0-R0+1)-b][j]-u[(int)(n*R0-R0+1)-b][j-1];
			v4=u[(int)(n*R0-R0+1)-b][j]-u[(int)(n*R0-R0+1)-b-1][j];
			v5=u[(int)(n*R0-R0+1)-b][j]-u[(int)(n*R0-R0+1)-b-1][j];
				        
			i=(int)(n*R0-R0+1)-b;
		
        	f[i][j]=K*pow((v4+v5)/(2.0*h),2)+K*pow((v1+v2)/(2.0*h),2);
		}

									for (j=2;j<=b;j++) //part 9
	{

	i=3*b+2;
	  
		    v1=u[i][j+1]-u[i][j];
		    v2=u[i][j]-u[i][j-1];
			v4=u[i+1][j]-u[i][j];
			v5=u[i+1][j]-u[i][j];
	        	        
		
	        	f[i][j]=K*pow((v4+v5)/(2.0*h),2)+K*pow((v1+v2)/(2.0*h),2);
			}


									for (i=2;i<=b;i++) //part 10
	{
		    v1=u[i][b+2]-u[i][b+1];
		    v2=u[i][b+2]-u[i][b+1];
			v4=u[i+1][b+1]-u[i][b+1];
			v5=u[i][b+1]-u[i-1][b+1];
			        
			j=b+1;
		
	        	f[i][j]=K*pow((v4+v5)/(2.0*h),2)+K*pow((v1+v2)/(2.0*h),2);
			}

										for (i=b+3;i<=3*b+1;i++) //part 10.2
	{
		    v1=u[i][b+2]-u[i][b+1];
		    v2=u[i][b+2]-u[i][b+1];
			v4=u[i+1][b+1]-u[i][b+1];
			v5=u[i][b+1]-u[i-1][b+1];
	       		        
			j=b+1;
		
	        	f[i][j]=K*pow((v4+v5)/(2.0*h),2)+K*pow((v1+v2)/(2.0*h),2);
	}


	// omit part 9.2 {i=1+b, (j=2; j<=b; j++)}
	// omit part 9.3 {i=2+b, (j=2; j<=b; j++)}

	

	

		for (i=(int)(n*R0-R0+1)-b+1;i<=(int)(n*R0-R0+1)-1;i++) //part 11
	{
		    v1=u[i][(int)(n*H-H+1)-b]-u[i][(int)(n*H-H+1)-b-1];
		    v2=u[i][(int)(n*H-H+1)-b]-u[i][(int)(n*H-H+1)-b-1];
			v4=u[i+1][(int)(n*H-H+1)-b]-u[i][(int)(n*H-H+1)-b];
			v5=u[i][(int)(n*H-H+1)-b]-u[i-1][(int)(n*H-H+1)-b];
			        
			j=(int)(n*H-H+1)-b;
		
	        	f[i][j]=K*pow((v4+v5)/(2.0*h),2)+K*pow((v1+v2)/(2.0*h),2);
			}



         
	  
			v1=u[1][(int)(n*H-H+1)]-u[1][(int)(n*H-H+1)-1];
		    v2=u[1][(int)(n*H-H+1)]-u[1][(int)(n*H-H+1)-1];
			v4=u[2][(int)(n*H-H+1)]-u[1][(int)(n*H-H+1)];
			v5=u[2][(int)(n*H-H+1)]-u[1][(int)(n*H-H+1)];
				       
			i=1;j=(int)(n*H-H+1); //part 12
		
			f[i][j]=K*pow((v4+v5)/(2.0*h),2)+K*pow((v1+v2)/(2.0*h),2);
			


			v1=u[(int)(n*R0-R0+1)][2]-u[(int)(n*R0-R0+1)][1];
		    v2=u[(int)(n*R0-R0+1)][2]-u[(int)(n*R0-R0+1)][1];
			v4=u[(int)(n*R0-R0+1)][1]-u[(int)(n*R0-R0+1)-1][1];
			v5=u[(int)(n*R0-R0+1)][1]-u[(int)(n*R0-R0+1)-1][1];
					        
			i=(int)(n*R0-R0+1);j=1; //part 13
		
			f[i][j]=K*pow((v4+v5)/(2.0*h),2)+K*pow((v1+v2)/(2.0*h),2);

	
            
	      
	
				i=3*b+2;j=1;
	v1=u[i][2]-u[i][1];
		    v2=u[i][2]-u[i][1];
			v4=u[i+1][1]-u[i][1];
			v5=u[i+1][1]-u[i][1];
			
		 //part 14
        	f[i][j]=K*pow((v4+v5)/(2.0*h),2)+K*pow((v1+v2)/(2.0*h),2);

		// omit part 14.2 i=1+b, j=1

        	// omit part 14.3 i=2+b, j=1
		
			
			v1=u[(int)(n*R0-R0+1)-b][(int)(n*H-H+1)]-u[(int)(n*R0-R0+1)-b][(int)(n*H-H+1)-1];
		    v2=u[(int)(n*R0-R0+1)-b][(int)(n*H-H+1)]-u[(int)(n*R0-R0+1)-b][(int)(n*H-H+1)-1];
			v4=u[(int)(n*R0-R0+1)-b][(int)(n*H-H+1)]-u[(int)(n*R0-R0+1)-b-1][(int)(n*H-H+1)];
			v5=u[(int)(n*R0-R0+1)-b][(int)(n*H-H+1)]-u[(int)(n*R0-R0+1)-b-1][(int)(n*H-H+1)]; 
			i=(int)(n*R0-R0+1)-b;j=(int)(n*H-H+1); //part 15
			f[i][j]=K*pow((v4+v5)/(2.0*h),2)+K*pow((v1+v2)/(2.0*h),2);
			
			
				v1=u[(int)(n*R0-R0+1)][(int)(n*H-H+1)-b]-u[(int)(n*R0-R0+1)][(int)(n*H-H+1)-b-1];
		    v2=u[(int)(n*R0-R0+1)][(int)(n*H-H+1)-b]-u[(int)(n*R0-R0+1)][(int)(n*H-H+1)-b-1];
			v4=u[(int)(n*R0-R0+1)][(int)(n*H-H+1)-b]-u[(int)(n*R0-R0+1)-1][(int)(n*H-H+1)-b];
			v5=u[(int)(n*R0-R0+1)][(int)(n*H-H+1)-b]-u[(int)(n*R0-R0+1)-1][(int)(n*H-H+1)-b];
			
			i=(int)(n*R0-R0+1);j=(int)(n*H-H+1)-b; //part 16
	        	f[i][j]=K*pow((v4+v5)/(2.0*h),2)+K*pow((v1+v2)/(2.0*h),2);

			
			
				i=1;j=b+1; //part 17
			v1=u[i][j+1]-u[i][j];
		    v2=u[i][j]-u[i][j-1];
			v4=u[i+1][j]-u[i][j];
			v5=u[i+1][j]-u[i][j];
	        	f[i][j]=K*pow((v4+v5)/(2.0*h),2)+K*pow((v1+v2)/(2.0*h),2);
			
				i=3*b+2;j=b+1; //part 18
			v1=u[i][j+1]-u[i][j];
		    v2=u[i][j]-u[i][j-1];
			v4=u[i+1][j]-u[i][j];
			v5=u[i][j]-u[i-1][j];
        	f[i][j]=K*pow((v4+v5)/(2.0*h),2)+K*pow((v1+v2)/(2.0*h),2);

		// omit part 18.2 i=1+b, j=1+b
		
		// omit part 18.3 i=2+b, j=1+b
		
			
			i=(int)(n*R0-R0+1)-b;j=(int)(n*H-H+1)-b; //part 19
			v1=u[i][j+1]-u[i][j];
		    v2=u[i][j]-u[i][j-1];
			v4=u[i+1][j]-u[i][j];
			v5=u[i][j]-u[i-1][j];
				  
			f[i][j]=K*pow((v4+v5)/(2.0*h),2)+K*pow((v1+v2)/(2.0*h),2);
	        
			
       return 0;
}
