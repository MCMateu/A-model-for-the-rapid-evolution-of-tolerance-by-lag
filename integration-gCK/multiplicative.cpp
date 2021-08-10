#include<iostream>
#include<cmath>
#include<unistd.h>
#include<iomanip>
#include <stdio.h>
#include <ctime>

using namespace std;

int main (void)
{

	//Declaration of variables		
	int i,j;
	int ciclos,flag,contador;
	double b,d,s,eta,sigmam,sigmam2;
	double aux1,aux2_1,aux2_2,aux3,taumax,max,tmax,t1,t2,pi;
	double ht,hx,hx2,t,integral,var,skew;
	double scale,sigma,period,Ta,Tf;

	double Norm,mean,Nspace,Ntime,standar;
	double raiz2;
	double raizpi;
	raiz2=sqrt(2.0);	

	FILE* f1;
	FILE* f2;
	f1=fopen("multiplicative-CK-distribution.txt","w");	
	f2=fopen("multiplicative-CK-cumulants.txt","w");	
	//Set Initial Values of the Variables************************************************************************
			d=3.6E-5; //death rate (by natural causes)
			b=2.40;  //birth rate 
			s=0.12; //enter-in-dormant-state-rate
				
			period=23.0;//Duration of a complete cycle

			ciclos=15; //Total Number of exposure cycles

			sigma=1.27;


			tmax=ciclos*period; //final time
			taumax=30.0; //right limit of the phenotype space

			
			Ta=3.0; //duration of the killing phase (antibiotic)
			Tf=period; //auxiliar variable (marks the end of each cycle)

			Nspace=taumax*1E2; ///Total number of space-bins
			Ntime=tmax*1E4; //Total number of time-bins

			ht=(tmax-0)/(1.0*Ntime);  //size of the time-bin
			hx=(taumax-0)/(1.0*Nspace);//size of the space bin
			hx2=hx*hx; //auxiliar variable, makes the code a bit faster
			cout<<"ht : "<<ht<<"   hx: "<<hx<<"  Von Neumann "<<ht/(hx*hx)<<endl;
			pi=3.14159265359; //number PI
			raizpi=sqrt(pi);
			
			sigmam=0.048; //mutation width (alpha in the main text)
			sigmam2=sigmam*sigmam; //auxiliar variable, makes the code a bit faster


				//here I declare the vectors for the frecuencies 
				double pGi[int(Nspace)],pDi[int(Nspace)],pGf[int(Nspace)],pDf[int(Nspace)];

	//***********************************************************************************************	

	
	
	//BOUNDART CONDITIONS (Diritlech)
	pGi[0]=0.0;
	pDi[0]=0.0;
	pGf[0]=0.0;
	pDf[0]=0.0;
	pGi[int(Nspace)]=0.0;
	pDi[int(Nspace)]=0.0;
	pGf[int(Nspace)]=0.0;
	pDf[int(Nspace)]=0.0;
	
	
	
	//PVI: (Initial Value Problem)
	Norm=0.0;
	mean=0.0;
	//I am using a trunqued gaussian
	for(i=0;i<15;i++)
	{
		pGi[i]=0.0;
		pDi[i]=0.0;
			
		pGf[i]=0.0;
		pDf[i]=0.0;
	
	}

	for(i=15;i<(Nspace);i++)
	{
		pGi[i]=(1.0/sqrt(sigma*sigma*pi*0.5))*exp(-pow(i*hx-0.32,2)/(2.0*sigma*sigma));
		pDi[i]=0.0;
		Norm=Norm+hx*(pGi[i]+pDi[i]); 
		//Compute the norm
		
		pGf[i]=0.0;
		pDf[i]=0.0;
	}	

	//Normalization
	for(i=1;i<(Nspace);i++)
	{
		pGi[i]=pGi[i]/Norm;
		pDi[i]=pDi[i]/Norm;
		mean=mean+hx*hx*i*(pGi[i]+pDi[i]);
	}
			
	cout<<"norm:  "<<Norm<<"  mean:  "<<mean<<endl;
	
	
	//setting initial time
	t=0.0;

	//this are a pair of flag variables
	//I use them to save the data each "flag" iterations
	//Each cycle has 23*ciclos/ht instants of time but I don't want to save them all 
	contador=46000-1;
	flag=46000;
	

	while(t<tmax)
	{		
		

		//KILLING PHASE***********************************
		eta=-1.0; //this variable controls the medium state
		while(t<=Ta){
			
			//Integral
			integral=0.0;
			for(i=1;i<(Nspace);i++)
			{
				integral=integral+hx*pGi[i];
			}	
			
			Norm=0.0; //Set to zero variable Norm to use in a loop		
			for(i=1;i<(Nspace);i++)
			{
				
				//Frecuency GROWING STATE
				pGf[i]=ht*(-(b+d+s)*pGi[i]+1.0/(hx*i)*pDi[i])+pGi[i]-ht*(b*eta-d)*integral*pGi[i];
				
				//Frecuency DORMANT STATE
				pDf[i]=ht*(-1.0/(hx*i)*pDi[i]+s*pGi[i]-(b*eta-d)*integral*pDi[i])+pDi[i];
				Norm=Norm+hx*(pGf[i]+pDf[i]); //Compute the norm (only to check its conservation)
				
			}
			
			
			//Re-start the values of the frecuencies for the next cycle
			
			for(i=1;i<(Nspace);i++)
			{
					pGi[i]=pGf[i];
					pDi[i]=pDf[i];			
			}
				
			//Save the data in files, only in the last cycle
			//the rest of cycles are necessary to reach the periodic steady-state	
			if((ciclos-1)*0*23<t && t<(ciclos)*23)
			{	

				contador++;

				if(contador==flag)
				{	
					//set sumatories to zero
					mean=0.0;
					var=0.0;
					skew=0.0;

					//Compute the mean
					for(i=1;i<(Nspace);i++)
					{
						mean=mean+hx*hx*i*(pGi[i]+pDi[i]);						
					}
					
					//Compute the variance
					for(i=1;i<(Nspace);i++)
					{
						
						var=var+hx*(pGi[i]+pDi[i])*(hx*i-mean)*(hx*i-mean);
						
					}
					//Compute the Skewness
					standar=sqrt(var);
					for(i=1;i<(Nspace);i++)
					{
						skew=skew+hx*(pGi[i]+pDi[i])*((hx*i-mean)/standar)*((hx*i-mean)/standar)*((hx*i-mean)/standar);
					}

						i=1;
						//File with the evolution of the distribution during the last cycle
						while(i<(Nspace))
						{
								fprintf(f1,"%f %f \n",i*hx,(hx*(pGi[i]+pDi[i]))); 	
								
								i=i+1;
						}
						contador=0;
						//File with the evolution of the comulants during the last cycle
						fprintf(f2,"%f %f %f %f %f \n",t-(ciclos-1)*23,mean,var,skew,Norm);

				}
				
			}
			
		
			
			//Increase the time*****************************************
			t=t+ht;
			
		}
		
			cout<<t<<endl;
		
		//GROWING PHASE***********************************
		eta=1.0;
		while(t<=Tf){

		cout<<t<<endl;			
			//Integral
			integral=0.0;
			for(i=1;i<(Nspace);i++)
			{
				integral=integral+hx*pGi[i];
			}	
		
			Norm=0.0;
			for(i=1;i<(Nspace);i++)
			{
				
				//Frecuency GROWING STATE
 				aux1=ht*((b-d)*pGi[i]+1.0/(hx*i)*pDi[i]);
				aux2_1=(1.0/hx2)*(pGi[i+1]-2.0*pGi[i]+pGi[i-1])*(hx*i)*(hx*i)+2.0*pGi[i]+4.0*(hx*i)*(pGi[i+1]-pGi[i-1])/(2.0*hx);
				aux2_2=1.0-exp(-1.0/(2.0*sigmam2))/erfc(-1.0/(raiz2*sigmam))*(1.0/(raiz2*sigmam)-exp(-1.0/(2.0*sigmam2))/erfc(-1.0/(raiz2*sigmam)));
				aux3=(raiz2/raizpi)*sigmam*(exp(-(i+1)*hx/(2.0*sigmam2))/erfc(-1.0/(raiz2*sigmam))*pGi[i+1]-exp(-(i)*hx/(2.0*sigmam2))/erfc(-1.0/(raiz2*sigmam))*pGi[i])/hx;
			
				pGf[i]=aux1-2.0*(b*eta-d)*ht*aux3+ht*(b*eta-d)*sigmam2*aux2_1*aux2_2-ht*(b*eta-d)*integral*pGi[i]+pGi[i];

				//Frecuency DORMANT STATE
				pDf[i]=ht*(-(1.0/(i*hx))*pDi[i]-(b*eta-d)*integral*pDi[i])+pDi[i];
			
				Norm=Norm+hx*(pGi[i]+pDi[i]); //Norm
			}
			
			//Restart for next cycle
			for(i=1;i<(Nspace);i++)
			{
				pGi[i]=pGf[i];
				pDi[i]=pDf[i];
				
			}
			

			if((ciclos-1)*0*23<t && t<(ciclos)*23)
			{	

				contador++;

				if(contador==flag)
				{	
					//set sumatories to zero
					mean=0.0;
					var=0.0;
					skew=0.0;

					//Compute the mean
					for(i=1;i<(Nspace);i++)
					{
						mean=mean+hx*hx*i*(pGi[i]+pDi[i]);						
					}
					
					//Compute the variance
					for(i=1;i<(Nspace);i++)
					{
						
						var=var+hx*(pGi[i]+pDi[i])*(hx*i-mean)*(hx*i-mean);
						
					}
					//Compute the Skewness
					standar=sqrt(var);
					for(i=1;i<(Nspace);i++)
					{
						skew=skew+hx*(pGi[i]+pDi[i])*((hx*i-mean)/standar)*((hx*i-mean)/standar)*((hx*i-mean)/standar);
					}

					i=1;
						//File with the evolution of the distribution during the last cycle
						while(i<(Nspace))
						{
								fprintf(f1,"%f %f \n",i*hx,(hx*(pGi[i]+pDi[i]))); 	
								
								i=i+1;
						}
						contador=0;
						//File with the evolution of the comulants during the last cycle
						fprintf(f2,"%f %f %f %f %f \n",t-(ciclos-1)*23,mean,var,skew,Norm);

				}
				
			}
			
			
			//Increase the time***************************************
			t=t+ht;
			
		}
		
		Ta=Ta+period;
		Tf=Tf+period;

		
		
	
	}

	cout<<"media final  : "<<mean<<"  Norma final : "<<Norm<<endl;

//do for [i=1:1500] {plot 'schr.txt' index i u 1:2 w l}

}
