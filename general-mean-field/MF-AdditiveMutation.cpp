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
	int i,j,k;
	int ciclos,flag,contador;
	double b,d,s,eta,sigmam,sigmam2;
	double aux1,aux2,taumax,max,tmax,t1,t2,PI;
	double ht,hx,hx2,t,integral1,integral2,var,skew;
	double scale,sigma,period,Ta,Tf;

	double Norm,mean,Nspace,Ntime,standar;
	double cte,variable,beta;
	double raiz2;	

	FILE* f1;
	FILE* f2;
	f1=fopen("additive-bigmut-distributions.txt","w");	
	f2=fopen("additive-bigmut-cumulants.txt","w");	
	//Set Initial Values of the Variables************************************************************************
			d=3.6E-5; //death rate (by natural causes)
			b=2.40;   //birth rate 
			s=0.12;   //enter-in-dormant-state-rate
				
			period=23.0; //Duration of a complete cycle

			ciclos=10; //Total Number of exposure cycles
			sigma=1.27;


			tmax=ciclos*period; //final time 
			taumax=20.0; //right limit of the phenotype space

			
			Ta=3.0; //duration of the killing phase (antibiotic)
			Tf=period; //auxiliar variable (marks the end of each cycle)

			
			Nspace=taumax*0.666E1; //Total number of space-bins
			Ntime=tmax*1E3; //Total number of time-bins

			ht=(tmax-0)/(1.0*Ntime); //size of the time-bin
			hx=(taumax-0)/(1.0*Nspace); //size of the space bin
			hx2=hx*hx; //auxiliar variable, makes the code a bit faster
			cout<<"ht : "<<ht<<"   hx: "<<hx<<"  Von Neumann "<<ht/(hx*hx)<<endl;
			PI=3.14159265359; //number PI

			sigmam=0.16; //mutation width (alpha in the main text)
			sigmam2=sigmam*sigmam;//auxiliar variable, makes the code a bit faster

			cte=1.0/sqrt(0.5*PI*sigmam2); //normalization constant of the mutation kernel
			raiz2=sqrt(2.0);
				//here I declare the vectors for the frecuencies
				double pGi[int(Nspace)],pi[int(Nspace)],pGf[int(Nspace)],pf[int(Nspace)];

	//***********************************************************************************************	

	
	
	//BOUNDART CONDITIONS (Diritlech)
	pGi[0]=0.0;
	pi[0]=0.0;
	pGf[0]=0.0;
	pf[0]=0.0;
	pGi[int(Nspace)]=0.0;
	pi[int(Nspace)]=0.0;
	pGf[int(Nspace)]=0.0;
	pf[int(Nspace)]=0.0;
	
	
	
	//PVI: (Initial Value Problem)
	Norm=0.0;
	mean=0.0;
	//I am using a trunqued gaussian
	for(i=1;i<(Nspace);i++)
	{
		pGi[i]=(1.0/sqrt(sigma*sigma*PI*0.5))*exp(-pow(i*hx,2)/(2.0*sigma*sigma));
		pi[i]=pGi[i];
		Norm=Norm+hx*(pi[i]); //Compute the norm
		
		pGf[i]=0.0;
		pf[i]=0.0;
	}	

	//Normalization
	for(i=1;i<(Nspace);i++)
	{
		pGi[i]=pGi[i]/Norm;
		pi[i]=pi[i]/Norm;
		mean=mean+hx*hx*i*(pi[i]);
	}
			
	
	

	cout<<"norm:  "<<Norm<<"  mean:  "<<mean<<endl;
	
	
	//setting initial time
	t=0.0;

	//this are a pair of flag variables
	//I use them to save the data each "flag" iterations
	//Each cycle has 23*ciclos/ht instants of time but I don't want to save them all 
	contador=500-1;
	flag=500;


	//Main loop
	while(t<tmax)
	{
		
		cout<<"Start cycle: "<<t<<endl;
		//KILLING PHASE***********************************
		eta=-1.0; //this variable controls the medium state
		while(t<=Ta){
			
			//Integral1
			integral1=0.0;
			for(k=1;k<(Nspace);k++)
			{
				integral1=integral1+hx*pGi[k];
				
			}	
			
			Norm=0.0; //Set to zero variable Norm to use in a loop		
			for(i=1;i<(Nspace);i++)
			{
				//Total frecuency (GROWING+DORMANT)
				pf[i]=ht*(-(b+d)*pGi[i]+(b+d)*pi[i]*integral1)+pi[i];
				//Frecuency GROWING STATE
				pGf[i]=ht*(-(b+d)*pGi[i]+(b+d)*pGi[i]*integral1+1.0/(i*hx)*(pi[i]-pGi[i])-s*pGi[i])+pGi[i];
				
				Norm=Norm+hx*(pf[i]); //Compute the norm (only to check its conservation)
				
			}
			
			
			//Re-start the values of the frecuencies for the next cycle
			for(i=1;i<(Nspace);i++)
			{
					pGi[i]=pGf[i];
					pi[i]=pf[i];			
			}

			//Save the data in files, only in the last cycle
			//the rest of cycles are necessary to reach the periodic steady-state	
			/*if((ciclos-1)*23<t && t<(ciclos)*23)
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
						mean=mean+hx*hx*i*(pi[i]);						
					}
					
					//Compute the variance
					for(i=1;i<(Nspace);i++)
					{
						
						var=var+hx*(pi[i])*(hx*i-mean)*(hx*i-mean);
						
					}
					//Compute the Skewness
					standar=sqrt(var);
					for(i=1;i<(Nspace);i++)
					{
						skew=skew+hx*(pi[i])*((hx*i-mean)/standar)*((hx*i-mean)/standar)*((hx*i-mean)/standar);
					}

						i=1;
						//File with the evolution of the distribution during the last cycle
						while(i<(Nspace))
						{
								fprintf(f1,"%f %f \n",i*hx,(hx*pi[i])); 	
								
								i=i+1;
						}
						contador=0;
						//File with the evolution of the comulants during the last cycle
						fprintf(f2,"%f %f %f %f %f \n",t-(ciclos-1)*23,mean,var,skew,Norm);

				}
				
			}*/
			
			//Increase the time*****************************************
			t=t+ht;
			
		}
		
		cout<<"Middle of the cycle: "<<t<<endl;
		
		//GROWING PHASE***********************************
		eta=1.0;
		while(t<=Tf){
			
			
			//Integral1
			integral1=0.0;
			for(i=1;i<(Nspace);i++)
			{
				integral1=integral1+hx*pGi[i];
			}	
			
			Norm=0.0;
			for(i=1;i<(Nspace);i++)
			{
				integral2=0.0;
				//Integral2
				for(k=1;k<(Nspace);k++)
				{
					variable=hx*hx*(i-k)*(i-k)/(2.0*sigmam2);
					beta=cte/(erfc(-hx*k/(raiz2*sigmam)))*exp(-variable);
					integral2=integral2+hx*pGi[k]*beta;	
				}
				
			//Total Frecuency
			pf[i]=ht*(2.0*b*integral2-(b-d)*pi[i]*integral1-(b+d)*pGi[i])+pi[i];
					
			//Growing State
			 pGf[i]=ht*(2.0*b*integral2+1.0/(i*hx)*(pi[i]-pGi[i])-(b-d)*pGi[i]*(integral1)-(b+d)*pGi[i])+pGi[i];				
		
				
				Norm=Norm+hx*(pf[i]); //Compute the norm
			}
			
			//Re-start the values of the frecuencies for the next cycle
			for(i=1;i<(Nspace);i++)
			{
				pGi[i]=pGf[i];
				pi[i]=pf[i];
			}
			
			

			/*if((ciclos-1)*23<t && t<(ciclos)*23)
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
						mean=mean+hx*hx*i*(pi[i]);						
					}
					
					//Compute the variance
					for(i=1;i<(Nspace);i++)
					{
						
						var=var+hx*(pi[i])*(hx*i-mean)*(hx*i-mean);
						
					}
					//Compute the Skewness
					standar=sqrt(var);
					for(i=1;i<(Nspace);i++)
					{
						skew=skew+hx*(pi[i])*((hx*i-mean)/standar)*((hx*i-mean)/standar)*((hx*i-mean)/standar);
					}

					i=1;
						//File with the evolution of the distribution during the last cycle
						while(i<(Nspace))
						{
								fprintf(f1,"%f %f \n",i*hx,(hx*pi[i])); 	
								
								i=i+1;
						}
						contador=0;
						//File with the evolution of the comulants during the last cycle
						fprintf(f2,"%f %f %f %f %f \n",t-(ciclos-1)*23,mean,var,skew,Norm);

				}
				
			}*/
			
			//Increase the time***************************************
			t=t+ht;
			
		}
		
		//I add a whole period-duration for the next cycle		
		Ta=Ta+period;
		Tf=Tf+period;

		
		
	
	}
		
	//I compute the final norm to check that the Norm has been conserved as expected
	cout<<"Final Mean  : "<<mean<<"  Final Norm : "<<Norm<<endl;


	
						i=1;
						//File with the evolution of the distribution during the last cycle
						while(i<(Nspace))
						{
								fprintf(f1,"%f %f \n",i*hx,(hx*pi[i])); 	
								
								i=i+1;
						}


}
