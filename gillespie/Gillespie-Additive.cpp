#include<iostream>
#include<cmath>
#include<unistd.h>
#include<iomanip>
#include <gsl/gsl_rng.h>
#include<gsl/gsl_randist.h> //Para usar la distribución gaussiana
#include <stdio.h>
#include <ctime>
#include <vector> //compilar con -std=c++11 

gsl_rng * r;
using namespace std;

//Función externa para la búsqueda binaria

int binarySearch(const vector<double> &array, int L, int R, double rand)
{
	int m;
	while((R-L)!=1)
	{
		m=(L+R)/2;
		if(rand<array[m])
		{
			R=m;
		}
		else
		{
			L=m;
		}
	}

		m=(L+R)/2;
		return m;
}



int main (void)
{

				
	int i,j,K,t,l,contador,w,aux,flag,ciclos,p,c,A_c,m;
	long long int Gtot,Ztot,Ntot;
	double time,Tf,Ta,d,b,z,u1,u2,x,SPmax,mu,sigma,Ta_0;
	double suma,s,t1,t2,scale,sum,period,alpha,media,varianza;
	double time_interval, interval,A_time_interval;
	int ciclo_flag,n,seed;
	int w2;
	double step,max;
	int ciclo_step,F_aux_count,K_aux_count,indice;
	double aa;
	int time_frames;
					//Rates
				/*****/ z=1.0;						
				/*****/	d=3.6E-5; //death by natural causes rate
				/*****/	b=2.40; //birth rate
				/*****/	s=0.12; //enter from dormancy rate
					/*****/	alpha=0.035; //mutation parameter		
					/*****/
					K=100000; //maximum carrying capacity

				//Simulation parameters
				period=23.0; //Time of a whole cycle
				sigma=1.27; //standard deviation of the initial state
				mu=0.0;
				ciclos=200; //total number of cycles
			

				vector <double> tau(K+1),ZD(K+1),ZG(K+1),SP(K+2);




	//Start GSL function to generate random numbers
	const gsl_rng_type * T;
					
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	gsl_rng_set(r,84565667); //change the seed of the generator after calling it



	double mean,variance,skew;


	//FILES
	
				FILE* f1;
				FILE* f2;	
			
			
			
				f1=fopen("time-tau-P-3h-3930-mut0.035.txt","w");
				f2=fopen("Ta-D-G-N-3h.txt","w");
			

				interval=0.5; //time interval of the measures of distribution 
				
				step=0.01; //size of the bins of the histogram

				max=20.0; //size of the phenotype space (it starts in 0)
	
				seed=1; //Number of temporal series
					
				time_frames=int(23/interval)+2; //add +1 to avoid the thero index of c++  vectors (in the vector declaration), +1 to include the t=0
							     //and +1 to		

				//ARRAYS TO EL SUPER-HISTOGRAM
				vector < vector<int> > histograma( ((time_frames+1)), vector <int> (int(max/step)+1) );
				vector <long long int> total_particulas((time_frames+1)),total_particulas2((time_frames+1));

				//vector<tipo> nombre(tamaño); para delcarar vectores de dimension 1
				//vector< vector<tipo> >(tamaño, vector<tipo>(tamaño, 0)); para declarar vectores de dimensión 2

				vector <double> tiempo((time_frames+1)),tiempo2((time_frames+1));
			
				vector <vector<double>> P(((time_frames+1)),vector <double> (int(max/step)+1));

				vector<vector<long long int>> particulas(((time_frames+1)), vector <long long int> (int(max/step)+1));
				vector<vector<long long int>> particulas2(((time_frames+1)), vector <long long int> (int(max/step)+1));
				
				vector <int> bacteria(int(max/step)+1);

				double error_particulas,error_ptotal,error_frecuencia,error_tiempo;


				vector <double> mean_total_particulas((time_frames+1)),mean_total_particulas2((time_frames+1));
				vector <vector <double>> mean_particulas(((time_frames+1)), vector <double> (int(max/step)+1));
				vector <vector <double>> mean_particulas2(((time_frames+1)),vector <double> (int(max/step)+1));
				
	//START THE HISTOGRAM-**************************************************************************
	for(i=1;i<=(time_frames);i++)
	{
		for(j=1;j<=int(max/step);j++)
		{
			histograma[i][j]=0;
			particulas[i][j]=0;
			particulas2[i][j]=0;
		}
		tiempo[i]=0.0;
		total_particulas[i]=0;
		total_particulas2[i]=0;
		
	}


	int n_promedio;
	n_promedio=0;
	
	//TEMPORAL SERIES LOOP***************************************************
	for(l=1;l<=seed;l++)
	{


		//INICIO NUEVA SERIE TEMPORAL

				//INICIO VARIABLES
				
				period=23.0; //PERIOD OF A WHOLE CYCLE
				sigma=1.27;
				mu=0.0;
				
						

				Ta=3.0; //DURATION OF THE KILLING PHASE
				Ta_0=Ta;
				Tf=period; //END OF THE COMPLETE CYCLE
			
				Ztot=0;//TOTAL NUMBER OF CELLS IN "dormant state" 
				Gtot=10000; //TOTAL NUMBER OF CELLS IN "GROWING state" 
				Ntot=Gtot+Ztot; //TOTAL NUMBER OF CELLS
				
				p=1;		
				c=0;
				A_c=0;
				contador=1;	
				
				time=0.0;
				
				
				ciclo_step=5; //each many cycles it is measured
				ciclo_flag=50; //cycle from which measures starts

				n=seed*((ciclos-ciclo_flag)/ciclo_step+1); //total number of meassures of each instant of time

				
				time_interval=(ciclo_flag-1)*(period)+Ta_0;
				A_time_interval=(ciclo_flag-1)*(period);		


				
				//INIZIAL POPULATION*************************************************************************
				for(i=1;i<=Ntot;i++)
				{
						
					do{	
						x=gsl_ran_gaussian (r, sigma);
						tau[i]=x+mu;
						ZG[i]=1.0/tau[i];
					}while(tau[i]<=0);
				}

			


		//GILLESPIE-CORE*********************************************************************************************************************************************
		t=1;

		time=0.0;
		m=0;
		do
		{
			//KILLING PHASE
			while(time<=Ta && Ntot!=0)
			{
						
				//COMPUTE THE TOTAL RATES
				SP[0]=d*Gtot;
				SP[1]=SP[0]+b*Gtot;
				SP[2]=SP[1]+s*Gtot;
				
				for(j=2;j<Ztot+2;j++)
				{
					SP[j+1]=SP[j]+ZD[j-1];
				}
				
				SPmax=SP[Ztot+2];
						
				//SELECTION AND EXECUTION OF THE REACTION
				u1=SPmax*(gsl_rng_uniform (r));
				
					//DEATH BY NATURAL CAUSES
					if(u1<=SP[0])
					{
						w=gsl_rng_uniform_int (r,Gtot)+1;
						ZG[w]=ZG[Gtot];
						Gtot=Gtot-1;
						Ntot=Gtot+Ztot;
					}

					//DEATH WHEN TRYING TO REPRODUCE
					else if(SP[0]<u1 && u1<=SP[1])
					{
						 w=gsl_rng_uniform_int (r,Gtot)+1;
						 ZG[w]=ZG[Gtot];
						 Gtot=Gtot-1;
						 Ntot=Gtot+Ztot;		
					}
					//Enter in "dormant state"
					else if(SP[1]<u1 && u1<=SP[2])
					{		
					
						 w=gsl_rng_uniform_int (r,Gtot)+1;
						 Ztot=Ztot+1;
						 ZD[Ztot]=ZG[w];	
						 ZG[w]=ZG[Gtot];
						 Gtot=Gtot-1;		
					}
					//Exit from "dormant state"
					else if(SP[2]<u1 && u1<=SPmax && Ztot!=0)
					{
						aux=binarySearch(SP,2,Ztot+2,u1);
						Gtot=Gtot+1;
						ZG[Gtot]=ZD[aux-1];	
						ZD[aux-1]=ZD[Ztot];
						Ztot=Ztot-1;
					}
				
				//increase time
				u2=gsl_rng_uniform_pos (r);
				time=time-log(u2)/SPmax;			
					
			
	//****************Save the data of the KILLING PHASE if proceeds***************************************************
			if(time>=A_time_interval)
			{
							
		
				for(j=1;j<=int(max/step);j++)
				{
					bacteria[j]=0;
				}
				
				A_time_interval=A_time_interval+interval;
				m++;
				
				tiempo[m]=tiempo[m]+(time-(period)*(ciclo_flag-1));
				tiempo2[m]=tiempo2[m]+(time-(period)*(ciclo_flag-1))*(time-(period)*(ciclo_flag-1));

					aa=0;
					indice=1;
					while(aa<max)
					{
						for(j=1;j<=Ztot;j++)
						{		
							if((1.0*aa)<z/ZD[j] && z/ZD[j]<(1.0*aa+step))
							{
								histograma[m][indice]++;
								bacteria[indice]++;	
							}
						}

						for(j=1;j<=Gtot;j++)
						{		
							if(1.0*aa<z/ZG[j] && z/ZG[j]<(1.0*aa+step))
							{
								histograma[m][indice]++;	
								bacteria[indice]++;
							}
						}
							
						aa=aa+step;
							
							particulas[m][indice]=particulas[m][indice]+bacteria[indice];	
							particulas2[m][indice]=particulas2[m][indice]+bacteria[indice]*bacteria[indice];						

						indice++;	
					}

					total_particulas[m]=total_particulas[m]+Ntot;
					total_particulas2[m]=total_particulas2[m]+(Ntot*Ntot);

					n_promedio++;

					fprintf(f2,"%f %lli %lli %lli \n",(time-(period)*(ciclo_flag-1)),Ztot,Gtot,Ntot);
			
			}
		
		}//end of the killing phase
		


		//fresh medium phase
		while(time<=Tf && Ntot!=0)
		{
			
			//Compute the total rates
			SP[0]=d*Gtot;
			SP[1]=SP[0]+b*Gtot;

			for(j=1;j<Ztot+1;j++)
			{
				SP[j+1]=SP[j]+ZD[j];
			}
			SPmax=SP[Ztot+1];
			//Select and execute a microscopic raction
			u1=SPmax*(gsl_rng_uniform (r));
			
				//death by natural causes
				if(u1<=SP[0])
				{
					w=gsl_rng_uniform_int (r,Gtot)+1;
					ZG[w]=ZG[Gtot];
					Gtot=Gtot-1;
					Ntot=Gtot+Ztot;
				}
			
				//Reproduction
				else if(SP[0]<u1 && u1<=SP[1])
				{
					
					if(Ntot!=K)
					{					
					
						w=gsl_rng_uniform_int (r,Gtot)+1;
						Gtot=Gtot+1;
						
							ZG[Gtot]=1.0/(1.0/ZG[w]+gsl_ran_gaussian_tail(r, -(1.0/ZG[w]), alpha));
									
						
						

							ZG[w]=1.0/(1.0/ZG[w]+gsl_ran_gaussian_tail(r, -(1.0/ZG[w]), alpha));		
						
						Ntot=Gtot+Ztot;
					}
					else
					{
						w=gsl_rng_uniform_int (r,Gtot)+1;
						w2=gsl_rng_uniform_int (r,Gtot)+1;
						
							ZG[w2]=1.0/(1.0/ZG[w]+gsl_ran_gaussian_tail(r, -(1.0/ZG[w]), alpha));
						
						

							ZG[w]=1.0/(1.0/ZG[w]+gsl_ran_gaussian_tail(r, -(1.0/ZG[w]), alpha));		
					
						
					}
			
				}
				
				//exit from "dormant state"
				else if(SP[1]<u1 && u1<=SPmax && Ztot!=0)
				{
					aux=binarySearch(SP,1,Ztot+1,u1);
					Gtot=Gtot+1;
					ZG[Gtot]=ZD[aux];
					ZD[aux]=ZD[Ztot];
					Ztot=Ztot-1;		
				}

				//increase time
				u2=gsl_rng_uniform_pos (r);
				time=time-log(u2)/SPmax;
				
		
			
	//****************save FRESH MEDIUM data**************************************************
			
			if(time>=time_interval)
			{
				
				for(j=1;j<=int(max/step);j++)
				{
					bacteria[j]=0;
				}
				
				time_interval=time_interval+interval;
				m++;
			
				tiempo[m]=tiempo[m]+(time-(period)*(ciclo_flag-1));
				tiempo2[m]=tiempo2[m]+(time-(period)*(ciclo_flag-1))*(time-(period)*(ciclo_flag-1));

					aa=0;
					indice=1;
					while(aa<max)
					{
						for(j=1;j<=Ztot;j++)
						{		
							if((1.0*aa)<z/ZD[j] && z/ZD[j]<(1.0*aa+step))
							{
								histograma[m][indice]++;
								bacteria[indice]++;	
							}
						}

						for(j=1;j<=Gtot;j++)
						{		
							if(1.0*aa<z/ZG[j] && z/ZG[j]<(1.0*aa+step))
							{
								histograma[m][indice]++;	
								bacteria[indice]++;
							}
						}
							
						aa=aa+step;
							
							particulas[m][indice]=particulas[m][indice]+bacteria[indice];	
							particulas2[m][indice]=particulas2[m][indice]+bacteria[indice]*bacteria[indice];						

						indice++;	
					}

					total_particulas[m]=total_particulas[m]+Ntot;
					total_particulas2[m]=total_particulas2[m]+(Ntot*Ntot);

					n_promedio++;
					fprintf(f2,"%f %lli %lli %lli \n",(time-(period)*(ciclo_flag-1)),Ztot,Gtot,Ntot);
			
			}
		


		
				
		} //end fresh medium
		
		cout<<l<<"  "<<t<<endl;
						
		Ta=Ta+period;
		Tf=Tf+period;
	
		
		if(contador==ciclo_flag)
		{
			ciclo_flag=ciclo_flag+ciclo_step;
			time_interval=(ciclo_flag-1)*(period)+Ta_0;
			A_time_interval=(ciclo_flag-1)*(period);

			p++;

			m=0;

		}
		contador++;
		t++; 
		

	}while(t<=ciclos && Ntot!=0);
	//end GILLESPIE-CORE****************************************************************************************************************************************
	
}//end temporal series loop 


for(i=1;i<=(time_frames);i++)
{
	
	tiempo[i]=tiempo[i]/n;
	tiempo2[i]=tiempo2[i]/n;
	cout<<"tiempor  "<<i<<"\t"<<tiempo[i]<<endl;
}


cout<<n<<endl;




//Compute the super-histogram**************************************************************************
for(i=1;i<=(time_frames);i++)
{

		error_tiempo=sqrt((tiempo2[i]-tiempo[i]*tiempo[i]))/sqrt(n);

	
		mean_total_particulas[i]=total_particulas[i]/(1.0*n);
		mean_total_particulas2[i]=total_particulas2[i]/(1.0*n);
				
		error_ptotal=sqrt((mean_total_particulas2[i]-mean_total_particulas[i]*mean_total_particulas[i]))/sqrt(n);

		aa=step;
		for(j=1;j<=int(max/step);j++)
		{
			mean_particulas[i][j]=particulas[i][j]/(1.0*n);
			mean_particulas2[i][j]=particulas2[i][j]/(1.0*n);

			error_particulas=sqrt((mean_particulas2[i][j]-mean_particulas[i][j]*mean_particulas[i][j]))/sqrt(n);
		
			error_frecuencia=pow(error_particulas/(mean_total_particulas[i]),2)+pow(error_ptotal*mean_particulas[i][j]/(mean_total_particulas[i]*mean_total_particulas[i]),2);
			error_frecuencia=sqrt(error_frecuencia);	
				
			
			P[i][j]=(1.0*histograma[i][j])/(1.0*mean_total_particulas[i]*n);
			fprintf(f1,"%f %f %f %f %f \n",tiempo[i],aa,P[i][j],error_tiempo,error_frecuencia);
			aa=aa+step;
		}
}



}
