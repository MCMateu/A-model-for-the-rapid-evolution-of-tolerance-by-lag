#include<iostream>
#include<cmath>
#include<unistd.h>
#include<iomanip>
#include <gsl/gsl_rng.h>
#include<gsl/gsl_randist.h> //Para usar la distribución gaussiana
#include <stdio.h>
#include <ctime>

gsl_rng * r;
using namespace std;

//Función externa para la búsqueda binaria

int binarySearch(double array[], int L, int R, double rand)
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
	//Declaración de variables

	int i,j,K,t,l,contador,w,aux,flag,ciclos;
	unsigned int Gtot,Ztot,Ntot;
	double time,Tf,Ta,d,b,z,u1,u2,x,SPmax,mu,sigma;
	double suma,s,t1,t2,sum,period,mutaciones;
	int w2;
	double step,max;
	double aa,N0;
		//Inicio variables 
			period=23.0; //Tiempo de un ciclo completo
			ciclos=8; //Número total de ciclos de exposición
			Ztot=0;     //número total de células en "dormant state"
			Gtot=100000; //número total de células en "growing state"			
			Ntot=Gtot+Ztot; //número total de células para un instante dado
			K=100000; //capacidad máxima de carga	
			sigma=1.27;
			mu=0.0;
			mutaciones=0.048;

			Tf=period;
			Ta=8.0;
			N0=Ntot;

				//Constantes que acompañan en los Rates
				d=3.6E-5;
				z=1.0;
				b=2.40;
				s=0.12;

	cout<<ciclos<<endl;		


	double tau[K+1],ZD[K+1],ZG[K+1],SP[K+2];
	int histograma[200+1];

	

	//Inicio funciones para la generación de número aleatorios
 	const gsl_rng_type * T;
	
	T = gsl_rng_default;
      	r = gsl_rng_alloc (T);
	gsl_rng_set(r,16977); //cambiar la semilla, pero despues de llamar al generador.


	//Estado inicial del sistema
	for(i=1;i<=Ntot;i++)
	{
		
		do{	
			x=gsl_ran_gaussian (r, sigma);
			tau[i]=x+mu;
			ZG[i]=z/tau[i];
		}while(tau[i]<=0);
		sum=sum+tau[i];

	}


	cout<<"t_lag medio inicial: "<<sum/Ntot<<endl;




	
	//Inicio variables para el bucle
	t=1;
	time=0.0;
	do
	{
		
		//CON Antibióticos
		while(time<=Ta && Ntot!=0)
		{
				
			//Cálculo de los rates
			SP[0]=d*Gtot;
			SP[1]=SP[0]+b*Gtot;
			SP[2]=SP[1]+s*Gtot;
			
			for(j=2;j<Ztot+2;j++)
			{
				SP[j+1]=SP[j]+ZD[j-1];
			}
			
			SPmax=SP[Ztot+2];
					
			//Selección y ejecución de la acción
			u1=SPmax*(gsl_rng_uniform (r));
			
				//Muerte por causas naturales
				if(u1<=SP[0])
				{
					w=gsl_rng_uniform_int (r,Gtot)+1;
					ZG[w]=ZG[Gtot];
					Gtot=Gtot-1;
					Ntot=Gtot+Ztot;
				}

				//Muerte al intentar reproducirse
				else if(SP[0]<u1 && u1<=SP[1])
				{
					//Si suponemos un rate constante, todas las células tienen igual probabilidad
					 w=gsl_rng_uniform_int (r,Gtot)+1;
					 ZG[w]=ZG[Gtot];
					 Gtot=Gtot-1;
					 Ntot=Gtot+Ztot;		
				}
				//Entrada en "dormant state"
				else if(SP[1]<u1 && u1<=SP[2])
				{		
					//Si suponemos un rate constante, todas las células tienen igual probabilidad
					 w=gsl_rng_uniform_int (r,Gtot)+1;
					 Ztot=Ztot+1;
					 ZD[Ztot]=ZG[w];	
					 ZG[w]=ZG[Gtot];
					 Gtot=Gtot-1;		
				}
				//Salida del "dormant state"
				else if(SP[2]<u1 && u1<=SPmax && Ztot!=0)
				{
					aux=binarySearch(SP,2,Ztot+2,u1);
					Gtot=Gtot+1;
					ZG[Gtot]=ZD[aux-1];	
					ZD[aux-1]=ZD[Ztot];
					Ztot=Ztot-1;
				}
			
			//Actualizo el tiempo
			u2=gsl_rng_uniform_pos (r);
			time=time-log(u2)/SPmax;


		}
			
		//SIN Antibióticos
		while(time<=Tf && Ntot!=0)
		{
			
			//Cálculo de los rates
			SP[0]=d*Gtot;
			SP[1]=SP[0]+b*Gtot;

			for(j=1;j<Ztot+1;j++)
			{
				SP[j+1]=SP[j]+ZD[j];
			}
			SPmax=SP[Ztot+1];
			//Selección y ejecución de la acción
			u1=SPmax*(gsl_rng_uniform (r));
			
				//Muerte por causas naturales
				if(u1<=SP[0])
				{
					w=gsl_rng_uniform_int (r,Gtot)+1;
					ZG[w]=ZG[Gtot];
					Gtot=Gtot-1;
					Ntot=Gtot+Ztot;
				}
				//Reproducción
				else if(SP[0]<u1 && u1<=SP[1])
				{
					
					if(Ntot!=K)
					{					
						//Si suponemos un rate constante, todas las células tienen igual probabilidad
						w=gsl_rng_uniform_int (r,Gtot)+1;
						Gtot=Gtot+1;
						ZG[Gtot]=1.0/(1.0/ZG[w]+gsl_ran_gaussian (r, (1.0/ZG[w])*mutaciones));
						ZG[w]=1.0/(1.0/ZG[w]+gsl_ran_gaussian (r, (1.0/ZG[w])*mutaciones));
						Ntot=Gtot+Ztot;
					}
					else
					{
						w=gsl_rng_uniform_int (r,Gtot)+1;
						w2=gsl_rng_uniform_int (r,Gtot)+1;
						ZG[w2]=1.0/(1.0/ZG[w]+gsl_ran_gaussian (r, (1.0/ZG[w])*mutaciones));
						ZG[w]=1.0/(1.0/ZG[w]+gsl_ran_gaussian (r, (1.0/ZG[w])*mutaciones));;
					}
				}
				
				//Salida del "dormant state"
				else if(SP[1]<u1 && u1<=SPmax && Ztot!=0)
				{
					aux=binarySearch(SP,1,Ztot+1,u1);
					Gtot=Gtot+1;
					ZG[Gtot]=ZD[aux];
					ZD[aux]=ZD[Ztot];
					Ztot=Ztot-1;		
				}

				//Actualizo el tiempo
				u2=gsl_rng_uniform_pos (r);
				time=time-log(u2)/SPmax;
				
				contador++;
				
		}
	

					
			
						
		Ta=Ta+period;
		Tf=Tf+period;
	
		t++; 
	
	}while(t<=ciclos && Ntot!=0);



	//***************************************************************
	//Para medir el t_lag medio final
	suma=0.0;
	
	for(j=1;j<=Ztot;j++)
	{		
		suma=suma+z/ZD[j];
	}

	for(j=1;j<=Gtot;j++)
	{		
		suma=suma+z/ZG[j];
	}
	//***************************************************************

	cout<<"t_lag medio final:   "<<suma/Ntot<<endl;


	//Inicio variables para el bucle
	contador=0;
	t=1;
	N0=Ntot;
	cout<<N0<<endl;
	time=0.0;
	flag=100;
	t1=clock();
	while(Ntot/N0>0.01)
	{
			
			//Cálculo de los rates
			SP[0]=d*Gtot;
			SP[1]=SP[0]+b*Gtot;
			SP[2]=SP[1]+s*Gtot;
			
			for(j=2;j<Ztot+2;j++)
			{
				SP[j+1]=SP[j]+ZD[j-1];
			}
			
			SPmax=SP[Ztot+2];
					
			//Selección y ejecución de la acción
			u1=SPmax*(gsl_rng_uniform (r));
			
				//Muerte por causas naturales
				if(u1<=SP[0])
				{
					w=gsl_rng_uniform_int (r,Gtot)+1;
					ZG[w]=ZG[Gtot];
					Gtot=Gtot-1;
					Ntot=Gtot+Ztot;
				}

				//Muerte al intentar reproducirse
				else if(SP[0]<u1 && u1<=SP[1])
				{
					//Si suponemos un rate constante, todas las células tienen igual probabilidad
					 w=gsl_rng_uniform_int (r,Gtot)+1;
					 ZG[w]=ZG[Gtot];
					 Gtot=Gtot-1;
					 Ntot=Gtot+Ztot;		
				}
				//Entrada en "dormant state"
				else if(SP[1]<u1 && u1<=SP[2])
				{		
					//Si suponemos un rate constante, todas las células tienen igual probabilidad
					 w=gsl_rng_uniform_int (r,Gtot)+1;
					 Ztot=Ztot+1;
					 ZD[Ztot]=ZG[w];	
					 ZG[w]=ZG[Gtot];
					 Gtot=Gtot-1;		
				}
				//Salida del "dormant state"
				else if(SP[2]<u1 && u1<=SPmax && Ztot!=0)
				{
					aux=binarySearch(SP,2,Ztot+2,u1);
					Gtot=Gtot+1;
					ZG[Gtot]=ZD[aux-1];	
					ZD[aux-1]=ZD[Ztot];
					Ztot=Ztot-1;
				}
			
			//Actualizo el tiempo
			u2=gsl_rng_uniform_pos (r);
			time=time-log(u2)/SPmax;
			
						
	
	}

	
	
	cout<<"MDK99"<<"\t"<<time<<"   "<<endl;


}
