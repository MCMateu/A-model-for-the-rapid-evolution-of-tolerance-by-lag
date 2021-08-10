#include<iostream>
#include<cmath>
#include<unistd.h>
#include<iomanip>
#include <stdio.h>
#include <ctime>

using namespace std;

int main (void)
{

	//files
	FILE* f1;
	f1=fopen("figureS4-additive-alphaA-0.035.txt","w");
	f1=fopen("figureS4-mutliplicaitve-alphaM-0.048.txt","w");

	//variable declaration
	int i,N;	
	double mu_A,mu_M,sigma_A,sigma_M;
	double auxA,auxM;
	double h,tau_max;
	double alpha1,alpha2;
	double pi;

	pi=3.14159265359;

	tau_max=100.0;
	N=10000;
	h=(tau_max-0.0)/N;
	alpha1=0.035;
	alpha2=0.048;

	for(i=0;i<=N;i++)
	{
		auxA=1.0/erfc(-h*i/(sqrt(2.0)*alpha1));
		auxM=1.0/erfc(-1.0/(sqrt(2.0)*alpha2));
		mu_A=2.0*alpha1*exp(-(h*h*i*i)/(2.0*alpha1*alpha1))*(auxA);
		mu_M=sqrt(2.0/pi)*i*h*alpha2*exp(-1.0/(2.0*alpha2*alpha2))*auxM;
		sigma_A=alpha1*alpha1*(1.0-auxA*exp(-i*i*h*h/(2.0*alpha1*alpha1))*(i*h/(sqrt(2.0)*alpha1)-exp(-i*i*h*h/(2.0*alpha1*alpha1))*auxA));
		sigma_M=alpha2*alpha2*h*h*i*i*(1.0-auxM*exp(-1.0/(sqrt(2.0)*alpha2*alpha2))*(1.0/(sqrt(2.0)*alpha2)-auxM*exp(-1.0/(sqrt(2.0)*alpha2*alpha2))));
		fprintf(f1,"%f %f20 %f \n",h*i,mu_M,sigma_M);
		fprintf(f2,"%f %f %f \n",h*i,mu_A,sigma_A);
	}
	
	
}
