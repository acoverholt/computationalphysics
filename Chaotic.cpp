// Drew Overholt
// Hamilton's Equation Solver

#include <fstream>
#include <iostream>
#include <math.h>

void func(double t, double y[], double fun[], int n)

{
	double lambda, nu;


		nu=30.11;
		lambda=2.2/nu;


		fun[0]=1.0/nu-lambda/sqrt(y[1])*cos(t+sqrt(y[1])*sin(y[0]))*sin(y[0]); 	//equation for q1 dot
		fun[1]=2*lambda*sqrt(y[1])*cos(t+sqrt(y[1])*sin(y[0]))*cos(y[0]);				//equation for Eta dot
}

void rk4(double y[], double h, int Ntot, int n, double t)
{
	int i,I;
	double h2,Coef;

	h2=0.5*h;		//Number of Steps
	Coef=1.0/6.0;

	double k1[n], k2[n], k3[n], k4[n], work[n], fun[n]; 	//Array of y functions, Ym(t)

	for (I=1;(I<=Ntot);I++) {
		t=t+h*I;
		for (i=0;(i<=n-1);i++) {
			(*func)(t,y,fun, n);	//Step 1
			k1[i]=h*fun[i];
		//	y[i]=y[i]+k1[i];						//For Euler's Method
		}
		for (i=0;(i<=n-1);i++) {
			work[i]=y[i]+0.5*k1[i];
		}
		for (i=0;(i<=n-1);i++) {
			(*func)(t+h2,work,fun, n);	//Step 2
			k2[i]=h*fun[i];
		}
		for (i=0;(i<=n-1);i++) {
			work[i]=y[i]+0.5*k2[i];
		}
		for (i=0;(i<=n-1);i++) {
			(*func)(t+h2,work,fun, n);	//Step 3
			k3[i]=h*fun[i];
		}
		for (i=0;(i<=n-1);i++) {
			work[i]=y[i]+k3[i];
		}
		for (i=0;(i<=n-1);i++) {
			(*func)(t+h,work,fun, n);	//Step 4
			k4[i]=h*fun[i];
		}
		for (i=0;(i<=n-1);i++) {				//Combining into final y(t)
			y[i]=y[i]+Coef*(k1[i]+2.0*(k2[i]+k3[i])+k4[i]);
		}
	}

}

int main()

{
	double h;
	int Ntot;
	int n=2;		//Dimension of the problem

	double y[n], tau, pi=acos(-1.0), low, height, etai;
	int count;

	count=0;
	low=45;
	height=5;

	etai=low*low;

	  std::ofstream op;

	    op.open("output");

	    do {
	    	y[0]=0.0;     // phi
			y[1]=etai;  //Eta

	    	h=0.01;		//Time-step value

	    	if (etai <= (low+height/4)*(low+height/4)) {
	    		tau=pi*1.0/2.0;		//Starting time
	    	}
	    	if (etai <= (low+height*3/4)*(low+height*3/4) & etai > (low+height/4)*(low+height/4)) {
	    		tau=pi*1.0/2.0+count*pi/10;
	    		count=count+1;
	    	}
	    	if (etai > (low+height*3/4)*(low+height*3/4) ) {
	    		tau=pi*3.0/2.0;
	    	}
	    	Ntot=1;	//Total number of steps

	do {
		rk4(y, h, Ntot, n, tau); // Runs Rk4

		if ( sqrt((y[0]-(int(y[0]/(2*pi))+.5)*2*pi)*(y[0]-(int(y[0]/(2*pi))+.5)*2*pi)) <= 0.001) {
			if (sqrt(y[1]) <= low+height & sqrt(y[1]) >=low) {
	    op << tau-int(tau/(2*pi))*2*pi << "\t" << sqrt(y[1]) << "\n"; //Prints to file
		}
		}

        tau=tau+h; //Steps forward in time once

	} while (y[0] < pi*180);

	etai=(sqrt(etai)+height/20)*(sqrt(etai)+height/20);
	    } while (etai < (low+height)*(low+height));

    op.close();
    printf ("Done");

}
