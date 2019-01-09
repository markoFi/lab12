#include <iostream>
#include <fstream>
#include <complex>

//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt, const int Nx, const double xmin);

void step(cmplx* const psi0, const double dx, const double dt, const int Nx,
          const double k, const double xmin);

void writeToFile(cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin, const double alpha,
                 const double lambda, const double omega, const double t);
//-----------------------------------
int main(){

	const int Nx = 300;
	const double xmin = -40;
  const double xmax = 40;
	const double Tend = 10*M_PI;
	const double dx = (xmax-xmin) / (Nx - 1);
	const double dt = dx / 100  ;
  double t = 0;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);

	const double lambda = 10;
  const double omega = 0.2;
  const double k = omega*omega;
  const double alpha = pow(k, 0.25);

  stringstream strm;

	cmplx* psi0 = new cmplx[Nx];

	init(psi0, alpha, lambda, dx, dt, Nx, xmin);

	writeToFile(psi0,"psi_0", dx,Nx,xmin, alpha, lambda, omega,t);

	for (int i = 1; i <= Na; i++) {
		for (int j = 1; j <= Nk-1; j++) {
          cout << j << endl;
         step(psi0, dx, dt, Nx, k, xmin);
         t+=dt;
		}
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin, alpha, lambda, omega,t);
	}
  cout << "t = " << t << endl;


  delete[] psi0;
	return 0;
}
//-----------------------------------
void step(cmplx* const psi0, const double dx, const double dt, const int Nx,
          const double k, const double xmin)
{
   // your code goes here
}
//-----------------------------------
void writeToFile(cmplx* const psi, const string s, const double dx,
                 const int Nx, const double xmin, const double alpha,
                 const double lambda, const double omega, const double t)
{
  ofstream out(s.c_str());
  double x, xi, xil;
  double h1, h2, h3;
  cmplx ana;
	for(int i=0; i<Nx; i++){
		x = xmin + i * dx;
    xi = alpha * x;
    xil = alpha * lambda;
    h1 = -0.5 * pow(xi - xil*cos(omega*t),2 );
    h2 = -omega*t/2 - xi * xil * sin(omega*t);
    h3 =  + 0.25 * xil*xil* sin(2*omega*t);
    ana = cmplx( h1 , h2 + h3  );
    ana = sqrt(alpha / sqrt(M_PI) ) * exp(ana);
		out << x << "\t" << norm(psi[i]) << "\t" << psi[i].real() << "\t" << psi[i].imag()
         << "\t" << norm(ana) << "\t" << ana.real() << "\t" << ana.imag() <<  endl;
	}
	out.close();
}
//-----------------------------------
void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt, const int Nx, const double xmin)
{
	const double x0 = dx*Nx * 0.5;
	for(int i=0;i<Nx; i++){
		double x = xmin + i*dx ;
		psi0[i] = sqrt(alpha/sqrt(M_PI)) * exp(- pow(alpha*(x-lambda),2)/2 );
	}
}
