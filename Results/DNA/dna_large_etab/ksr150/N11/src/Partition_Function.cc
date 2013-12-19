#include <cmath>
#include <string>
#include <sstream>
#include <complex>
#include <vector>
#include "Potential.hh"
#include "Matrix.hh"
#include "Vector.hh"
#include "Eigenfunctions.hh"
#include "Eigenvalues.hh"

#define PI 3.1415926535897932384626433832795

using namespace std;

extern double L;
extern int m;
extern double Delta;
extern int N;

extern int smax;
extern int tmax;

extern Vector t00_eval;
extern Matrix t00_evec;

extern vector<vector<complex<double> > > sVt00;
extern vector<vector<complex<double> > > CsVt00;
extern vector<vector<complex<double> > > s1_xi_s2;
extern vector<vector<double> > t_eta_t2;


complex<double> EXP_sVt00(int s, int t)
{
	complex<double> _EXP_sVt00 = (0,0);
	for(int i=-m;i<=m;i++){
		_EXP_sVt00 = _EXP_sVt00 + Delta*PSI_S(s,Delta*i)*exp(-Potential(Delta*i)/2)*t00_evec.get(t,i+m);
	}
	return _EXP_sVt00;
}

complex<double> EXP_CsVt00(int s, int t)
{
	complex<double> _EXP_CsVt00 = (0,0);
	for(int i=-m;i<=m;i++){
		_EXP_CsVt00 = _EXP_CsVt00 + Delta*conj(PSI_S(s,Delta*i))*exp(-Potential(Delta*i)/2)*t00_evec.get(t,i+m);
	}
	return _EXP_CsVt00;
}

vector<double> Partition_Function(int Extension){
	vector<double> _PF_DATA (2);
	double Sum_PF = 0;
	double Sum_DPF = 0;
    #pragma omp parallel for reduction(+:Sum_PF,Sum_DPF)
    for(int s=-smax;s<=smax;s++){
        complex<double> s_complex(0.0,(-(4.*PI*s)/L));
        for(int t=0;t<tmax;t++){
            Sum_PF = Sum_PF + pow(lambda_s(s),(N-1))
            *pow(t00_eval.get(t),(N-1))
            *real(
                  g(s,Extension*Delta)
                  *sVt00[s+smax][t]
                  *sVt00[s+smax][t]
                  );
            Sum_DPF = Sum_DPF + pow(lambda_s(s),(N-1))
            *pow(t00_eval.get(t),(N-1))
            *real(
                  s_complex
                  *g(s,Extension*Delta)
                  *pow(sVt00[s+smax][t],2)
                  );
        }
    }
    printf("***Intact calculation for extension %3i complete.***\n",Extension);

	_PF_DATA[0] = Sum_PF;
	_PF_DATA[1] = Sum_DPF;
	return _PF_DATA;
}

double ZNJ_ETA(int Extension,int nj){ // This calculates the exact ZNJ. The calculates include an addition nested loop over t'.
	double Sum_ZNJ_ETA = 0;
	#pragma omp parallel for reduction(+:Sum_ZNJ_ETA)
	for(int s1=-smax;s1<=smax;s1++){
		for(int t1=0;t1<tmax;t1++){
			for(int t2=0;t2<tmax;t2++){
				Sum_ZNJ_ETA = Sum_ZNJ_ETA + pow(lambda_s(s1),(N-1))
				*pow(t00_eval.get(t2),(N-nj))
				*pow(t00_eval.get(t1),(nj-1))
				*real(
					g(s1,Extension*Delta)
					*sVt00[s1+smax][t1]
					*sVt00[s1+smax][t2]
					*t_eta_t2[t1][t2]
				);
			}
		}
	}
	return Sum_ZNJ_ETA;
}

double ZNJ_XI(int Extension,int nj){ // This calculates the exact ZNJ. The calculates include an addition nested loop over t'.
	double Sum_ZNJ_XI = 0;
	#pragma omp parallel for reduction(+:Sum_ZNJ_XI)
	for(int s1=-smax;s1<=smax;s1++){
		for(int s2=-smax;s2<=smax;s2++){
			for(int t1=0;t1<tmax;t1++){
				Sum_ZNJ_XI = Sum_ZNJ_XI + pow(t00_eval.get(t1),(N-1))
				*pow(lambda_s(s1),(N-nj))
				*pow(lambda_s(s2),(nj-1))
				*real(  
					g(s1,Extension*Delta)
					*sVt00[s1+smax][t1]
					*CsVt00[s2+smax][t1]
					*s1_xi_s2[s2+smax][s1+smax]
				);
			}
		}
	}
	return Sum_ZNJ_XI;
}

vector<double> ZNJ_XI_ETA(int Extension,int nj){ // This calculates the exact ZNJ. The calculates include an addition nested loop over t'.
    
    vector<double> _ZNJ_DATA (2);
	double Sum_ZNJ_XI = 0;
    double Sum_ZNJ_ETA = 0;
#pragma omp parallel for reduction(+:Sum_ZNJ_XI,Sum_ZNJ_ETA)
	for(int s1=-smax;s1<=smax;s1++){
		for(int t1=0;t1<tmax;t1++){
            
            for(int s2=-smax;s2<=smax;s2++){
                Sum_ZNJ_XI = Sum_ZNJ_XI + pow(t00_eval.get(t1),(N-1))
				*pow(lambda_s(s1),(N-nj))
				*pow(lambda_s(s2),(nj-1))
				*real(
                      g(s1,Extension*Delta)
                      *sVt00[s1+smax][t1]
                      *CsVt00[s2+smax][t1]
                      *s1_xi_s2[s2+smax][s1+smax]
                      );
            }
			for(int t2=0;t2<tmax;t2++){
				Sum_ZNJ_ETA = Sum_ZNJ_ETA + pow(lambda_s(s1),(N-1))
				*pow(t00_eval.get(t2),(N-nj))
				*pow(t00_eval.get(t1),(nj-1))
				*real(
                      g(s1,Extension*Delta)
                      *sVt00[s1+smax][t1]
                      *sVt00[s1+smax][t2]
                      *t_eta_t2[t1][t2]
                      );
			}
            
		}
	}
    
    _ZNJ_DATA[0] = Sum_ZNJ_XI;
	_ZNJ_DATA[1] = Sum_ZNJ_ETA;

	return _ZNJ_DATA;
}



