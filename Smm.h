#ifndef SMM_H
#define SMM_H
#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>
#include <vector>
#define _USE_MATH_DEFINES

class Smm
{
    
public:
 
    double calculate_score();
    double score;
    void save(std::string,std::vector<double>);
    //constructor
    Smm(int,int,double,double,double);  
private:
/*std::vector<std::complex<double>> & operator=(const std::vector<std::complex<double>>& target)
{
    std::vector<std::complex<double>> tempv;
 for (int i=1;i<target.size()+1;i++)
 {
   tempv.push_back(target[i]);
 }
 return (tempv);
 
} */
/*std::vector<double> & operator=(const std::vector<double>& target)
{
    std::vector<double> tempv;
 for (int i=1;i<target.size()+1;i++)
 {
   tempv.push_back(target[i]);
 }
 return (tempv);
 
}*/

double gap;
double cep;
double q;
double c;
double wavelength; 
double freq;
double pi;
double omega;
double tau_FWHM;
double P_average;
int rep_rate;
double E_pulse;
double P_peak;
double cycle;
double tspan;
int t_point;
double dT;
double b;
std::vector<double> T;
int time_trace;
std::vector<double> t_trace;
int n_trace;
double M2;
double n0;
double eps0;
double coeff;
double f;
double D;
double NA;
double radius;
double I0;
double E0;
int tt;
std::vector<std::complex<double>> E;
std::vector<double> E_real;
std::vector<std::complex<double>> E_origin;
std::vector<std::complex<double>> I_origin;
std::complex<double> P_total;
double phi;
double c_fn;
std::vector<std::complex<double>> E_p;
std::vector<double> P;
double m_e;
std::vector<double> x;
std::vector<double> v;
std::vector<double> survive;
double a_temp;


};

#endif