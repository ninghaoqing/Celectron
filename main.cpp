#include <iostream>
#include <fstream>
#include <complex>
#include<math.h>
#include <cmath>
#include <vector>
#include "Smm.h"
//#include <mpi.h>

using namespace std;
double pi=3.1415926;
int main()
{

 Smm *smm_sample=new Smm(8,0,pi/4.0*3.0,5.0e-9,0);
 double score=smm_sample->calculate_score();
 std::cout<<score<<std::endl;
 return 0;
}
