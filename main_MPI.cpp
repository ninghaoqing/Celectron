#include <iostream>
#include <fstream>
#include <complex>
#include<math.h>
#include <cmath>
#include <vector>
#include "Smm.h"
#include <mpi.h>

using namespace std;
double pi=3.1415926;
int main()
{
 // MPI initialization
 int cep_num=4;
 int gap_num=16;
 double cep_min=0;
 double cep_max=pi;
 double gap_min=1e-10;
 double gap_max=10e-9;
 double *cep;
 double *gap;
 double *p_gap,*p_cep;
 int total_length;
 total_length=cep_num*gap_num;
 cep=(double*)malloc(sizeof(double)*cep_num);
 gap=(double*)malloc(sizeof(double)*gap_num);
 p_gap=(double*)malloc(sizeof(double)*total_length);
 p_cep=(double*)malloc(sizeof(double)*total_length);
 double d_cep;
 double d_gap;
 d_gap=(gap_max-gap_min)/(gap_num-1);
 d_cep=(cep_max-cep_min)/(cep_num-1);
 cep[0]=cep_min;
 for (int i=1;i<cep_num;i++)
 {
  cep[i]=cep[i-1]+d_cep; 
 // std::cout<<cep[i]<<std::endl;
 }
 
 gap[0]=gap_min;
 for (int i=1;i<gap_num;i++)
 {
   gap[i]=gap[i-1]+d_gap;
   //std::cout<<cep[0]<<" "<<gap[i]<<std::endl;
 }
 int index;
 index=0; 
 
 for (int i=0;i<cep_num;i++)
  {
   for (int j=0;j<gap_num;j++)
     {
      index=i*gap_num+j;
      p_gap[index]=gap[j];
      p_cep[index]=cep[i];
      
    //  std::cout<<cep[i]<<" "<<gap[j]<<" "<<sizeof(p_gap)<<" "<<p_gap[index]<<" "<<p_cep[index]<<std::endl;
      index++;
     }
  }
 
 MPI_Init(NULL ,NULL);
 int rank;
 MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 int size;
 MPI_Comm_size(MPI_COMM_WORLD, &size);
 int task;
 task=ceil(total_length/size);
 if (size<2){

 std::cout<<"less than 2 worker"<<endl;
 MPI_Abort(MPI_COMM_WORLD,1);
 }


//double *h=NULL;
      //  if (rank==0){
        //     h=(double*)malloc(sizeof(double)*task*size);
        //} 
        double *sub_score=(double*)malloc(sizeof(double)*task);
        //assert(sub_h!=NULL);
       // MPI_Scatter(h,task,MPI_DOUBLE,sub_h,task,MPI_DOUBLE,0,MPI_COMM_WORLD);
        std::cout<<"I'm rank "<<rank<<", I'm in charge of ";
        for (int i=0;i<task;i++)
          {  
             int index=(rank)*task+i; 
              Smm *smm_sample=new Smm(8,0,p_cep[index],p_gap[index],0);
              double tempscore=smm_sample->calculate_score();    
              sub_score[i]=tempscore;
              free(smm_sample);
             std::cout<<index<<" with score "<<sub_score[i]<<" ";
           }
        std::cout<<std::endl;
	double *score=NULL;
         if (rank==0)
          {
        score=(double*)malloc(sizeof(double)*task*size);
         }
        MPI_Gather(sub_score,task,MPI_DOUBLE,score,task,MPI_DOUBLE,0,MPI_COMM_WORLD);
        if (rank==0)
         {
             std::ofstream fout;
 	     fout.open("data.dat");
             fout.precision(10); 
	     for (int i=0;i<task*size;i++)
		{
 		 fout<<i<<"	"<<p_cep[i]<<"	"<<p_gap[i]<<"	"<<score[i]<<std::endl;
                 cout<<i<<"	"<<p_cep[i]<<"	"<<p_gap[i]<<"	"<<score[i]<<std::endl;
		}
  	     fout.close(); 
         }
       if (rank==0)
         {
           free(score);
           //free(h);
         }
           free(sub_score);

          MPI_Barrier(MPI_COMM_WORLD);
          MPI_Finalize();     

return 0;
}
