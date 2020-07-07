#include <vector>
#include "stdio.h"
#include "mpi.h"
 
int main(int argc,char **argv)
{
	int size,rank;
	static int max=16;
	static int task=4;
	
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
         
        double *h=NULL;
        if (rank==0){
             h=(double*)malloc(sizeof(double)*task*size);
        } 
        double *sub_h=(double*)malloc(sizeof(double)*task);
        //assert(sub_h!=NULL);
        //MPI_Scatter(h,task,MPI_DOUBLE,sub_h,task,MPI_DOUBLE,0,MPI_COMM_WORLD);
        std::cout<<"I'm rank "<<rank<<", I'm in charge of ";
        for (int i=0;i<task;i++)
          {  
             int index=(rank)*task+i; 
             sub_h[i]=index;
             std::cout<<sub_h[i]<<" ";
           }
         std::cout<<std::endl;
	double *score=NULL;
         if (rank==0)
          {
          score=(double*)malloc(sizeof(double)*task*size);
         }
        MPI_Gather(sub_h,task,MPI_DOUBLE,score,task,MPI_DOUBLE,0,MPI_COMM_WORLD);
        if (rank==0)
         {
              for (int i=1;i<task*size+1;i++)
          {std::cout<<score[i]<<std::endl;
          }
         }
       if (rank==0)
         {
           free(score);
           free(h);
         }
           free(sub_h);

          MPI_Barrier(MPI_COMM_WORLD);
          MPI_Finalize();          
return 0;
}
