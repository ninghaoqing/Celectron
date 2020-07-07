#include "Smm.h"
#include<vector>
#include<math.h>
void Smm::save(std::string name,std::vector<double> data )
{

std::ofstream fout;
fout.open(name);
for (int i=0;i<data.size();i++)
{
 fout<<data[i]<<'\n';

}
fout.close();



}
Smm::Smm(int time1,int time2,double cep,double gap,double V_DC)
{
      
        //this->pi=3.14159265358979323846;
        this->pi=3.1415926;
       this->gap=gap;
        this->cep=cep;
        //unit charge
        this->q=1.602e-19;
        //electron mass
        this->m_e=9.109e-31;
        //light speed
        this->c=3.0e8;
        // wavelength
        this->wavelength=800e-9;
        // frequency
        this->freq=this->c/this->wavelength;
        //pi
        //this->pi=M_PI;
        // angular frequency
        this->omega=2*this->pi*this->freq;
        // Full width at half maximum
        this->tau_FWHM=6e-15;
        // average std::pow 
        this->P_average=3.5-3;  //3.5 mW
        // repetition rate of laser
        this->rep_rate=90e6;  //90 Mhz
        // this->rep_rate=4e3
        //  
        // pulse energy 
        this->E_pulse=this->P_average*1/this->rep_rate; 
        //peak std::pow
        this->P_peak=this->E_pulse*0.94/this->tau_FWHM;
        //this->t_point=1001
        // optical cycle
        this->cycle=1/this->freq;
        //size of the pulse simulation box 
        this->tspan=20*this->tau_FWHM;
        //size of the points for pulse simulation
        this->t_point=8192*time1;
        double temppoint;
        //time interval of pulse
        this->dT=this->tspan/(t_point-1);
        //  time array for pulse
        temppoint=-tspan/2+dT;
        this->T=std::vector<double>(t_point,0);
        this->t_trace=std::vector<double>(t_point,0);
        for (int i=0;i<t_point;i++)
        {
          T[i]=temppoint;
          temppoint=temppoint+this->dT;

        }
        //this->T=np.linspace(-this->tspan/2,this->tspan/2,num=this->t_point);
        // ratio between total simulation time/ pulse simulation time 
      
        //total time box
        this->n_trace=this->t_trace.size();            //DDDDDDDDDDD
        // beam quality
        this->M2=1.0;
        // refractive index
        this->n0=1.0;
        //dielectric constant
        this->eps0=8.854e-12;
        //std::pow to intensity
        this->coeff=this->c*this->n0*this->eps0/2;
        // focal length
        this->f=30e-3;
        // beam diameter before lens
        this->D=4e-3;
        //numerical aperture
        this->NA=this->D/2.0/this->f;
        // beam radius after lens
        this->radius=2.0*std::pow(this->M2,2)/this->pi*this->wavelength/this->NA/10.0;
        // peak intensity  
        this->I0=this->P_peak/(this->pi*std::pow(this->radius,2)); 
        // peak field strength 
        this->E0=std::sqrt(this->I0*(1.0/this->coeff)); 
        // cep=this->pi/2
        // cep=this->pi/4*3
        this->P_total=0.0;
               this->b=6.83;
        this->phi=5.5;
        this->c_fn=double(1.0/60000);
        this->E_origin=std::vector<std::complex<double>>(t_point,0);
        this->I_origin=std::vector<std::complex<double>>(t_point,0);
        this->E=std::vector<std::complex<double>>(t_point,0);
        this->E_p=std::vector<std::complex<double>>(t_point,0);
        this->P=std::vector<double>(t_point,0);

        for (int i=0;i<this->E_origin.size();i++)
        {this->E_origin[i]=this->E0*std::exp(std::complex<double>(0.0, 1.0)*this->omega*this->T[i]+std::complex<double>(0.0, 1.0)*cep)*std::exp(-2.0*log(2.0)*std::pow(this->T[i]/this->tau_FWHM,2));
        this->I_origin[i]=this->coeff*std::pow(abs(this->E_origin[i]),2);
          this->P_total=this->P_total+this->I_origin[i]*this->pi*std::pow(this->radius,2)*this->dT;
         

             //enhancement
          this->E[i]=this->E_origin[i]*std::complex<double>(2.0,0);
          this->E_p[i]=this->E[i]*1e-9+V_DC;
          if (abs(this->E_p[i])>1)
                {
                  this->P[i]=this->c_fn*(std::pow(std::real(this->E_p[i]),2))*std::exp(this->b*-std::pow(this->phi,(3.0/2))/abs(this->E_p[i]));
                  //std::cout<<c_fn<<std::endl;
                  }
            else{this->P[i]=0;}
            E_real.push_back(real(E[i]));
        }
        //save("P.txt",P);
        //save("E.txt",E_real);
        //std::cout<< this->E0*std::exp(std::complex<double>(0.0, 1.0)*this->omega*this->T[60000]/*+std::complex<double>(0.0, 1.0)*cep*/)<<std::endl;
        
        //std::cout<<P[60000]<<std::endl;
        // this->E0=this->E0*20 
        
// this->E=this->E*0
// plot(this->T,std::real(this->E))
 
// this->E=10
// V_DC=0
       
       
// this->E_p=V_DC*ones(size(this->E))
// plot(this->E_p)
//plot(this->P)
//figure
// this->P=this->b*-this->phi^(3/2)./std::real(this->E)
// display(this->P)
// this->P=std::real(this->E).^(-1)
// hold on
// L_total=5e-9
// L_point=this->t_point
// L=linspace(0,L_total,l_point)
      
          this->x=std::vector<double>(T.size(),0);
          this->v=std::vector<double>(T.size(),0);
          this->survive=std::vector<double>(T.size(),0);
// rho_x(1)=1
        //print(this->survive) 
        this->survive[0]=1;
        this->a_temp=this->q*std::real(this->E[0])/this->m_e;
        tt=1;
        int i=tt;
        if (this->a_temp>0)
             {this->x[i]=-gap/2;}
        else{
             this->x[i]=gap/2;
             }
        this->survive[i]=1;
        this->score=0;            
}

double Smm::calculate_score()
{
for (tt=1;tt<t_point;tt++)
  { //    display(tt)
       a_temp=q*std::real(E[tt])/m_e; 
       for (int i=0;i<t_point;i++)
       { 
//         display(survive(i,tt));
          if (survive[i]>0)
//           display("Yes");
               
            {x[i]=x[i]+v[i]*dT;
            v[i]=v[i]+a_temp*dT;
            if (abs(x[i])>gap/2)
             {    //std::cout<<P[i]<<" "<<dT<<std::endl;   
                survive[i]=0;
                if (x[i]>gap/2)
                 { score=score-P[i]*dT;
                 }
                else
                {
                  score=score+P[i]*dT;  
                }
             }
            else
            {
             survive[i]=1;
             }
            }
       }
     int i=tt;
     //std::cout<<i<<std::endl;
     if (a_temp>0)
     {
         x[i]=-gap/2;
     }
     else
     {
         x[i]=gap/2;
     }
      survive[i]=1;
   
  }
       for (int i=0;i<t_point;i++)
        {
           if (survive[i]>0)
           {
               if (v[i]!=0)
               {
                   if (v[i]>0)
                   { 
                       score =score-P[i]*dT;
                   }
                   else
                   { 
                       score=score+P[i]*dT;
                   }
               }
           }
        } 
return score;
}
