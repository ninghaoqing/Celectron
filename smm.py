import numpy as np
import matplotlib.pyplot as plt


class SMM(object):
    def __init__(self,time1,time2,cep,gap,V_DC):
        self.gap=gap
        self.cep=cep
        #unit charge
        self.q=1.602e-19
        #light speed
        self.c=3e8
        # wavelength
        self.wavelength=800e-9
        # frequency
        self.freq=self.c/self.wavelength
        #pi
        self.pi=np.pi
        # angular frequency
        self.omega=2*self.pi*self.freq
        # Full width at half maximum
        self.tau_FWHM=6e-15
        # average np.power 
        self.P_average=3.5-3  #3.5 mW
        # repetition rate of laser
        self.rep_rate=90e6  #90 Mhz
        # self.rep_rate=4e3
        #  
        # pulse energy 
        self.E_pulse=self.P_average*1/self.rep_rate 
        #peak np.power
        self.P_peak=self.E_pulse*0.94/self.tau_FWHM
        #self.t_point=1001
        # optical cycle
        self.cycle=1/self.freq
        #size of the pulse simulation box 
        self.tspan=20*self.tau_FWHM
        #size of the points for pulse simulation
        self.t_point=8192*time1
        #time interval of pulse
        self.dT=self.tspan/self.t_point
        #  time array for pulse
        self.T=np.linspace(-self.tspan/2,self.tspan/2,num=self.t_point)
        # ratio between total simulation time/ pulse simulation time 
        self.time_trace=time2
        #total time box
        self.t_trace=np.linspace(-self.tspan/2,self.tspan*(self.time_trace-1)*self.tspan+self.tspan/2,num=self.time_trace*self.t_point)
        #
        self.n_trace=self.t_trace.size            #DDDDDDDDDDD
        # beam quality
        self.M2=1
        # refractive index
        self.n0=1
        #dielectric constant
        self.eps0=8.854e-12
        #np.power to intensity
        self.coeff=self.c*self.n0*self.eps0/2
        # focal length
        self.f=30e-3
        # beam diameter before lens
        self.D=4e-3
        #numerical aperture
        self.NA=self.D/2/self.f
        # beam radius after lens
        self.radius=2*np.power(self.M2,2)/self.pi*self.wavelength/self.NA/10 
        # peak intensity  
        self.I0=self.P_peak/(self.pi*np.power(self.radius,2)) 
        # peak field strength 
        self.E0=np.sqrt(self.I0*(1/self.coeff)) 
        # cep=self.pi/2
        # cep=self.pi/4*3
        self.E_origin=self.E0*np.exp(1j*self.omega*self.T+1j*cep)*np.exp(-2*np.log(2)*np.power(self.T/self.tau_FWHM,2))
        self.I_origin=self.coeff*np.power(abs(self.E_origin),2)
        self.P_total=sum(self.I_origin)*(self.pi*np.power(self.radius,2))*self.dT
        # self.E0=self.E0*20 
        self.E=self.E_origin*2
# self.E=self.E*0
# plot(self.T,np.real(self.E))
        self.b=6.83
        self.phi=5.5
        self.c_fn=1/60000
# self.E=10
# V_DC=0
        self.E_p=self.E*1e-9+V_DC
        self.P=np.zeros(self.E_p.size)
# self.E_p=V_DC*ones(size(self.E))
# plot(self.E_p)
        for i in range(1,self.E.size):
            if abs(self.E_p[i])>1:
                self.P[i]=self.c_fn*(np.power(np.real(self.E_p[i]),2))*np.exp(self.b*-np.power(self.phi,(3/2))/abs(self.E_p[i]))
            else:
                self.P[i]=0
#plot(self.P)
#figure
# self.P=self.b*-self.phi^(3/2)./np.real(self.E)
# display(self.P)
# self.P=np.real(self.E).^(-1)
# hold on
# L_total=5e-9
# L_point=self.t_point
# L=linspace(0,L_total,l_point)
        self.m_e=9.109e-31
        self.x=np.zeros(self.T.size)
        self.v=np.zeros(self.T.size)
        self.survive=np.zeros((self.T.size))
# rho_x(1)=1
        #print(self.survive) 
        self.survive[1]=1
        self.a_temp=self.q*np.real(self.E[1])/self.m_e
        tt=1
        i=tt
        if self.a_temp>0:
             self.x[i]=-2.5e-9
        else:
             self.x[i]=2.5e-9
        self.survive[i]=1
        self.score=0

    def calculate_score(self):
        for tt in range(2,self.n_trace):
               print(tt)
               if tt<self.T.size :
#     rho[i]=1
                     self.a_temp=self.q*np.real(self.E[tt])/self.m_e    
               else:
                     self.a_temp=0
#     self.survive[i][tt]=True
               for i in range(1,self.T.size):    
#          display(self.survive[i][tt])
                    if (self.survive[i]>0):
#           display("Yes")
                     self.x[i]=self.x[i]+self.v[i]*self.dT
                     self.v[i]=self.v[i]+self.a_temp*self.dT
                    if abs(self.x[i])>self.gap/2:
                       self.survive[i]=0
                       if (self.x[i]>=self.gap/2):
                            self.score=self.score-self.P[i]*self.dT
                       else:
                            self.score=self.score+self.P[i]*self.dT             
                    else:
                     self.survive[i]=1
        i=tt
        if self.a_temp>0:
            self.x[i]=-self.gap/2
        else:
            self.x[i]=self.gap/2
        self.survive[i]=1
        return self.score
#      self.v[i][tt]=self.v(i,tt-1)+self.a_temp*self.dT
      
#  plot(self.x(1,self.T),self.T,'self.x')
#   hold on
#   pause

#display(self.score/self.dT)
#  L=max(max(self.x))
#  R=min(min(self.x))
#  resolution=8192
#  zz=np.zeros(resolution,self.T.size)
# # zz=np.zeros(10,10)
# for i=1:1:self.T.size
#     for tt=1:1:self.T.size
#         if self.survive[i][tt]
# #         plot(self.x[i][tt],tt,'self.x')
#          tempx=round(((self.x[i][tt]+(L-R)/2)/(L-R))*resolution)+1
#          self.x[i][tt]=tempx
# #          display(tempx)
# #        tempx=abs(round(self.x[i][tt]*1e9))+1
#          zz(tempx,tt)=self.v[i][tt]
# #          display(tempx)
# #          display(zz(tempx,tt))
# #          pause
# #         hold on 
#         end
#     end
# end
# # plot(self.x(1,:))
# self.x=linspace(-2.5e-9,2.5e-9,resolution)
#  imagesc(self.T,self.x,zz)
#  polarmap(1024)
# figure
# # yyaxis right
# plot(self.T,np.real(self.E))




