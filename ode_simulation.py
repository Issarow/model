from numpy import *
from pylab import *

#define parameters
pi=0.02       
beta=30
mu= 0.02     
mut=0.075       
p=0.05       
x=0.21  
c=0.0053 
q=0.03  
r=0.0
y= 0.0026     
k=0.68
a=0.01    
z=0.99    
g=0.0
S=[100000]
Ln=[0]
Lr=[0]
I=[1]
R=[0]
T=[0]
N=100001
dt=1.01
#force=[0]

#write ODE equations
for i in range (1,200):
    dS=pi*N-(beta*I[i-1]/N+mu)*S[i-1]+r*R[i-1]
    
    dLn=(1-p)*beta*I[i-1]*S[i-1]/N-(c+ x*beta*I[i-1]/N + mu)*Ln[i-1]+g*R[i-1]
    
    dLr=x*(1-q)*beta*I[i-1]*Ln[i-1]/N-(y+mu)*Lr[i-1]+ z*R[i-1]
    
    dI=p*beta*I[i-1]*S[i-1]/N + (c + x*q*beta*I[i-1]/N)*Ln[i-1]+ y*Lr[i-1] -(k+mu+mut)*I[i-1]+a*R[i-1]
    
    dR=k*I[i-1]-(r+mu+a+z+g)*R[i-1]

#Identify states
    s =S[i-1]+dS*dt
    L1=Ln[i-1]+dLn*dt
    L2 =Lr[i-1]+dLr*dt
    ii =I[i-1]+dI*dt
    rr =R[i-1]+dR*dt
    #N =S+Ln+Lr+I+R
    t =T[i-1]+dt
    S.append(s)
    Ln.append(L1)
    Lr.append(L2)
    I.append(ii)
    R.append(rr)
    T.append(t)
###########################################################3
for i in range(1,200):
  #s=p*beta*I[i-1]*S[i-1]*dt/N
 #l1= c*Ln[i-1]*dt #+ 
 l1=q*beta*I[i-1]*Ln[i-1]*dt/N
  #l2 = y*Lr[i-1]*dt
print l1 #print the last one

#for i in range(0,50):
##    Force = q*x*Ln[i-1]*dt
#    df = y*Lr[i-1]*dt
#    df2  = f*Ln[i-1]*dt
##    f=force[i-1]+df*dt
#    force.append(df2)
#print sum(force)    # print in each year
##############################################################
#Plot figures
ylim(0,105000)
xlim(0,200)
xlabel('Time (years)',fontsize=18)
ylabel('Population',fontsize=18)
plot (T,S,'g',label="S",lw=4)
plot (T,Ln,'b:',label="L1",lw=4)
plot (T,Lr,'k-.',label="L2",lw=4)
plot (T,I,'r--',label="I",lw=4)
plot (T,R,'m.',label="T",lw=4)
####################################
#xlabel('FOI')
##plot (force,S,'m',label="S",lw=6)
#plot (force,Lr,'g',label="L2",lw=6)
#plot (force,Ln,'r',label="L1",lw=6)
#plot (force,I,'y',label="I",lw=6)
#plot (force,R,'b',label="R",lw=6)

#plot ((force,Lr)+(force,Ln))
################################
#print (T,I)
#grid()
legend()
show()
