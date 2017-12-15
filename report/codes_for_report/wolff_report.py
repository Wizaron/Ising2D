import random
import numpy as np

def initialize(nlist,L,N):                    
    for i in xrange(0,N):
        nlist[i][0] = (i-1)%N
        nlist[i][1] = (i-L)%N
        nlist[i][2] = (i+1)%N
        nlist[i][3] = (i+L)%N   

def wolff_step(spin,nlist,N,Padd):
    
    stack = np.zeros(N)   # duzelt
    js = random.randrange(0,N)
    stack[0] = js
    sp=1
    
    oldspin=spin[js]
    newspin = -1*spin[js]
    spin[js] = newspin
    cluster_size=0
    
    while sp:
        sp = sp - 1
        current = stack[sp]
        for i in xrange (4):
            if spin[nlist[current][i]] == oldspin: 
                if random.uniform(0, 1) < Padd :
                    stack[sp] = nlist[current][i]
                    sp = sp + 1 
                    spin[nlist[current][i]] = newspin
                    cluster_size = cluster_size + 1
                    
    return cluster_size

def measure(spin,nlist,N,Nobs):
    
    observable=np.zeros(Nobs)
    
    observable[0] = np.abs((1.0* np.sum(spin) / N)   )                                                                         # magnetization 
    observable[1] = -1.0* np.sum(spin[i]* (spin[nlist[i][0]] + spin[nlist[i][1]]) for i in xrange (0,N) ) / N          # energy
    observable[2] = observable[0]**2                                                                                   # m^2
    observable[3] = observable[1]**2                                                                                   # e^2    
    
    return observable


def autocorrelation(Nsw,J,L,N,H,filename,Nobs):
    
    autocorrelation= [[] for i in xrange(Nobs)] 
    observable=load_observable1(Nsw,J,L,H,filename,Nobs)
    
    for i in xrange(0,Nobs):
        fm=np.fft.rfft(observable[i][:]-np.mean(observable[i][:]))/np.sqrt(len(observable[i][:]))
        fm2=np.abs(fm)**2
        cm=np.fft.irfft(fm2, len(observable[i][:])) 
        cm_2= cm / cm[0]
        autocorrelation[i] = cm_2
    
    np.savetxt('%s/autocorrelation/autocorrelation_fnct_magnet_J%.2f_L%i_H%.2f' %(filename,J,L,H), autocorrelation[0])    
    np.savetxt('%s/autocorrelation/autocorrelation_fnct_energy_J%.2f_L%i_H%.2f' %(filename,J,L,H), autocorrelation[1])    
    np.savetxt('%s/autocorrelation/autocorrelation_fnct_magnet2_J%.2f_L%i_H%.2f' %(filename,J,L,H), autocorrelation[2])    
    np.savetxt('%s/autocorrelation/autocorrelation_fnct_energy2_J%.2f_L%i_H%.2f' %(filename,J,L,H), autocorrelation[3])  

def correlation_time(J,L,H,filename,Nobs):
    
    correlation = np.zeros(Nobs) 
    autocorrelation= [[] for i in xrange(Nobs)]  
      
    autocorrelation[0]=np.loadtxt('%s/autocorrelation/autocorrelation_fnct_magnet_J%.2f_L%i_H%.2f' %(filename,J,L,H))    
    autocorrelation[1]=np.loadtxt('%s/autocorrelation/autocorrelation_fnct_energy_J%.2f_L%i_H%.2f' %(filename,J,L,H))    
    autocorrelation[2]=np.loadtxt('%s/autocorrelation/autocorrelation_fnct_magnet2_J%.2f_L%i_H%.2f' %(filename,J,L,H))    
    autocorrelation[3]=np.loadtxt('%s/autocorrelation/autocorrelation_fnct_energy2_J%.2f_L%i_H%.2f' %(filename,J,L,H))  

    for i in xrange(0,Nobs):
        log_cm_2 =[]
        j=0
        
        while autocorrelation[i][j] > 0 :        
            log_cm_2.append(np.log(autocorrelation[i][j]))
            j = j+1
        
        x = np.linspace(0,len(log_cm_2)-1,len(log_cm_2))
        p = np.polyfit(x,log_cm_2,1)
        tau = -1 / p[0] 
        correlation[i] = 2*int(np.ceil(tau))
        
    print correlation
    
    corr_time = int(max(correlation[i] for i in xrange(0,Nobs)))
    
    k_handle = file('%s/correlation_times' %filename, 'a')
    np.savetxt(k_handle, np.column_stack(np.append([J,H,L], corr_time)),fmt='%20.2e')
    k_handle.close()
    
    return corr_time

def ising2dmc(J,L,N,H,Nsw,max_correlation,filename,spin_read_filename,spin_write_filename,Nobs):
    #  J = J(interaction energy) * beta(= 1/kT)  
    # H external uniform magnetic field
    
    spin =  np.loadtxt('%s' %spin_read_filename) 
    nlist = np.zeros([N,4], dtype=np.int64)
    Padd= 1 - np.exp(-2*J)
    initialize(nlist,L,N)

    datafile = file('%s/observable_J%.2f_L%i_H%.2f' %(filename, J,L,H), 'a')
        
    cluster=0

    for i in xrange(1,Nsw+1):
        for j in xrange(1,max_correlation+1):
            for k in xrange(1,N+1):
                a=wolff_step(spin,nlist,N,Padd)
                cluster = cluster + a
        observable=measure(spin,nlist,N,Nobs)
        np.savetxt(datafile, observable)
        
    mean_cluster = 1.0 * cluster / (Nsw*max_correlation*N*N)
    np.savetxt('%s' %spin_write_filename,spin)
    a_handle = file('%s/cluster_size' %filename, 'a')
    np.savetxt(a_handle, np.column_stack(np.append([J, H, L], mean_cluster)) , fmt='%19.5e')
    a_handle.close()
    datafile.close()
    
# main program

f_handle = open('uncorrelated/properties_and_error','a')
f_handle.write('#            J           /        H           /       L          / Magnetic Susceptibility / Specific Heat /    Suscept Blocking   /   Heat Blocking   /  Suscept Jackknife  /  Heat Jackknife  \n')
f_handle.close()

f_handle = open('uncorrelated/mean_observable_and_error','a')
f_handle.write('#            J         /        H           /       L          /   Magnetization   /  Internal Energy  /     Magnet^2     /      Energy^2       / Magnet Blocking  / Energy Blocking  / Magnet^2 Blocking /  Energy^2 Blocking  / Magnet Jackknife  / Energy Jackknife  /  Magnet^2 Jackknife  /Energy^2 Jackknife     \n')
f_handle.close()

f_handle = open('correlated/correlation_times','a')
f_handle.write('#              J           /        H           /        L          /   Correlation Time \n')
f_handle.close()


Nobserv=4 
L=8
H=0.0
N=L*L

Temp= 1.2
J= 1.0/Temp

init_spin = np.ones(N, dtype=np.int64)  
np.savetxt('%s/spin' %('spin'),init_spin)   

Nsw_equil1=1000

ising2dmc(J,L,N,H,Nsw=Nsw_equil1,max_correlation=1,filename='equil',spin_read_filename='spin/spin',spin_write_filename='spin/spin',Nobs=Nobserv)

Nsw_corre1=2000

ising2dmc(J,L,N,H,Nsw=Nsw_corre1,max_correlation=1,filename='correlated',spin_read_filename='spin/spin',spin_write_filename='spin/spin',Nobs=Nobserv)

autocorrelation(Nsw_corre1,J,L,N,H,'correlated',Nobserv)

tau=correlation_time(J,L,H,'correlated',Nobserv)
print tau

Ndata=1000

ising2dmc(J,L,N,H,Nsw=Ndata,max_correlation=tau,filename='uncorrelated',spin_read_filename='spin/spin',spin_write_filename='spin/spin',Nobs=Nobserv)

observable_uncorre=load_observable1(Ndata,J,L,H,'uncorrelated',Nobserv)

length=len(observable_uncorre[0])

print length

analyze(Nobserv,observable_uncorre,length,J,L,N,H)

