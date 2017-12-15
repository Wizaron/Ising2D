# -*- coding: utf-8 -*-
import random
import numpy as np
import math

def initialize(nlist,L,N):              
    
    for i in xrange(0,N):
        nlist[i][0] = (i-1)%N
        nlist[i][1] = (i-L)%N
        nlist[i][2] = (i+1)%N
        nlist[i][3] = (i+L)%N   

def MCstep(init_n,spin,nlist,J,N,accept_ratio_spin,accept_ratio_n,E):
    
    epsilon = 2.0*J*N   
    
    
    new_n = random.randrange(init_n-1,init_n+2)    
    
    if new_n == 0:
        new_n=new_n+1
    
    cond = J*(epsilon - E )
    
    if new_n == init_n:
        init_n=new_n
        accept_ratio_n = accept_ratio_n +1
    
    elif new_n == init_n+1 and cond >= init_n+1:
        init_n=new_n
        accept_ratio_n = accept_ratio_n +1
        
    elif new_n == init_n-1 and cond <= init_n:
        init_n=new_n
        accept_ratio_n = accept_ratio_n +1
    
    elif new_n == init_n+1 and cond < init_n+1:
        if random.uniform(0,1) < cond/(init_n+1):
            init_n=new_n
            accept_ratio_n = accept_ratio_n +1
    
    else:
        if random.uniform(0,1) < init_n/cond:
            init_n=new_n
            accept_ratio_n = accept_ratio_n +1
    
    js = random.randrange(0,N) 
    
    MF = (spin[nlist[js][0]] + spin[nlist[js][1]] + spin[nlist[js][2]] + spin[nlist[js][3]])
    
    DE = 2*spin[js]*MF       
                                                
    if DE <= 0:         
        accept_ratio_spin = accept_ratio_spin +1
        spin[js] = -1 * spin[js]
        return spin[js], DE,accept_ratio_spin,init_n,accept_ratio_n
        
    else:   
         
        accept = pow((epsilon- (DE + E)) / (epsilon - E),init_n)
        
        if random.uniform(0, 1) < accept :   
            accept_ratio_spin = accept_ratio_spin +1
            spin[js] = -1 * spin[js]
            return spin[js], DE,accept_ratio_spin,init_n,accept_ratio_n
            
        else:
            return 0,0,accept_ratio_spin,init_n,accept_ratio_n

def ising2dmc(J,L,N,H,Nsw,max_correlation,filename,spin_read_filename,spin_write_filename,n_read_filename,n_write_filename,Nobs):
    #  J = J(interaction energy) * beta(= 1/kT)  
    # H external uniform magnetic field
    
    spin =  np.loadtxt('%s' %spin_read_filename) 
    
    init_n = int(np.loadtxt('%s' %n_read_filename)) 
    
    nlist = np.zeros([N,4], dtype=np.int64)

    accept_ratio_spin = np.zeros(1)
    accept_ratio_n = np.zeros(1)

    initialize(nlist,L,N)

    datafile = file('%s/observable_J%.2f_L%i_H%.2f' %(filename, J,L,H), 'a')
    datafile2= file('%s/n_J%.2f_L%i_H%.2f' %(filename, J,L,H), 'a')   
    datafile3= file('%s/n2_J%.2f_L%i_H%.2f' %(filename, J,L,H), 'a')
    observable=init_measure(spin,nlist,N,Nobs)
    
    E= -1.0* np.sum(spin[i]* (spin[nlist[i][0]] + spin[nlist[i][1]]) for i in xrange (0,N))  
    
    for i in xrange(1,Nsw+1):
        for j in xrange(1,max_correlation+1):
            for k in xrange(1,N+1):
                s_k, DE,accept_ratio_spin[0],init_n,accept_ratio_n[0] = MCstep(init_n,spin,nlist,J,N,accept_ratio_spin[0],accept_ratio_n[0],E)
                measure(s_k,DE,observable,N,Nobs)
                E = E + DE
        np.savetxt(datafile, observable)  
        n=[init_n]    
        n2=[init_n*init_n] 
        np.savetxt(datafile2, n)  
                                        
        np.savetxt(datafile3, n2)
        
    accept_ratio_spin[0] = 1.0*accept_ratio_spin[0] / (Nsw*N*max_correlation)
    accept_ratio_n[0] = 1.0*accept_ratio_n[0] / (Nsw*N*max_correlation)
    
    np.savetxt('%s' %spin_write_filename,spin)
    np.savetxt('%s' %n_write_filename,n)  
                                               
    
    k_handle = file('%s/acceptance_ratio_spin' %filename, 'a')
    np.savetxt(k_handle, np.column_stack(np.append([J,H,L], accept_ratio_spin)),fmt='%20.5e')
    k_handle.close()
    
    l_handle= file('%s/acceptance_ratio_n' %filename, 'a')
    np.savetxt(l_handle, np.column_stack(np.append([J,H,L], accept_ratio_n)),fmt='%20.5e')
    l_handle.close()
    
    datafile.close()

# Main Program

f_handle = open('uncorrelated/properties_and_error','a')
f_handle.write('#            J           /        H           /       L          / Magnetic Susceptibility / Specific Heat /    Suscept Blocking   /   Heat Blocking   /  Suscept Jackknife  /  Heat Jackknife  \n')
f_handle.close()

f_handle = open('uncorrelated/mean_observable_and_error','a')
f_handle.write('#            J         /        H           /       L          /   Magnetization   /  Internal Energy  /     Magnet^2     /      Energy^2       / Magnet Blocking  / Energy Blocking  / Magnet^2 Blocking /  Energy^2 Blocking  / Magnet Jackknife  / Energy Jackknife  /  Magnet^2 Jackknife  /Energy^2 Jackknife     \n')
f_handle.close()

f_handle = open('correlated/correlation_times','a')
f_handle.write('#              J           /        H           /        L          /   Correlation Time \n')
f_handle.close()

f_handle = open('equil/acceptance_ratio_spin','a')
f_handle.write('#            J           /        H           /        L          /   Acceptance Ratio \n')
f_handle.close()

f_handle = open('equil/acceptance_ratio_n','a')
f_handle.write('#            J           /        H           /        L          /   Acceptance Ratio \n')
f_handle.close()

f_handle = open('correlated/acceptance_ratio_spin','a')
f_handle.write('#            J           /        H           /        L          /   Acceptance Ratio \n')
f_handle.close()

f_handle = open('correlated/acceptance_ratio_n','a')
f_handle.write('#            J           /        H           /        L          /   Acceptance Ratio \n')
f_handle.close()

f_handle = open('uncorrelated/acceptance_ratio_spin','a')
f_handle.write('#            J           /        H           /        L          /   Acceptance Ratio \n')
f_handle.close()

f_handle = open('uncorrelated/acceptance_ratio_n','a')
f_handle.write('#            J           /        H           /        L          /   Acceptance Ratio \n')
f_handle.close()

Nobserv=4 
L=32
T=3.0
J=1/T
H=0.0
N=L*L

init_spin = np.ones(N, dtype=np.int64)  
np.savetxt('%s/spin' %('spin'),init_spin) 

init_n = [int(J*N)]
np.savetxt('%s/n' %('n'),init_n)

Nsw_equil1=3000

ising2dmc(J,L,N,H,Nsw=Nsw_equil1,max_correlation=1,filename='equil',spin_read_filename='spin/spin',spin_write_filename='spin/spin',n_read_filename='n/n',n_write_filename='n/n',Nobs=Nobserv,Num_n=num_n,Num_spin=num_spin)

Nsw_corre1=4000

ising2dmc(J,L,N,H,Nsw=Nsw_corre1,max_correlation=1,filename='correlated',spin_read_filename='spin/spin',spin_write_filename='spin/spin',n_read_filename='n/n',n_write_filename='n/n',Nobs=Nobserv,Num_n=num_n,Num_spin=num_spin)

tau=autocorrelation(Nsw_corre1,J,L,N,H,'correlated',Nobserv)

print tau

Ndata=3000

ising2dmc(J,L,N,H,Nsw=Ndata,max_correlation=tau,filename='uncorrelated',spin_read_filename='spin/spin',spin_write_filename='spin/spin',n_read_filename='n/n',n_write_filename='n/n',Nobs=Nobserv,Num_n=num_n,Num_spin=num_spin)

observable_uncorre=load_observable1(Ndata,J,L,H,'uncorrelated',Nobserv)

length=len(observable_uncorre[0])

print length

analyze(Nobserv,observable_uncorre,length,J,L,N,H)

n_uncorre=np.zeros([2,Ndata])

n_uncorre[0]=np.loadtxt('%s/n_J%.2f_L%i_H%.2f' %('uncorrelated' , J, L,H))
n_uncorre[1]=np.loadtxt('%s/n2_J%.2f_L%i_H%.2f' %('uncorrelated' , J, L,H))

length=len(observable_uncorre[0])

print length

epsilon=0
n_mean=np.zeros(2)
for i in xrange(2):
    n_mean[i]=np.mean(n_uncorre[i])

energy_n= (epsilon - 1.0*n_mean[0]/J)/N
energy2_n= (epsilon - 1.0*n_mean[0]/J)/(N*N)
c_n= 1.0*(n_mean[1] - n_mean[0]*n_mean[0] - n_mean[0])/N 


print energy_n, energy2_n, c_n



