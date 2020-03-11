#!/usr/bin/env python
# -*- coding: ASCII -*-

import subprocess
from subprocess import Popen
import os
##import dircache
import string
from string import join
import math
##from time import strftime
import random
from random import randint
from random import uniform
from random import gauss
from random import gammavariate
from random import betavariate
from math import sqrt
from sys import argv
import numpy as np
import datetime
import numpy.ma as ma
from scipy.optimize import minimize
from scipy.optimize import fmin
import gzip
import multiprocessing as mp

np.warnings.filterwarnings('ignore')

filenamein=argv[1]
filenameinplink=argv[2]#'SNP_CDX.freq'
FST=float(argv[3])
nbthreads=int(argv[4])


def kin(k,IBS,mask_matrix):
    k3=np.array([1-(k[0]+k[1]),k[0],k[1]])
    ll=-np.sum(ma.masked_array(np.log(np.sum(IBS*k3,axis=1)),mask=mask_matrix))
    pen=0
    if k3[0]<0:
        pen+=1
    if k3[0]>1:
        pen+=1
    if k3[1]<0:
        pen+=1
    if k3[1]>1:
        pen+=1
    if k3[2]<0:
        pen+=1
    if k3[2]>1:
        pen+=1
    if 4*k3[2]*k3[0]>k3[1]**2:
        pen+=1
    if np.isinf(ll)==True:
        pen+=1
    if pen>0:
        ll=10E10
    return ll



def GLkin(k,GL,IBS,mask_matrix):
    k3=np.array([1-(k[0]+k[1]),k[0],k[1]])
    ll=-np.sum(ma.masked_array(np.log(GL[0]*np.sum(IBS[0]*k3,axis=1)+GL[1]*np.sum(IBS[1]*k3,axis=1)+GL[2]*np.sum(IBS[2]*k3,axis=1)+GL[3]*np.sum(IBS[3]*k3,axis=1)+GL[4]*np.sum(IBS[4]*k3,axis=1)+GL[5]*np.sum(IBS[5]*k3,axis=1)+GL[6]*np.sum(IBS[6]*k3,axis=1)+GL[7]*np.sum(IBS[7]*k3,axis=1)+GL[8]*np.sum(IBS[8]*k3,axis=1)),mask=mask_matrix))
    pen=0
    if k3[0]<0:
        pen+=1
    if k3[0]>1:
        pen+=1
    if k3[1]<0:
        pen+=1
    if k3[1]>1:
        pen+=1
    if k3[2]<0:
        pen+=1
    if k3[2]>1:
        pen+=1
    if 4*k3[2]*k3[0]>k3[1]**2:
        pen+=1
    if np.isinf(ll)==True:
        pen+=1
    if pen>0:
        ll=10E10
    
    return ll

def Mij(freq,fst,i):
    return (1.0-fst)*freq+i*fst

#filenamein=argv[1]#'sim80_10K.vcf'
if filenamein[-2:]=='gz':
    file = gzip.open(filenamein)    
else:
    file = open(filenamein)

snp_dic={}
snp_index={}
snp_count=0
snp_data={}

x=file.readline()
while x<>'':
    if x[1]<>'#':
        if x[0]=='#':
            head=string.split(x[:-1],'\t')[9:]
            unrel_ind=[]
            for g in range(len(head)):
                unrel_ind.append(g)
                
        else:
            y=string.split(x[:-1])
            snp_dic[y[2]]=[y[3],y[4]]
            snp_index[snp_count]=y[2]
            snp_data[y[2]]=[y[8],y[9:]]
            snp_count+=1


    x=file.readline()


snp_index_rev={}
for i in snp_index:
    snp_index_rev[snp_index[i]]=i

###record individual GLs in an array
unrel_dic={}
unrel_dic_mask={}
for g in range(len(unrel_ind)):
    unrel_dic[g]=np.zeros((3,snp_count),dtype='float64')
    unrel_dic_mask[g]=np.zeros((snp_count),dtype='int32')
    

for g in range(snp_count):
    key,ind_data=snp_data[snp_index[g]]

    GT_key=-9
    GQ_key=-9
    PL_key=-9

    INFO=string.split(key,':')

    for gg in range(len(INFO)):
        if INFO[gg]=='GT':
            GT_key=gg
        if INFO[gg]=='GQ':
            GQ_key=gg
        if INFO[gg]=='PL':
            PL_key=gg

    for gg in range(len(ind_data)):
        PL=string.split(ind_data[gg],':')[PL_key]
        GQ=string.split(ind_data[gg],':')[GQ_key]

        if GQ=='.':
            unrel_dic[gg][:,g]=-9
            unrel_dic_mask[gg][g]=1
        elif int(GQ)<1:
            unrel_dic[gg][0][g]=-9
            unrel_dic[gg][1][g]=-9
            unrel_dic[gg][2][g]=-9
            unrel_dic_mask[gg][g]=1
        else:
            PL=np.array(string.split(PL,','),dtype='float64') 
            unrel_dic[gg][:,g]=10**(-PL/10)


print 'calculating allele frequencies'


#####New code for calculating AFs from plink'
nbSNPs=snp_count

AF=np.zeros(nbSNPs,float)
AF_mask=np.zeros(nbSNPs,float)


Popen.wait(Popen('plink2 --bfile '+filenameinplink+' --freq --out '+filenameinplink,shell=True))

file=open(filenameinplink+'.frq','r')
data=file.read()
data=string.split(data,'\n')
if data[-1]=='':
    del(data[-1])

snpfreq_dic={}

for g in range(1,len(data)):
    k=string.split(data[g])
    snpfreq_dic[k[1]]=[k[2],k[3],float(k[4])]


for g in range(nbSNPs):
    key=snp_index[g]
    temp=snpfreq_dic[key]
    ref,alt=snp_dic[snp_index[g]]
    if (ref==temp[1]) and (alt==temp[0]):
        freq=1-temp[2]
    elif (ref==temp[0]) and (alt==temp[1]):
        freq=temp[2]
    else:
        sdfsdf
    AF[g]=freq
    if freq<0.05:
        AF_mask[g]=1
    if freq>0.95:
        AF_mask[g]=1

print 'calculating prestored IBS|IBD'

##For every SNP (as each has a different allele frequency) calculate the P(IBS=x|IBD=z) for all combinations. We only do it for the 3 IBD values, i.e. assuming no inbreeding in the two individuals

##Makes a matrix to store these values. One matrix for each possible genotype combination we might observe. If we consider the reference as always p, there are 9 possible combinations.
##In a standard table like Mulligan or Thompson, some of these are collapsed (i.e. ppqq is the same as qqpp, as p has no meaning  has no meaning in that context)

#IBS0=np.zeros((nbSNPs,3),float)   #IBD0,IBD1,IBD2
ppqq=np.zeros((nbSNPs,3),float)
qqpp=np.zeros((nbSNPs,3),float)

#IBS1=np.zeros((nbSNPs,3),float)   #IBD0,IBD1,IBD2
pppq=np.zeros((nbSNPs,3),float)
pqpp=np.zeros((nbSNPs,3),float)
pqqq=np.zeros((nbSNPs,3),float)
qqpq=np.zeros((nbSNPs,3),float)

#IBS2=np.zeros((nbSNPs,3),float)   #IBD0,IBD1,IBD2
pppp=np.zeros((nbSNPs,3),float)
pqpq=np.zeros((nbSNPs,3),float)
qqqq=np.zeros((nbSNPs,3),float)


#populates matrix with actual probabilities, now using Anderson and Weir2007 Genetics values.
for g in range(len(AF)):
    p=AF[g]
    q=1.0-p

    IBD0_den=(1.0+2.0*FST)*(1.0+FST)*(1.0-FST)
    IBD1_den=(1.0+FST)*(1.0-FST)
    IBD2_den=(1.0-FST)
    
    #IBS0
    ppqq[g]=[(Mij(p,FST,1)*Mij(p,FST,0)*Mij(q,FST,1)*Mij(q,FST,0))/IBD0_den,0,0]   
    qqpp[g]=[(Mij(q,FST,1)*Mij(q,FST,0)*Mij(p,FST,1)*Mij(p,FST,0))/IBD0_den,0,0]
    
    #IBS1
    pppq[g]=[(2*Mij(p,FST,2)*Mij(p,FST,1)*Mij(p,FST,0)*Mij(q,FST,0))/IBD0_den,(Mij(p,FST,1)*Mij(p,FST,0)*Mij(q,FST,0))/IBD1_den,0]
    pqpp[g]=[(2*Mij(p,FST,2)*Mij(p,FST,1)*Mij(p,FST,0)*Mij(q,FST,0))/IBD0_den,(Mij(p,FST,1)*Mij(p,FST,0)*Mij(q,FST,0))/IBD1_den,0]
    pqqq[g]=[(2*Mij(q,FST,2)*Mij(q,FST,1)*Mij(q,FST,0)*Mij(p,FST,0))/IBD0_den,(Mij(q,FST,1)*Mij(q,FST,0)*Mij(p,FST,0))/IBD1_den,0]
    qqpq[g]=[(2*Mij(q,FST,2)*Mij(q,FST,1)*Mij(q,FST,0)*Mij(p,FST,0))/IBD0_den,(Mij(q,FST,1)*Mij(q,FST,0)*Mij(p,FST,0))/IBD1_den,0]
    
    #IBS2
    pppp[g]=[(Mij(p,FST,3)*Mij(p,FST,2)*Mij(p,FST,1)*Mij(p,FST,0))/IBD0_den,(Mij(p,FST,2)*Mij(p,FST,1)*Mij(p,FST,0))/IBD1_den,(Mij(p,FST,1)*Mij(p,FST,0))/IBD2_den]
    pqpq[g]=[(4*Mij(p,FST,1)*Mij(p,FST,0)*Mij(q,FST,1)*Mij(q,FST,0))/IBD0_den,(Mij(p,FST,0)*Mij(q,FST,0)*(Mij(p,FST,1)+Mij(q,FST,1)))/IBD1_den,(2*Mij(p,FST,0)*Mij(q,FST,0))/IBD2_den]
    qqqq[g]=[(Mij(q,FST,3)*Mij(q,FST,2)*Mij(q,FST,1)*Mij(q,FST,0))/IBD0_den,(Mij(q,FST,2)*Mij(q,FST,1)*Mij(q,FST,0))/IBD1_den,(Mij(q,FST,1)*Mij(q,FST,0))/IBD2_den]


#Store in a single matrix, one dimension for each of the 9 possible genotype combinations we might observe
    
IBS_all=np.array([ppqq,qqpp,pppq,pqpp,pqqq,qqpq,pppp,pqpq,qqqq])  #AATT,TTAA,AAAT,ATAA,ATTT,TTAT,AAAA,ATAT,TTTT



#pairwise analysis

k_combs=[]   #define parameter space i.e. what k1 and k2 coefficients (k0 is whatever is left over) we will look at. Currently set up as a grid moving 1% at a time, though bound by certain impossible constraints

k1_lis=[]
for g in range(101):
    k1_lis.append(g/100.0)

k2_lis=[]
for g in range(101):
    k2_lis.append(g/100.0)


for g in range(len(k1_lis)):
    k1=k1_lis[g]
    for gg in range(len(k2_lis)):
        k2=k2_lis[gg]       
        k0=(1.0-(k1+k2))
        if k1+k2<=1:
            if k0+k1+k2==1.0:
                if 4*k2*k0<k1**2:
                    k_combs.append([k0,k1,k2])

k_combs.sort()
k_combs=np.array(k_combs)

pw=[]   ##here we do every combinationpw
for g in range(len(head)):
    for gg in range(g+1,len(head)):
#        if string.split(head[g],'_')[0]==string.split(head[gg],'_')[0]:
        pw.append([g,gg])
                  

####workout multiprocessing batches
nb_process=len(pw)

batches=[]
for g in range(0,nb_process,nbthreads):
    batches.append(pw[g:g+nbthreads])


filenameout=filenamein+'.relateF_optim_'+filenameinplink+'_PL_LD_MP_AW_'+str(FST)
file=open(filenameout,'w') 
out='Ind1\tInd2\tZ0ag\tZ1ag\tZ2ag\tPI_HATag\tnbSNP\n'
file.write(out)
file.close()



print 'starting pairwise IBD computations\n'

print out[:-1]

#iterate through each pairwise comparison
def batch_lcmlkin(x,splits,nbSNPs,AF_mask,unrel_dic_mask,unrel_dic,filenameout,filenameinplink,snp_index_rev,IBS_all,output):

    ind1,ind2=splits[x]

    PIBS=np.zeros((nbSNPs,9),float)   ###matrix for all possible pairs of genotypes for every SNP. This will eventually store genotype likelihoods based on the calling likelihood (for example based on read depth)

    ###A matrix to denote SNPs we may want to mask 
    mask_mat=np.zeros(nbSNPs,int)
    temp_mask=AF_mask+unrel_dic_mask[ind1]+unrel_dic_mask[ind2]
    mask_mat[np.where(temp_mask>0)[0]]=1
    
    ###We now calculate the probability of observing all possible two genotype combination by multiplying their likelihoods together
    #I0=(l1[0]*l2[2])+(l1[2]*l2[0])  #AA,TT    TT,AA
    PPQQ=unrel_dic[ind1][0]*unrel_dic[ind2][2]
    QQPP=unrel_dic[ind1][2]*unrel_dic[ind2][0]

    #I1=(l1[0]*l2[1])+(l1[1]*l2[0])+(l1[1]*l2[2])+(l1[2]*l2[1])       #AA,AT    AT,AA   AT,TT   TT,AT
    PPPQ=unrel_dic[ind1][0]*unrel_dic[ind2][1]
    PQPP=unrel_dic[ind1][1]*unrel_dic[ind2][0]
    PQQQ=unrel_dic[ind1][1]*unrel_dic[ind2][2]
    QQPQ=unrel_dic[ind1][2]*unrel_dic[ind2][1]

    #I2=(l1[0]*l2[0])+(l1[1]*l2[1])+(l1[2]*l2[2])    #AA,AA    AT,AT   TT,TT
    PPPP=unrel_dic[ind1][0]*unrel_dic[ind2][0]
    PQPQ=unrel_dic[ind1][1]*unrel_dic[ind2][1]
    QQQQ=unrel_dic[ind1][2]*unrel_dic[ind2][2]

    ###Store these pairwise likelihoods in an array
    PIBS=np.zeros((snp_count,9),dtype='float64')
    PIBS[:,0]=PPQQ
    PIBS[:,1]=QQPP
    PIBS[:,2]=PPPQ
    PIBS[:,3]=PQPP
    PIBS[:,4]=PQQQ
    PIBS[:,5]=QQPQ
    PIBS[:,6]=PPPP
    PIBS[:,7]=PQPQ
    PIBS[:,8]=QQQQ  #AATT,TTAA,AAAT,ATAA,ATTT,TTAT,AAAA,ATAT,TTTT 

    fileout_snpldtest=open(filenameout+'_'+str(x)+'.snpldtest','w')
    for gg in range(len(mask_mat)):
        if mask_mat[gg]==0:
            fileout_snpldtest.write(snp_index[gg]+'\n')
    fileout_snpldtest.close()      
    
    Popen.wait(Popen('plink2 --bfile '+filenameinplink+' --extract '+filenameout+'_'+str(x)+'.snpldtest --make-bed --out '+filenameout+'_'+str(x)+'.snpldtest > dump_'+str(x),shell=True))
    Popen.wait(Popen('plink2 --bfile '+filenameout+'_'+str(x)+'.snpldtest --indep-pairwise 50 5 0.8 --out '+filenameout+'_'+str(x)+'.snpldtest > dump_'+str(x),shell=True))

    filein_snpldtest=open(filenameout+'_'+str(x)+'.snpldtest.prune.in','r')
    data=filein_snpldtest.read()
    data=string.split(data,'\n')
    if data[-1]=='':
        del(data[-1])

    Popen.wait(Popen('rm '+filenameout+'_'+str(x)+'.snpldtest*',shell=True))
    Popen.wait(Popen('rm dump_'+str(x),shell=True))

    
    snp_use=[]
    for gg in range(len(data)):
        snp_use.append(snp_index_rev[data[gg]])

    mask_mat[snp_use]=0
       
    called_SNPs=len(snp_use)


    if called_SNPs>1000:
        ###We transpose the matrices
        PIBSt=PIBS.transpose()

        ###identify the most likely genotype combination
        BestGT=PIBS.argmax(axis=1)

        BestIBS=np.zeros((nbSNPs,3),float)


        ###For each SNP, given the best genotype combination, pulls out the appropriate P(IBS|IBD) for all three IBS possibilities
        for gg in range(nbSNPs):
            BestIBS[gg]=IBS_all[BestGT[gg]][gg]


        IBS_all2=IBS_all[:,snp_use]
        BestIBS2=BestIBS[snp_use]
        PIBSt2=PIBSt[:,snp_use]
        mask_mat2=mask_mat[snp_use]

        res=[]
        for gg in range(3):
            ok=0
            while ok==0:
                if gg==0:
                    k1,k2=0.0,0.0
                else:
                    k1,k2=random.random(),random.random()
                if k1+k2<=1.0:
                    k0=1-(k1+k2)
                    if 4*k2*k0<=k1**2:
                        k=np.array([k1,k2])
                        if GLkin(k,PIBSt2,IBS_all2,mask_mat2)<>10E10:
                            ok=1
            temp=fmin(GLkin,k,args=(PIBSt2,IBS_all2,mask_mat2),xtol=0.01,ftol=0.01,maxiter=None,maxfun=None,full_output=1, disp=0, retall=0, callback=None)                  
            res.append([temp[1],temp[0]])

        out=head[ind1]+'\t'+head[ind2]

        try:
            res.sort()
            out=out+'\t'+str(round(1-(res[0][1][0]+res[0][1][1]),2))+'\t'+str(round(res[0][1][0],2))+'\t'+str(round(res[0][1][1],2))+'\t'+str(round(0.5*res[0][1][0]+res[0][1][1],3))
        except:
            out=out+'\t-9\t-9\t-9\t-9'
            
        out=out+'\t'+str(called_SNPs)+'\n'
        
        #print out[:-1]
    else:
        out=head[ind1]+'\t'+head[ind2]
        out=out+'\t-9\t-9\t-9\t-9'
        out=out+'\t'+str(called_SNPs)+'\n'
        #print out[:-1]
      
    output.put([ind1,ind2,out])
    
##    file=open(filenameout,'a') 
##    file.write(out)
##    file.close()



for G in range(len(batches)):
    nbthreads2=len(batches[G])

    ###queue for parallelism output
    output = mp.Queue()

    # Setup a list of processes
    processes = [mp.Process(target=batch_lcmlkin, args=(x,batches[G],nbSNPs,AF_mask,unrel_dic_mask,unrel_dic,filenameout,filenameinplink,snp_index_rev,IBS_all,output)) for x in range(nbthreads2)]

    # Run processes
    for p in processes:
        p.start()
     
    # Exit the completed processes
    for p in processes:
        p.join()

    results = [output.get() for p in processes]

    results.sort()

    out_all=''
    for GG in range(len(results)):
        print results[GG][2][:-1]
        out_all=out_all+results[GG][2]

    file=open(filenameout,'a') 
    file.write(out_all)
    file.close()
