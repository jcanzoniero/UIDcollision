#### Simulation of UID collisions from discrete Uniform distribution ##############
# Written by Jenna VanLiere Canzoniero using R version 3.1.2 (2014-10-31)

# Simulates results of next-generation sequencing run performed with barcoding and processed
# into consensus sequences by family.  Examines only 1 locus of mutation.  
# Returns a numerical vector with components in order:
# emutFreq1, error1, emutFreq2, error2, ncfreq, uUIDs, nreads

# This simulation was initially described in the paper 
# "The impact of collisions on the ability to detect mutant alleles using
# barcode-type next generation sequencing techniques" 
# by Jenna VanLiere Canzoniero, Karen Cravero, and Ben Ho Park
# For further information, please see Using UIDcollision_Uniform,r by Jenna VanLiere Canzoniero


UIDerror<-function(nsamp,mutFreq,nUID,seed=1,conFrac=0.95,PCRcyc=25, PCRsuc=0.5){
  set.seed(seed) #sets seed random num generator
  
  mutFreq<-round(nsamp*mutFreq)/nsamp #set the mutation freq such that it is whole # of sample molecules
  #create vector of sample molecules where 0 is mutant, 1 is normal
  samples<-rep(c(0,1), c(nsamp*mutFreq, nsamp-(nsamp*mutFreq))) 
  #mix up the order
  samples<-sample(samples)
  
  #vector of binomial distrubuted random numbers indicating number of successful PCR cycles
  amp<- rbinom(nsamp,PCRcyc,PCRsuc)
  # PCR amplification
  samples<-samples*amp
  #total number of reads for sequencing
  nreads<-sum(amp)
  
  #vector of actual UIDs used, sampled with replacement (discrete uniform distribution)
  UIDs<-sample(1:nUID,nsamp,replace=TRUE)
  
  #dataframe of UIDs and samples
  data<-data.frame(U=UIDs,S=samples, A=amp)
  #aggregates data by UID and generates mean of samples
  combined<-aggregate(data,list(data$U), sum)
  uUIDs<-sum(combined$A!=0) #number of actual UIDs used (removes those with amplification 0)
  mutant<-sum((combined$S/combined$A)<(1-conFrac), na.rm=TRUE) #number of mutant consensus seq
  normal<-sum((combined$S/combined$A)>conFrac, na.rm=TRUE) #number of normal consensus seq
  nocall<-uUIDs-(mutant+normal) #number of nocalls
  emutFreq1<-mutant/(mutant+normal) #mutation freq estimated from data, assuming nocalls removed
  emutFreq2<-mutant/uUIDs #mutation freq estimated from data, nocalls included in denominator
  ncfreq<-nocall/uUIDs #no call frequency
  error1<-abs(emutFreq1-mutFreq)/mutFreq #error for estimated mutation freq 
  error2<-abs(emutFreq2-mutFreq)/mutFreq 
  return (c(emutFreq1, error1, emutFreq2, error2, ncfreq, uUIDs, nreads))
  
}


##### ADDITIONAL FUNCTIONS ####

#function to estimate probablity that at least 2 ms have the same u (i.e. at least 2 people
#have same birthday, or at least 2 molecules have same UID)
estP2same<-function(m,u){1-exp(-m*(m-1)/(2*u))}

#expected number of samples that will have a non-unique UID in n samples from d UIDs
ecol<-function(n,d){n-d+(d*((d-1)/d)^n)}


############### EXAMPLE RUN #############################
#Below is example of how above function was used to generate simulation data described in paper.
# The below call results in an output file with header line describing input parameters
#and then table with columns
#mutFreq nsamp nUID p2same expCol emutFreq1 err1 emutFreq2 err2 ncfreq uUIDs nreads

#Name of output file
of<-"UniformSimExample.txt"
#range of actual mutation frequencies
mutFreq<-c(0.0001, 0.001,0.01,0.1,0.2,0.5)
# range of number of sample molecules
nsamp<-c(10000,25000,100000)
# range of number of UIDs
nUID<-c(4^10, 4^12)
# each combination of above parameters will be simulated 100 times, with seed set to the simulation count
seed<-1:100
i<-1
#Calculate the number of rows
nr<-length(mutFreq)*length(nsamp)*length(nUID)*length(seed)
# initialized the data frame
out<-data.frame(matrix(NA, nrow=nr, ncol=12))
names(out)<-c("mutFreq", "nsamp", "nUID", "p2same", "expCol", "emutFreq1", "err1", "emutFreq2", "err2", "ncfreq", "uUIDs", "nreads")

# For each combination of parameters and each seed/repetition, call UIDerror to do simulation
# then put results into out data frame
for (m in mutFreq) {
  for (ns in nsamp){
    for (u in nUID){
      p<-estP2same(ns,u)
      ec<-ecol(ns,u)
      for (s in seed){
        result<-UIDerror(nsamp=ns,mutFreq=m,nUID=u,seed=s, PCRsuc=0.25, conFrac=0.95)
        out[i,1]<-m
        out[i,2]<-ns
        out[i,3]<-u
        out[i,4]<-p
        out[i,5]<-ec
        out[i,6:12]<-result
        
        i<-i+1
      }
    }
    
  }
}
# put header on outfile
sink(out)
cat(c("Uniform PCRsuc 0.25 conFrac 0.95 \n"))
sink()
# write data to outfile
write.table(out, file=of, quote=FALSE, col.names=TRUE, append=FALSE)

