#### Simulation of UID collisions with GC preference ##############
# Written by Jenna VanLiere Canzoniero using R version 3.1.2 

# Simulates results of next-generation sequencing run performed with barcoding and processed
# into consensus sequences by family.  Examines only 1 locus of mutation.  (0.3, 0.3, 0.2, 0.2)
# Assumes that nucleotides are added to UID according to probability vector 
# Returns a numerical vector with components in order:
# emutFreq1, error1, emutFreq2, error2, ncfreq, uUIDs, nreads

# This simulation was initially described in the paper 
# "The impact of collisions on the ability to detect mutant alleles using
# barcode-type next generation sequencing techniques" 
# by Jenna VanLiere Canzoniero, Karen Cravero, and Ben Ho Park
# For additional information, please see Using UIDcollision_GCpref,r by Jenna VanLiere Canzoniero

#in this version, nUID must be a factor of 4
UIDerrorGCpref<-function(nsamp,mutFreq,nUID, seed=1,conFrac=0.95,PCRcyc=25, PCRsuc=0.2, minMembers=2){
  set.seed(seed) #sets seed random num generator
  
  mutFreq<-round(nsamp*mutFreq)/nsamp #set the mutation freq such that it is whole # of sample molecules
  #create vector of sample molecules where 0 is mutant, 1 is normal
  samples<-rep(c(0,1), c(nsamp*mutFreq, nsamp-(nsamp*mutFreq))) 
  samples<-sample(samples)
  
  #vector of binomial distributed random numbers indicating num successful PCR cycles
  amp<- rbinom(nsamp,PCRcyc,PCRsuc)
  samples<-samples*amp
  nreads<-sum(amp)
  
  #Creates a vector of actual UIDs used
  # 4 possible nucleotides giving preferential weight to G and C probabilities (0.3, 0.3, 0.2, 0.2)
  #If different probabilities are desired, then adjust this vector accordingly
  # G=1, C=2, A=3, T=4
  nt<-c(1,1,1,2,2,2,3,3,4,4)
  nu<-log(nUID,base=4)
  u<-sample(nt,(nu*nsamp),replace=TRUE)
  mUID<-matrix(u ,ncol=nu, nrow=nsamp)
  UIDs<-as.numeric(apply(mUID,1,paste,collapse=""))
  
  
  #dataframe of UIDs and samples
  data<-data.frame(U=UIDs,S=samples, A=amp)
  #aggregates data by UID and generates mean of samples
  combined<-aggregate(data,list(data$U), sum)
  #remove UID families that have fewer than minimum reads
  combined<-combined[combined$A>minMembers, ]
  uUIDs<-dim(combined)[1] #number of actual UIDs used (removes those with amplification 0)
  mutant<-sum((combined$S/combined$A)<(1-conFrac)) #number of mutant consensus seq
  normal<-sum((combined$S/combined$A)>conFrac) #number of normal consensus seq
  nocall<-uUIDs-(mutant+normal) #number of nocalls
  emutFreq1<-mutant/(mutant+normal) #mutation freq estimated from data, assuming nocalls removed
  emutFreq2<-mutant/uUIDs #mut freq estimated from data, including nocalls in denom
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

# Output file
of<-"UIDerrorGCprefExampple.txt"
mutFreq<-c(0.0001, 0.001,0.01,0.1,0.2) #mutation frequency
nsamp<-c(10000,50000,100000,1e6,3e6) #number of sample molecules
taglength<-c(10,12, 14) #number of nucleotides in the UID
seed<-1:100 #number of repeats
i<-1
nr<-length(mutFreq)*length(nsamp)*length(taglength)*length(seed)
# initialized the data frame
GCpref<-data.frame(matrix(NA, nrow=nr, ncol=12))
names(GCpref)<-c("mutFreq", "nsamp", "nUID", "p2same", "expCol", "emutFreq1", "err1", "emutFreq2", "err2", "ncfreq", "uUIDs", "nreads")
for (m in mutFreq) {
  for (ns in nsamp){
    for (t in taglength){
      u<-4^t 
      p<-estP2same(ns,u)
      ec<-ecol(ns,u)
      for (s in seed){
        result<-UIDerrorGCpref(nsamp=ns,mutFreq=m,nUID=u, seed=s, conFrac=0.95)
        GCpref[i,1]<-m
        GCpref[i,2]<-ns
        GCpref[i,3]<-u
        GCpref[i,4]<-p
        GCpref[i,5]<-ec
        GCpref[i,6:12]<-result
        
        i<-i+1
      }
    }
    
  }
}

# put header on outfile
sink(of)
cat(c("GC pref 0.3 0.3 0.2 0.2 conFrac 0.95 PCRcyc 25 minMembers 2 \n"))
sink()
# Write data to outfile
write.table(GCpref, file=of, quote=FALSE, col.names=TRUE, append=TRUE)
