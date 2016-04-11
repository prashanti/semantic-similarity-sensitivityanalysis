from __future__ import division

def load_sizes(filename,dbprofiles):
    infile=open(filename)
    infile.next()
    sizes=set()
    for line in infile:
        queryid,matchid=line.split("\t")[0],line.split("\t")[2]
        sizes.add((len(dbprofiles[queryid]),len(dbprofiles[matchid])))
    infile.close()
    print "Size len", len(sizes)
    return sizes
def create_noise_profiles(size,annotationpool,numprofiles):
    queryprofiles=[]
    randomindices=[]
    randomprofileindices=dict()
    selected=set()      
    poolsize=len(annotationpool)
    realindices=list(range(0,len(annotationpool)))
    while (len(queryprofiles)<numprofiles):
        profilesize=0
        selected=set()
        tempprofiles=[]
        while profilesize<size:
            randomindex=random.randint(0,poolsize-1)
            if (randomindex not in selected):
                tempprofiles.append(annotationpool[randomindex])
                selected.add(randomindex)
                profilesize+=1
        queryprofiles.append(tempprofiles)
    return queryprofiles

def noise_bp_asym_resnik(profile1,profile2,icdict,ancestors,bpsimilaritydict):
    finalsim=0
    bestmatchiclist=[]
    termmatchic=[]
    matchdata=[]
    for term1 in profile1:
        termmatchic=[]
        for term2 in profile2:
            termtuple=tuple(sorted((term1,term2)))
            if termtuple in bpsimilaritydict:
                termmatchic.append(bpsimilaritydict[termtuple])
            
            else:
                micaic=getmicaic(term1,term2,ancestors,icdict)
                termmatchic.append(micaic)
                bpsimilaritydict[termtuple]=micaic
        bestmatchiclist.append(np.max(termmatchic))
    
    scorelist=",".join(str(x) for x in bestmatchiclist)
    return scorelist,bpsimilaritydict


def noise_bp_sym_resnik(profile1,profile2,icdict,ancestors,bpsimilaritydict):
    finalsim=0
    bestmatchiclist1=[]
    bestmatchiclist2=[]
    termmatchic=[]
    matchdata=[]
    for term1 in profile1:
        termmatchic=[]
        for term2 in profile2:
            termtuple=tuple(sorted((term1,term2)))
            if termtuple in bpsimilaritydict:
                termmatchic.append(bpsimilaritydict[termtuple])
            
            else:
                micaic=getmicaic(term1,term2,ancestors,icdict)
                termmatchic.append(micaic)
                bpsimilaritydict[termtuple]=micaic
        
        bestmatchiclist1.append(np.max(termmatchic))
    

    termmatchic=[]
    matchdata=[]
    for term1 in profile2:
        termmatchic=[]
        for term2 in profile1:
            termtuple=tuple(sorted((term1,term2)))
            termmatchic.append(bpsimilaritydict[termtuple])
        bestmatchiclist2.append(np.max(termmatchic))
    scorelist=",".join(str(x) for x in bestmatchiclist1)+"\t"+",".join(str(x) for x in bestmatchiclist2)
    return scorelist,bpsimilaritydict
def noise_ap_resnik(profile1,profile2,icdict,ancestors,similaritydict):
    finalsim=0
    
    termmatch=[]
    for term1 in profile1:
        termmatch=[]
        for term2 in profile2:
            termtuple=tuple(sorted((term1,term2)))
            if termtuple in similaritydict:
                score=similaritydict[termtuple]
                
            else:
                
                micaic=getmicaic(term1,term2,ancestors,icdict)
                score=micaic
                similaritydict[termtuple]=score
            termmatch.append(score)
        
    
    
    scorelist=",".join(str(x) for x in termmatch)
    return scorelist,similaritydict

def noise_bp_sym_lin(profile1,profile2,icdict,ancestors,similaritydict):
    finalsim=0
    bestmatch1=[]
    bestmatch2=[]
    termmatch=[]
    for term1 in profile1:
        termmatch=[]
        for term2 in profile2:
            termtuple=tuple(sorted((term1,term2)))
            if termtuple in similaritydict:
                score=similaritydict[termtuple]
                termmatch.append(score)
            else:
                
                micaic=getmicaic(term1,term2,ancestors,icdict)
               
                score=(2*micaic)/(icdict[term1]+icdict[term2])
                similaritydict[termtuple]=score
            termmatch.append(score)
        bestmatch1.append(max(termmatch))
    
    
    termmatch=[]
    for term1 in profile2:
        termmatch=[]
        for term2 in profile1:
            termtuple=tuple(sorted((term1,term2)))
            score=similaritydict[termtuple]
            termmatch.append(score)
        bestmatch2.append(max(termmatch))
    
    scorelist=",".join(str(x) for x in bestmatch1)+"\t"+",".join(str(x) for x in bestmatch2)
    return scorelist,similaritydict
def noise_bp_asym_lin(profile1,profile2,icdict,ancestors,similaritydict):
    finalsim=0
    bestmatch1=[]
    termmatch=[]
    for term1 in profile1:
        termmatch=[]
        for term2 in profile2:
            termtuple=tuple(sorted((term1,term2)))
            if termtuple in similaritydict:
                score=similaritydict[termtuple]
                
            else:
                
                micaic=getmicaic(term1,term2,ancestors,icdict)
               
                score=(2*micaic)/(icdict[term1]+icdict[term2])
                similaritydict[termtuple]=score
            termmatch.append(score)
        bestmatch1.append(max(termmatch))
    
    
    
    scorelist=",".join(str(x) for x in bestmatch1)
    return scorelist,similaritydict
def noise_ap_lin(profile1,profile2,icdict,ancestors,similaritydict):
    finalsim=0
    
    termmatch=[]
    for term1 in profile1:
        termmatch=[]
        for term2 in profile2:
            termtuple=tuple(sorted((term1,term2)))
            if termtuple in similaritydict:
                score=similaritydict[termtuple]
                
            else:
                
                micaic=getmicaic(term1,term2,ancestors,icdict)
                score=(2*micaic)/(icdict[term1]+icdict[term2])
                similaritydict[termtuple]=score
            termmatch.append(score)
        
    
    
    scorelist=",".join(str(x) for x in termmatch)
    return scorelist,similaritydict


def noise_bp_sym_jiang(profile1,profile2,icdict,ancestors,similaritydict):
    finalsim=0
    bestmatch1=[]
    bestmatch2=[]
    termmatch=[]
    for term1 in profile1:
        termmatch=[]
        for term2 in profile2:
            termtuple=tuple(sorted((term1,term2)))
            if termtuple in similaritydict:
                score=similaritydict[termtuple]
                
            else:
                
                micaic=getmicaic(term1,term2,ancestors,icdict)
               
                score=1/(1+(icdict[term1]+icdict[term2]-(2*micaic)))
                similaritydict[termtuple]=score
            termmatch.append(score)
        bestmatch1.append(max(termmatch))
    
    
    termmatch=[]
    for term1 in profile2:
        termmatch=[]
        for term2 in profile1:
            termtuple=tuple(sorted((term1,term2)))
            score=similaritydict[termtuple]
            termmatch.append(score)
        bestmatch2.append(max(termmatch))
    
    scorelist=",".join(str(x) for x in bestmatch1)+"\t"+",".join(str(x) for x in bestmatch2)
    return scorelist,similaritydict
def noise_bp_asym_jiang(profile1,profile2,icdict,ancestors,similaritydict):
    finalsim=0
    bestmatch1=[]
    termmatch=[]
    for term1 in profile1:
        termmatch=[]
        for term2 in profile2:
            termtuple=tuple(sorted((term1,term2)))
            if termtuple in similaritydict:
                score=similaritydict[termtuple]
                termmatch.append(score)
            else:
                
                micaic=getmicaic(term1,term2,ancestors,icdict)
               
                score=1/(1+(icdict[term1]+icdict[term2]-(2*micaic)))
                similaritydict[termtuple]=score
            termmatch.append(score)
        bestmatch1.append(max(termmatch))
    
    
    
    scorelist=",".join(str(x) for x in bestmatch1)
    return scorelist,similaritydict
def noise_ap_jiang(profile1,profile2,icdict,ancestors,similaritydict):
    finalsim=0
    
    termmatch=[]
    for term1 in profile1:
        termmatch=[]
        for term2 in profile2:
            termtuple=tuple(sorted((term1,term2)))
            if termtuple in similaritydict:
                score=similaritydict[termtuple]
                
            else:
                
                micaic=getmicaic(term1,term2,ancestors,icdict)
                score=1/(1+(icdict[term1]+icdict[term2]-(2*micaic)))
                similaritydict[termtuple]=score
            termmatch.append(score)
        
    
    
    scorelist=",".join(str(x) for x in termmatch)
    return scorelist,similaritydict


def noise_bp_sym_jaccard(profile1,profile2,ancestors,similaritydict):
    finalsim=0
    bestmatchsimj1=[]
    bestmatchsimj2=[]
    termmatchsimj=[]
    for term1 in profile1:
        termmatchsimj=[]
        for term2 in profile2:
            termtuple=tuple(sorted((term1,term2)))
            if termtuple in similaritydict:
                simj=similaritydict[termtuple]
                termmatchsimj.append(simj)
            else:
                simj=getsimj(term1,term2,ancestors)
                similaritydict[termtuple]=simj
            termmatchsimj.append(simj)
        bestmatchsimj1.append(max(termmatchsimj))
    
    
    termmatchsimj=[]
    for term1 in profile2:
        termmatchsimj=[]
        for term2 in profile1:
            termtuple=tuple(sorted((term1,term2)))
            simj=similaritydict[termtuple]
            termmatchsimj.append(simj)
        bestmatchsimj2.append(max(termmatchsimj))
    
    scorelist=",".join(str(x) for x in bestmatchsimj1)+"\t"+",".join(str(x) for x in bestmatchsimj2)
    return scorelist,similaritydict
def noise_bp_asym_jaccard(profile1,profile2,ancestors,similaritydict):
    finalsim=0
    bestmatchsimj=[]
    for term1 in profile1:
        termmatchsimj=[]
        for term2 in profile2:
            termtuple=tuple(sorted((term1,term2)))
            if termtuple in similaritydict:
                simj=similaritydict[termtuple]
                
            else:
                simj=getsimj(term1,term2,ancestors)
                similaritydict[termtuple]=simj
            termmatchsimj.append(simj)
        bestmatchsimj.append(max(termmatchsimj))
    
    scorelist=",".join(str(x) for x in bestmatchsimj)
    return scorelist,similaritydict
def noise_ap_jaccard(profile1,profile2,ancestors,similaritydict):
    finalsim=0
    termmatchsimj=[]
    for term1 in profile1:
        for term2 in profile2:
            termtuple=tuple(sorted((term1,term2)))
            if termtuple in similaritydict:
                simj=similaritydict[termtuple]
            else:
                simj=getsimj(term1,term2,ancestors)
                similaritydict[termtuple]=simj
            termmatchsimj.append(simj)
        
    
    scorelist=",".join(str(x) for x in termmatchsimj)
    return scorelist,similaritydict

def noise_groupwise_jaccard(profile1,profile2,ancestors):
    ancestors1=set()
    ancestors2=set()
    for term in profile1:
        ancestors1=set.union(ancestors1,ancestors[term])
    for term in profile2:
        ancestors2=set.union(ancestors2,ancestors[term])
    common=set.intersection(ancestors1,ancestors2)
    
    if len(common) > 0:
        union=set.union(ancestors1,ancestors2)  
        simj=len(common)/len(union)

    else:
        simj=0  
    return simj
def noise_simGIC(profile1,profile2,ancestors,icdict):
    
    try:
        ancestors1=set()
        ancestors2=set()
        commonic=0
        unionic=0
        for term in profile1:
            ancestors1=set.union(ancestors1,ancestors[term])
        for term in profile2:

            ancestors2=set.union(ancestors2,ancestors[term])
        
        common=set.intersection(ancestors1,ancestors2)
        
        if len(common) > 0:
            union=set.union(ancestors1,ancestors2) 
            for term in common:
                commonic=commonic+icdict[term]
            for term in union:
                unionic=unionic+icdict[term] 
            simgic=commonic/unionic    
        else:
            simgic=0
         
        return simgic
    except Exception, e:
        print "Error in simGIC",e
    

def getsimj(term1,term2,ancestors):
    if len(set.union(ancestors[term1],ancestors[term2])) >0:
        simj=len(set.intersection(ancestors[term1],ancestors[term2]))/len(set.union(ancestors[term1],ancestors[term2]))
    else:
        simj=0
    return simj
    
def noise_comparisons(sizepair,runs,metriccombo,filename,dbprofiles,icdict,ancestors,annotationpool):
    print "here",sizepair,runs
    try:
        results=[]
        similaritydict=dict()
        
        
        if metriccombo=="bp_asym_aic_resnik" or metriccombo=="bp_asym_pic_resnik" :
            for i  in range(0,runs):
                noiseprofile1=create_noise_profiles(sizepair[0],annotationpool,1)[0]
                noiseprofile2=create_noise_profiles(sizepair[1],annotationpool,1)[0]

                scorelist,similaritydict=noise_bp_asym_resnik(noiseprofile1,noiseprofile2,icdict,ancestors,similaritydict)
                results.append(scorelist)

        if metriccombo=="bp_sym_aic_resnik" or metriccombo=="bp_sym_pic_resnik":
            for i  in range(0,runs):
                noiseprofile1=create_noise_profiles(sizepair[0],annotationpool,1)[0]
                noiseprofile2=create_noise_profiles(sizepair[1],annotationpool,1)[0]

                scorelist,similaritydict=noise_bp_sym_resnik(noiseprofile1,noiseprofile2,icdict,ancestors,similaritydict)
                results.append(scorelist)



        elif metriccombo=="aic_simGIC" or metriccombo =="pic_simGIC":
         
            for i  in range(0,runs):
                noiseprofile1=create_noise_profiles(sizepair[0],annotationpool,1)[0]
                
                noiseprofile2=create_noise_profiles(sizepair[1],annotationpool,1)[0]
                
                score=noise_simGIC(noiseprofile1,noiseprofile2,ancestors,icdict)
                
                results.append(score)



        elif metriccombo=="groupwise_jaccard":
            
         
            for i  in range(0,runs):
                noiseprofile1=create_noise_profiles(sizepair[0],annotationpool,1)[0]
                
                noiseprofile2=create_noise_profiles(sizepair[1],annotationpool,1)[0]
                
                score=noise_groupwise_jaccard(noiseprofile1,noiseprofile2,ancestors)
                
                results.append(score)

        

        elif metriccombo=="bp_sym_aic_lin" or metriccombo=="bp_sym_pic_lin":
            
            for i  in range(0,runs):
                noiseprofile1=create_noise_profiles(sizepair[0],annotationpool,1)[0]
                
                noiseprofile2=create_noise_profiles(sizepair[1],annotationpool,1)[0]
                
                scorelist,similaritydict=noise_bp_sym_lin(noiseprofile1,noiseprofile2,icdict,ancestors,similaritydict)
                
                results.append(scorelist)

        elif metriccombo=="bp_sym_aic_jiang" or metriccombo=="bp_sym_pic_jiang":
            for i  in range(0,runs):
                noiseprofile1=create_noise_profiles(sizepair[0],annotationpool,1)[0]
                noiseprofile2=create_noise_profiles(sizepair[1],annotationpool,1)[0]
                scorelist,similaritydict=noise_bp_sym_jiang(noiseprofile1,noiseprofile2,icdict,ancestors,similaritydict)
                results.append(scorelist)

        elif metriccombo=="bp_asym_aic_jiang" or metriccombo=="bp_asym_pic_jiang":
            for i  in range(0,runs):
                noiseprofile1=create_noise_profiles(sizepair[0],annotationpool,1)[0]
                noiseprofile2=create_noise_profiles(sizepair[1],annotationpool,1)[0]
                scorelist,similaritydict=noise_bp_asym_jiang(noiseprofile1,noiseprofile2,icdict,ancestors,similaritydict)
                results.append(scorelist)
        
        elif metriccombo=="bp_asym_aic_lin" or metriccombo=="bp_asym_pic_lin":
            for i  in range(0,runs):
                noiseprofile1=create_noise_profiles(sizepair[0],annotationpool,1)[0]
                noiseprofile2=create_noise_profiles(sizepair[1],annotationpool,1)[0]
                scorelist,similaritydict=noise_bp_asym_lin(noiseprofile1,noiseprofile2,icdict,ancestors,similaritydict)
                results.append(scorelist)

        


        elif metriccombo=="bp_sym_jaccard":
            for i  in range(0,runs):
                
                noiseprofile1=create_noise_profiles(sizepair[0],annotationpool,1)[0]
                noiseprofile2=create_noise_profiles(sizepair[1],annotationpool,1)[0]
                scorelist,similaritydict=noise_bp_sym_jaccard(noiseprofile1,noiseprofile2,ancestors,similaritydict)
               
                results.append(scorelist)
        

        elif metriccombo=="bp_asym_jaccard":    
            
                for i  in range(0,runs):
                    noiseprofile1=create_noise_profiles(sizepair[0],annotationpool,1)[0]
                    noiseprofile2=create_noise_profiles(sizepair[1],annotationpool,1)[0]
                    scorelist,similaritydict=noise_bp_asym_jaccard(noiseprofile1,noiseprofile2,ancestors,similaritydict)
                    results.append(scorelist)

        elif metriccombo=="ap_aic_jiang":
            for i  in range(0,runs):
                noiseprofile1=create_noise_profiles(sizepair[0],annotationpool,1)[0]
                noiseprofile2=create_noise_profiles(sizepair[1],annotationpool,1)[0]
                scorelist,similaritydict=noise_ap_jiang(noiseprofile1,noiseprofile2,icdict,ancestors,similaritydict)
                results.append(scorelist)

        elif metriccombo=="ap_aic_lin":
            for i  in range(0,runs):
                noiseprofile1=create_noise_profiles(sizepair[0],annotationpool,1)[0]
                noiseprofile2=create_noise_profiles(sizepair[1],annotationpool,1)[0]
                scorelist,similaritydict=noise_ap_lin(noiseprofile1,noiseprofile2,icdict,ancestors,similaritydict)
                results.append(scorelist)

        elif metriccombo=="ap_aic_resnik":
            for i  in range(0,runs):
                noiseprofile1=create_noise_profiles(sizepair[0],annotationpool,1)[0]
                noiseprofile2=create_noise_profiles(sizepair[1],annotationpool,1)[0]
                scorelist,similaritydict=noise_ap_resnik(noiseprofile1,noiseprofile2,icdict,ancestors,similaritydict)
                results.append(scorelist)
        
        elif metriccombo=="ap_jaccard":
            for i  in range(0,runs):
                noiseprofile1=create_noise_profiles(sizepair[0],annotationpool,1)[0]
                noiseprofile2=create_noise_profiles(sizepair[1],annotationpool,1)[0]
                scorelist,similaritydict=noise_ap_jaccard(noiseprofile1,noiseprofile2,ancestors,similaritydict)
                results.append(scorelist)


        print "Done"
        
        
        return results

    except Exception, e:
        print "Error",e
     


def load_randomprofiles(granularity):
    randomprofiles=dict()
    annotationpool=[]
    infile=open("../data/RandomProfiles2016.txt")
    for line in infile:
        profileid,annotation=line.strip().split("\t")
        if profileid not in randomprofiles:
            randomprofiles[profileid]=[]
        randomprofiles[profileid].append(annotation)
        annotationpool.append(annotation)
    infile.close()
    return randomprofiles,annotationpool

def load_ancestors(granularity):
    ancestors=dict()
    if granularity =='E':
        infile=open("../data/ESubsumers.txt")
    elif granularity=='EA':
        infile=open("../data/EASubsumers.txt")
    elif granularity =='EQ':
        infile=open("../data/EQSubsumers.txt")
    for line in infile:
        term,subsumer=line.strip().replace(">","").replace("<","").split("\t")
        if term not in ancestors:
            ancestors[term]=set()
            ancestors[term].add(term)
        if subsumer !="owl:Thing":      
            ancestors[term].add(subsumer)
    infile.close()
    return ancestors
def getmicaic(term1,term2,ancestors,icdict):
    micaic=0
    mica=""
  
    commonancestors=set.intersection(ancestors[term1],ancestors[term2])
   
    lcslist=[icdict[anc] for anc in commonancestors]
    
    
    if len(lcslist)>0:
        micaic=np.max(lcslist)     
        return micaic
    else:
        return 0
def log_results(result):
    results.append(result)

def writetofile(outfilename,results):
    outfile=open(outfilename,'w')
    if "simGIC" in outfilename or "Groupwise" in outfilename:
        noisedist=[]
        for arr in results:
            for score in arr:
                noisedist.append(score)
        print len(noisedist)
        json.dump(noisedist,outfile)
            
    else:
        for arr in results:
            for scorelist in arr:
                outfile.write(scorelist+"\n")
    outfile.close()


def main():
    resultsfile=sys.argv[1]
    metriccombo=sys.argv[2]

    if "E_Decay" in resultsfile:
        granularity='E'
    elif "EA_Decay" in resultsfile:
        granularity='EA'
    else:
        granularity ='EQ'
    print "granularity",granularity
    if "_AIC_" in resultsfile:
        icflag="AIC"
        icdict=json.load(open("../data/"+granularity+"_AnnotationIC.txt"))
    else:
        icflag ="PIC"
        icdict=json.load(open("../data/"+granularity+"_ProfileIC.txt"))
    
    



    if "simGIC" in resultsfile or "Groupwise" in resultsfile:
        outfilename=resultsfile.replace("/Decay","/Noise").replace("E_","E_NoiseDist_").replace("EA_","EA_NoiseDist_").replace("EQ_","EQ_NoiseDist_")

        #E_NoiseDist50_Decay_Quartile50_ProfileSize10__AIC_simGIC_Results.tsv
        
    else:
        outfilename=resultsfile.replace("BP/Decay","Noise").replace("AP/Decay","Noise").replace("E_","E_Noise_").replace("EA_","EA_Noise_").replace("EQ_","EQ_Noise_")
    print outfilename
    
    dbprofiles,annotationpool=load_randomprofiles(granularity)
    ancestors=load_ancestors(granularity)
    sizes=load_sizes(resultsfile,dbprofiles)
    pool = Pool(processes=6)
    runs=int(5100/len(sizes))
    
    manager = multiprocessing.Manager()
    sharedancestors = manager.dict(ancestors)
    sharedicdict = manager.dict(icdict)

    print "Loaded all requirements"
    for i in range(0,len(sizes)):
        print i
        pool.apply_async(noise_comparisons, args =(list(sizes)[i],runs,metriccombo,resultsfile,dbprofiles,sharedicdict,sharedancestors,annotationpool,), callback = log_results)       
    pool.close()
    pool.join()

    writetofile(outfilename,results)
    if "simGIC" not in resultsfile and "Groupwise" not in resultsfile: 
        os.system("python analyze_decay.py "+outfilename+" 50  Noise" )


if __name__ == "__main__":
    import os
    import json
    import random
    import pandas as pd
    from copy import deepcopy
    import multiprocessing
    from multiprocessing import Process, Pool, cpu_count
    from operator import itemgetter
    from scipy.stats import rankdata
    import numpy as np
    import math
    from scipy import stats
    import sys
    results=[]
    main()      

