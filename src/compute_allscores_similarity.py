from __future__ import division
def getmicaic(term1,term2,ancestors,icdict):
    micaic=0
    mica=""
    lcslist=[]
    if term1 in ancestors and term2 in ancestors:
        commonancestors=set.intersection(ancestors[term1],ancestors[term2])
        lcslist=[icdict[anc] for anc in commonancestors]
    
    
    if len(lcslist)>0:
        micaic=np.max(lcslist)     
        return micaic
    else:
        return 0




def load_ancestors(granularity):
    ancestors=dict()
    if granularity =='E':
        infile=open("../data/ESubsumers.txt")
        
    elif granularity=='EA':
        infile=open("../data/EASubsumers.txt")
        
    elif granularity =='EQ':
        infile=open("../data/EQSubsumers.txt")
        
    for line in infile:
        term,subsumer=line.strip().replace("<","").replace(">","").split("\t")
        if term not in ancestors:
            ancestors[term]=set()
            ancestors[term].add(term)
        if subsumer !="owl:Thing":      
            ancestors[term].add(subsumer)
    infile.close()
    return ancestors

def load_randomprofiles(granularity):
    randomprofiles=dict()
    annotationpool=[]
    if granularity=='E':
        infile=open("../data/RandomProfiles2016_E.txt")
    else:
        infile=open("../data/RandomProfiles2016.txt")
    for line in infile:
        profileid,annotation=line.strip().split("\t")
        if profileid not in randomprofiles:
            randomprofiles[profileid]=[]
        randomprofiles[profileid].append(annotation)
        annotationpool.append(annotation)
    infile.close()
    return randomprofiles,annotationpool

def calculate_bestpairs_symmetric_resnik(profile1,profile2,icdict,ancestors,bpsimilaritydict):
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
            
    return bestmatchiclist1,bestmatchiclist2,bpsimilaritydict

def load_hrss_precomputation(granularity,queryprofilesize,decaytype):

    

    distancedict=dict()
    micadict=dict()
    infile=open("../data/"+decaytype+"_"+granularity+"_ProfileSize_"+str(queryprofilesize)+"_DistanceHRSS.tsv")
    for line in infile:
        data=line.split("\t")
        term1,term2,distance=data[0],data[1],data[2]
        if int(distance)!=0:
            distancedict[tuple(sorted((term1,term2)))]=int(distance)
    infile.close()
    

    infile=open("../data/"+decaytype+"_"+granularity+"_MICAforHRSS_ProfileSize"+str(queryprofilesize)+".txt")

    for line in infile:
        key,mica=line.strip().split("\t")
        micadict[key]=mica
    infile.close()

    return distancedict,micadict

def load_subclasses(granularity,icdict):
    subclasses=dict()
    infile=open("../data/"+granularity+"_Subclasses.txt")
    for line in infile:

        term,subclass=line.strip().split("\t")
        
        if term not in subclasses:
            subclasses[term]=set()
            subclasses[term].add(term)
        if subclass !="owl:Nothing" and subclass in icdict:        
            subclasses[term].add(subclass)
    infile.close()
    return subclasses

def calculate_hrss(profile1,profile2,distancedict,icdict,ancestors,subclasses,micadict):
    bp_hrssscores=[]
    bp_revhrssscores=[]
    ap_hrssscores=[]
    scores=dict()
    revscores=dict()
    for term1 in profile1:
        for term2 in profile2: 
            
            pair=sorted((term1,term2))
            key=pair[0]+","+pair[1]
            mica=micadict[key]
            if mica in icdict:
                micaic=icdict[mica]
            else:
                micaic=0
            term1dist=0
            term2dist=0
            if tuple(sorted((term1,mica))) in distancedict:
                term1dist=distancedict[tuple(sorted((term1,mica)))]
            
            if tuple(sorted((term2,mica))) in distancedict:
                term2dist=distancedict[tuple(sorted((term2,mica)))]
            

            gamma=int(term1dist)+int(term2dist)
            mil1ic=mil2ic=0
            leaficlist=[]
            if term1 in subclasses:
                for term in subclasses[term1]:
                    if term in icdict:
                        leaficlist.append(icdict[term])
            
            
            
            if len(leaficlist)>0:
                mil1ic=max(leaficlist)
            else:
                mil1ic=icdict[term1]
            leaficlist=[]
            
            if term2 in subclasses:
                for term in subclasses[term2]:
                    if term in icdict:
                        leaficlist.append(icdict[term])

            if len(leaficlist)>0:
                mil2ic=max(leaficlist)
            else:
                mil2ic=icdict[term2]

    
            alpha=micaic
            
            beta=(mil1ic-icdict[term1]+mil2ic-icdict[term2])/2


            if alpha==0:
                hrss=0
            else:
                hrss=(1/(1+gamma))*(alpha/(alpha+beta))
            
            if term1 not in scores:
                scores[term1]=dict()
            if term2 not in scores[term1]:
                scores[term1][term2]=hrss

            if term2 not in revscores:
                revscores[term2]=dict()
            if term1 not in revscores[term2]:
                revscores[term2][term1]=hrss


            ap_hrssscores.append(hrss)
    for term1 in scores:
        bp_hrssscores.append(max(scores[term1].iteritems(), key=itemgetter(1))[1])

    for term2 in revscores:
        bp_revhrssscores.append(max(revscores[term2].iteritems(), key=itemgetter(1))[1])

            
    return ap_hrssscores,bp_hrssscores,bp_revhrssscores
    
def calculate_allpairs_resnik(profile1,profile2,icdict,ancestors,similaritydict):
    finalsim=0
    matchiclist=[]
    for term1 in profile1:
        for term2 in profile2:
            termtuple=tuple(sorted((term1,term2)))
            if termtuple in similaritydict:
                matchiclist.append(similaritydict[termtuple])
            
            else:
                micaic=getmicaic(term1,term2,ancestors,icdict)
                matchiclist.append(micaic)
                similaritydict[termtuple]=micaic
        
    return similaritydict,matchiclist

def calculate_allpairs_lin(profile1,profile2,icdict,ancestors,similaritydict):
    finalsim=0
    matchiclist=[]
    for term1 in profile1:
        for term2 in profile2:
            termtuple=tuple(sorted((term1,term2)))
            if termtuple in similaritydict:
                matchiclist.append(similaritydict[termtuple])
            
            else:
                micaic=getmicaic(term1,term2,ancestors,icdict)
                lin=(2*micaic)/(icdict[term1]+icdict[term2])
                matchiclist.append(lin)
                similaritydict[termtuple]=lin
        
    return similaritydict,matchiclist

def calculate_allpairs_jiang(profile1,profile2,icdict,ancestors,similaritydict):
    finalsim=0
    matchiclist=[]
    for term1 in profile1:
        for term2 in profile2:
            termtuple=tuple(sorted((term1,term2)))
            if termtuple in similaritydict:
                matchiclist.append(similaritydict[termtuple])
            
            else:
                micaic=getmicaic(term1,term2,ancestors,icdict)
                jiang=1/(1+(icdict[term1]+icdict[term2]-(2*micaic)))
                matchiclist.append(jiang)
                similaritydict[termtuple]=jiang
        
    return similaritydict,matchiclist


def calculate_allpairs_jaccard(profile1,profile2,ancestors,similaritydict):
    
    matchsimj=[]
    for term1 in profile1:
        
        for term2 in profile2:
            termtuple=tuple(sorted((term1,term2)))
            if termtuple in similaritydict:
                matchsimj.append(similaritydict[termtuple])
            else:
                simj=getsimj(term1,term2,ancestors)
                similaritydict[termtuple]=simj
                matchsimj.append(simj)
        
    return similaritydict,matchsimj



def calculate_bestpairs_asymmetric_resnik(profile1,profile2,icdict,ancestors,bpsimilaritydict):
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
    return bpsimilaritydict,bestmatchiclist

def calculate_bestpairs_symmetric_lin(profile1,profile2,icdict,ancestors,bpsimilaritydict):
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
                lin=(2*micaic)/(icdict[term1]+icdict[term2])
                termmatchic.append(lin)
                bpsimilaritydict[termtuple]=lin
        
        bestmatchiclist1.append(np.max(termmatchic))
    

    termmatchic=[]
    matchdata=[]
    for term1 in profile2:
        termmatchic=[]
        for term2 in profile1:
            termtuple=tuple(sorted((term1,term2)))
            termmatchic.append(bpsimilaritydict[termtuple])
        bestmatchiclist2.append(np.max(termmatchic))
            
    return bestmatchiclist1,bestmatchiclist2,bpsimilaritydict


def calculate_bestpairs_symmetric_jiang(profile1,profile2,icdict,ancestors,bpsimilaritydict):
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
                jiang=1/(1+(icdict[term1]+icdict[term2]-(2*micaic)))
                termmatchic.append(jiang)
                bpsimilaritydict[termtuple]=jiang
        
        bestmatchiclist1.append(np.max(termmatchic))
    

    termmatchic=[]
    matchdata=[]
    for term1 in profile2:
        termmatchic=[]
        for term2 in profile1:
            termtuple=tuple(sorted((term1,term2)))
            termmatchic.append(bpsimilaritydict[termtuple])
        bestmatchiclist2.append(np.max(termmatchic))
            
    return bestmatchiclist1,bestmatchiclist2,bpsimilaritydict

def calculate_bestpairs_asymmetric_jiang(queryprofile,profile2,icdict,ancestors,bpsimilaritydict):
    finalsim=0
    bestmatchiclist=[]
    termmatchic=[]
    matchdata=[]
    for term1 in queryprofile:
        termmatchic=[]
        for term2 in profile2:
            termtuple=tuple(sorted((term1,term2)))
            if termtuple in bpsimilaritydict:
                termmatchic.append(bpsimilaritydict[termtuple])
            
            else:
                micaic=getmicaic(term1,term2,ancestors,icdict)
                
                jiang=1/(1+(icdict[term1]+icdict[term2]-(2*micaic)))
                termmatchic.append(jiang)
                bpsimilaritydict[termtuple]=jiang
        bestmatchiclist.append(np.max(termmatchic))
    return bpsimilaritydict,bestmatchiclist


def calculate_bestpairs_asymmetric_lin(profile1,profile2,icdict,ancestors,bpsimilaritydict):
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
                lin=(2*micaic)/(icdict[term1]+icdict[term2])
                termmatchic.append(lin)
                bpsimilaritydict[termtuple]=lin
        bestmatchiclist.append(np.max(termmatchic))
    return bpsimilaritydict,bestmatchiclist


def calculate_bestpairs_asymmetric_jaccard(profile1,profile2,ancestors,similaritydict):
    finalsim=0
    bestmatchsimj=[]
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
        bestmatchsimj.append(max(termmatchsimj))
    return similaritydict,bestmatchsimj
    

def calculate_simgic(queryprofile,profile2,ancestors,icdict):
    
    ancestors1=set()
    ancestors2=set()
    commonic=0
    unionic=0
    for term in queryprofile:
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
    

def calculate_groupwise_jaccard(profile1,profile2,ancestors):
    
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

def calculate_bestpairs_symmetric_jaccard(profile1,profile2,ancestors,similaritydict):
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
    
    
    return bestmatchsimj1,bestmatchsimj2,similaritydict



def getsimj(term1,term2,ancestors):
    if len(set.union(ancestors[term1],ancestors[term2])) >0:
        simj=len(set.intersection(ancestors[term1],ancestors[term2]))/len(set.union(ancestors[term1],ancestors[term2]))
    else:
        simj=0
    return simj

def experiment_simgic(dbprofiles,ancestors,queryprofilesize,annotationpool,metric,icdict,icflag,granularity,decaytype):
    

    if decaytype=="RandomReplacement":
        decayedprofiles=json.load(open("../data/DecayedProfilesDict_Size"+str(queryprofilesize)+".txt"))
        outfile="../results/FullDistribution/RandomReplacement/"+metric+"/"+granularity+"_ProfileSize"+str(queryprofilesize)+"__"+icflag+"_"+metric+"_Results.tsv"
        out=open(outfile,'w')
        out.write("Number of annotations replaced\tQuery ID\tDatabase ID\tScore List\n")

        
    else:
        outfile="../results/FullDistribution/AncestralReplacement/"+metric+"/"+granularity+"_ProfileSize"+str(queryprofilesize)+"__"+icflag+"_"+metric+"_Results.tsv"
        decayedprofiles=json.load(open("../data/"+granularity+"_Ancestral_DecayedProfilesDict_Size"+str(queryprofilesize)+".txt"))
        out=open(outfile,'w')
        out.write("Number of decay degrees\tQuery ID\tDatabase ID\tScore List\n")



    subclassset=set()
    for queryid in decayedprofiles:
            for numreplaced in sorted([int(x) for x in decayedprofiles[queryid]]):
                print queryid,numreplaced
                queryprofile=decayedprofiles[queryid][str(numreplaced)]        
                for dbprofileid in dbprofiles:
                    profile2=dbprofiles[dbprofileid]
                    score=calculate_simgic(queryprofile,profile2,ancestors,icdict)
                    
                    out.write(str(numreplaced)+"\t"+queryid+"\t"+dbprofileid+"\t"+str(score)+"\n")
    
    out.close()
    return outfile


def experiment_groupwise_jaccard(dbprofiles,ancestors,queryprofilesize,annotationpool,metric,granularity,decaytype):
    


    if decaytype=="RandomReplacement":
        decayedprofiles=json.load(open("../data/DecayedProfilesDict_Size"+str(queryprofilesize)+".txt"))
        outfile="../results/FullDistribution/RandomReplacement/"+metric+"/"+granularity+"_ProfileSize"+str(queryprofilesize)+"_"+metric+"_Results.tsv"
        out=open(outfile,'w')
        out.write("Number of annotations replaced\tQuery ID\tDatabase ID\tScore List\n")

        
    else:
        decayedprofiles=json.load(open("../data/"+granularity+"_Ancestral_DecayedProfilesDict_Size"+str(queryprofilesize)+".txt"))
        outfile="../results/FullDistribution/AncestralReplacement/"+metric+"/"+granularity+"_ProfileSize"+str(queryprofilesize)+"_"+metric+"_Results.tsv"
        out=open(outfile,'w')
        out.write("Number of decay degrees\tQuery ID\tDatabase ID\tScore List\n")



    subclassset=set()
    for queryid in decayedprofiles:
            for numreplaced in sorted([int(x) for x in decayedprofiles[queryid]]):
                print queryid,numreplaced
                queryprofile=decayedprofiles[queryid][str(numreplaced)]        
                for dbprofileid in dbprofiles:
                    profile2=dbprofiles[dbprofileid]
                    score=calculate_groupwise_jaccard(queryprofile,profile2,ancestors)
                    
                    out.write(str(numreplaced)+"\t"+queryid+"\t"+dbprofileid+"\t"+str(score)+"\n")
    
    out.close()
    return outfile



def experiment_allpairs(dbprofiles,ancestors,queryprofilesize,annotationpool,icdict,icflag,metric,granularity,decaytype):
    similaritydict=dict()
    if decaytype=="RandomReplacement":
        decayedprofiles=json.load(open("../data/DecayedProfilesDict_Size"+str(queryprofilesize)+".txt"))
        outfile="../results/FullDistribution/RandomReplacement/"+metric+"/AP/"+granularity+"_ProfileSize"+str(queryprofilesize)+"_AP_"+icflag+ "_"+metric+"_Results.tsv"
        out=open(outfile,'w')
        out.write("Number of annotations replaced\tQuery ID\tDatabase ID\tScore List\n")
    else:
        decayedprofiles=json.load(open("../data/"+granularity+"_Ancestral_DecayedProfilesDict_Size"+str(queryprofilesize)+".txt"))
        outfile="../results/FullDistribution/AncestralReplacement/"+metric+"/AP/"+granularity+"_ProfileSize"+str(queryprofilesize)+"_AP_"+icflag+ "_"+metric+"_Results.tsv"
        out=open(outfile,'w')
        out.write("Number of decay degrees\tQuery ID\tDatabase ID\tScore List\n")
    

    for queryid in decayedprofiles:
            for numreplaced in sorted([int(x) for x in decayedprofiles[queryid]]):
                print queryid,numreplaced
                queryprofile=decayedprofiles[queryid][str(numreplaced)]        
                for dbprofileid in dbprofiles:
                    profile2=dbprofiles[dbprofileid]
                    if metric == "Resnik":
                        similaritydict,matchlist=calculate_allpairs_resnik(queryprofile,profile2,icdict,ancestors,similaritydict)
                    elif metric == "Jaccard":
                        bpsimilaritydict,matchlist=calculate_allpairs_jaccard(queryprofile,profile2,ancestors,similaritydict)
                    elif metric =="Lin":
                        similaritydict,matchlist=calculate_allpairs_lin(queryprofile,profile2,icdict,ancestors,similaritydict)
                    elif metric =="Jiang":
                        similaritydict,matchlist=calculate_allpairs_jiang(queryprofile,profile2,icdict,ancestors,similaritydict)    
                    out.write(str(numreplaced)+"\t"+queryid+"\t"+dbprofileid+"\t"+','.join(str(x) for x in matchlist)+"\n")
    out.close()
    return outfile



def experiment_bestpairs_asymmetric(dbprofiles,ancestors,queryprofilesize,annotationpool,icdict,icflag,metric,granularity,decaytype):
    bpsimilaritydict=dict()
    if decaytype=="RandomReplacement":
        decayedprofiles=json.load(open("../data/DecayedProfilesDict_Size"+str(queryprofilesize)+".txt"))
        outfile="../results/FullDistribution/RandomReplacement/"+metric+"/BP/"+granularity+"_ProfileSize"+str(queryprofilesize)+"_BPAsym_"+icflag+ "_"+metric+"_Results.tsv"
        out=open(outfile,'w')
        out.write("Number of annotations replaced\tQuery ID\tDatabase ID\tScore List\n")
    else:
        decayedprofiles=json.load(open("../data/"+granularity+"_Ancestral_DecayedProfilesDict_Size"+str(queryprofilesize)+".txt"))
        outfile="../results/FullDistribution/AncestralReplacement/"+metric+"/BP/"+granularity+"_ProfileSize"+str(queryprofilesize)+"_BPAsym_"+icflag+ "_"+metric+"_Results.tsv"
        out=open(outfile,'w')
        out.write("Number of decay degrees\tQuery ID\tDatabase ID\tScore List\n")
    

    for queryid in decayedprofiles:
            for numreplaced in sorted([int(x) for x in decayedprofiles[queryid]]):
                print queryid,numreplaced
                queryprofile=decayedprofiles[queryid][str(numreplaced)]        
                for dbprofileid in dbprofiles:
                    profile2=dbprofiles[dbprofileid]
                    if metric == "Resnik":
                        bpsimilaritydict,bestmatchlist=calculate_bestpairs_asymmetric_resnik(queryprofile,profile2,icdict,ancestors,bpsimilaritydict)
                    elif metric == "Jaccard":
                        bpsimilaritydict,bestmatchlist=calculate_bestpairs_asymmetric_jaccard(queryprofile,profile2,ancestors,bpsimilaritydict)
                    elif metric =="Lin":
                        bpsimilaritydict,bestmatchlist=calculate_bestpairs_asymmetric_lin(queryprofile,profile2,icdict,ancestors,bpsimilaritydict)
                    elif metric =="Jiang":
                        bpsimilaritydict,bestmatchlist=calculate_bestpairs_asymmetric_jiang(queryprofile,profile2,icdict,ancestors,bpsimilaritydict)    
                    out.write(str(numreplaced)+"\t"+queryid+"\t"+dbprofileid+"\t"+','.join(str(x) for x in bestmatchlist)+"\n")
    out.close()
    return outfile


def experiment_hrss(dbprofiles,queryprofilesize,ancestors,icflag,icdict,subclasses,granularity,decaytype):



    distancedict,micadict=load_hrss_precomputation(granularity,queryprofilesize,decaytype)
  
    print "loaded distance and mica"

    if decaytype=="RandomReplacement":
        decayedprofiles=json.load(open("../data/DecayedProfilesDict_Size"+str(queryprofilesize)+".txt"))
        bpsymoutfile="../results/FullDistribution/RandomReplacement/HRSS/BP/"+granularity+"_ProfileSize"+str(queryprofilesize)+"_BPSym_"+icflag+ "_HRSS_Results.tsv"
        bpasymoutfile="../results/FullDistribution/RandomReplacement/HRSS/BP/"+granularity+"_ProfileSize"+str(queryprofilesize)+"_BPAsym_"+icflag+ "_HRSS_Results.tsv"
        apoutfile="../results/FullDistribution/RandomReplacement/HRSS/AP/"+granularity+"_ProfileSize"+str(queryprofilesize)+"_AP_"+icflag+ "_HRSS_Results.tsv"

        bpsym=open(bpsymoutfile,'w')

        bpasym=open(bpasymoutfile,'w')

        ap=open(apoutfile,'w')

        ap.write("Number of annotations replaced\tQuery ID\tDatabase ID\tScore List\n")

        bpasym.write("Number of annotations replaced\tQuery ID\tDatabase ID\tScore List\n")

        bpsym.write("Number of annotations replaced\tQuery ID\tDatabase ID\tScore List1\tScoreList2\n")





    


    else:
        decayedprofiles=json.load(open("../data/"+granularity+"_Ancestral_DecayedProfilesDict_Size"+str(queryprofilesize)+".txt"))
        bpasymoutfile="../results/FullDistribution/AncestralReplacement/HRSS/BP/"+granularity+"_ProfileSize"+str(queryprofilesize)+"_BPAsym_"+icflag+ "_HRSS_Results.tsv"

        bpsymoutfile="../results/FullDistribution/AncestralReplacement/HRSS/BP/"+granularity+"_ProfileSize"+str(queryprofilesize)+"_BPSym_"+icflag+ "_HRSS_Results.tsv"

        apoutfile="../results/FullDistribution/AncestralReplacement/HRSS/AP/"+granularity+"_ProfileSize"+str(queryprofilesize)+"_AP_"+icflag+ "_HRSS_Results.tsv"



        bpsym=open(bpsymoutfile,'w')

        bpasym=open(bpasymoutfile,'w')

        ap=open(apoutfile,'w')

        ap.write("Number of decay degrees\tQuery ID\tDatabase ID\tScore List\n")

        bpasym.write("Number of decay degrees\tQuery ID\tDatabase ID\tScore List\n")

        bpsym.write("Number of decay degrees\tQuery ID\tDatabase ID\tScore List1\tScoreList2\n")

    
    for queryid in decayedprofiles:
        for numreplaced in sorted([int(x) for x in decayedprofiles[queryid]]):
            print queryid,numreplaced
      
            queryprofile=decayedprofiles[queryid][str(numreplaced)]
            for dbprofileid in dbprofiles:
                dbprofile=dbprofiles[dbprofileid]
                ap_hrss,bp_hrss,bp_revhrss=calculate_hrss(queryprofile,dbprofile,distancedict,icdict,ancestors,subclasses,micadict)

                


                ap.write(str(numreplaced)+"\t"+queryid+"\t"+dbprofileid+"\t"+','.join(str(x) for x in ap_hrss)+"\n")

                bpasym.write(str(numreplaced)+"\t"+queryid+"\t"+dbprofileid+"\t"+','.join(str(x) for x in bp_hrss)+"\n")

                bpsym.write(str(numreplaced)+"\t"+queryid+"\t"+dbprofileid+"\t"+','.join(str(x) for x in bp_hrss)+    "\t"+','.join(str(x) for x in bp_revhrss)+"\n")
                
    print "finished writing"
    distancedict=dict()
    micadict=dict()    
    ap.close()
    bpasym.close()
    bpsym.close()
    return apoutfile,bpsymoutfile,bpasymoutfile 



def experiment_bestpairs_symmetric(dbprofiles,ancestors,queryprofilesize,annotationpool,icdict,icflag,metric,granularity,decaytype):
    
    bpsimilaritydict=dict()
    if decaytype == "RandomReplacement":

        
        decayedprofiles=json.load(open("../data/DecayedProfilesDict_Size"+str(queryprofilesize)+".txt"))
        outfile="../results/FullDistribution/RandomReplacement/"+metric+"/BP/"+granularity+"_ProfileSize"+str(queryprofilesize)+"_BPSym_"+icflag+ "_"+metric+"_Results.tsv"
        out=open(outfile,'w')
        out.write("Number of annotations replaced\tQuery ID\tDatabase ID\tScore List1\tScore list2\n")


    else:
        decayedprofiles=json.load(open("../data/"+granularity+"_Ancestral_DecayedProfilesDict_Size"+str(queryprofilesize)+".txt"))
        outfile="../results/FullDistribution/AncestralReplacement/"+metric+"/BP/"+granularity+"_ProfileSize"+str(queryprofilesize)+"_BPSym_"+icflag+ "_"+metric+"_Results.tsv"
        out=open(outfile,'w')
        out.write("Number of degrees of decay\tQuery ID\tDatabase ID\tScore List1\tScore list2\n")


    
    for queryid in decayedprofiles:
            for numreplaced in sorted([int(x) for x in decayedprofiles[queryid]]):
                print queryid,numreplaced
                queryprofile=decayedprofiles[queryid][str(numreplaced)]        
                for dbprofileid in dbprofiles:
                    profile2=dbprofiles[dbprofileid]
                    if metric =="Resnik":
                        bestmatchlist1,bestmatchlist2,bpsimilaritydict=calculate_bestpairs_symmetric_resnik(queryprofile,profile2,icdict,ancestors,bpsimilaritydict)
                    elif metric=="Jaccard":
                        bestmatchlist1,bestmatchlist2,bpsimilaritydict=calculate_bestpairs_symmetric_jaccard(queryprofile,profile2,ancestors,bpsimilaritydict)
                    elif metric =="Lin":
                        bestmatchlist1,bestmatchlist2,bpsimilaritydict=calculate_bestpairs_symmetric_lin(queryprofile,profile2,icdict,ancestors,bpsimilaritydict)

                    elif metric =="Jiang":
                        bestmatchlist1,bestmatchlist2,bpsimilaritydict=calculate_bestpairs_symmetric_jiang(queryprofile,profile2,icdict,ancestors,bpsimilaritydict)
    

                    
                    out.write(str(numreplaced)+"\t"+queryid+"\t"+dbprofileid+"\t"+','.join(str(x) for x in bestmatchlist1)+"\t"+','.join(str(x) for x in bestmatchlist2)+"\n")
    out.close()
    return outfile



def run_ap(icflag,metric):
    queryprofilesize=int(sys.argv[1])
    granularity=sys.argv[2]
    decaytype=sys.argv[4]
    
    dbprofiles,annotationpool=load_randomprofiles(granularity)
    ancestors=ancestors=load_ancestors(granularity)



    if icflag=="PIC":
        icdict=json.load(open("../data/"+granularity+"_ProfileIC.txt"))
        outfile=experiment_allpairs(dbprofiles,ancestors,queryprofilesize,annotationpool,icdict,"PIC",metric,granularity,decaytype)
    elif icflag=="AIC":
        icdict=json.load(open("../data/"+granularity+"_AnnotationIC.txt"))
        outfile=experiment_allpairs(dbprofiles,ancestors,queryprofilesize,annotationpool,icdict,"AIC",metric,granularity,decaytype)
    else:
        icdict=dict()
        outfile=experiment_allpairs(dbprofiles,ancestors,queryprofilesize,annotationpool,icdict,"",metric,granularity,decaytype)
    os.system("python analyze_decay.py "+ outfile+" 50 Real")




def run_bp_sym_lin(icflag):
    queryprofilesize=int(sys.argv[1])
    granularity=sys.argv[2]
    decaytype=sys.argv[4]
    
    dbprofiles,annotationpool=load_randomprofiles(granularity)
    ancestors=ancestors=load_ancestors(granularity)



    if icflag=="PIC":
        icdict=json.load(open("../data/"+granularity+"_ProfileIC.txt"))
        outfile=experiment_bestpairs_symmetric(dbprofiles,ancestors,queryprofilesize,annotationpool,icdict,"PIC","Lin",granularity,decaytype)
    else:
        icdict=json.load(open("../data/"+granularity+"_AnnotationIC.txt"))
        outfile=experiment_bestpairs_symmetric(dbprofiles,ancestors,queryprofilesize,annotationpool,icdict,"AIC","Lin",granularity,decaytype)
    os.system("python analyze_decay.py "+ outfile+" 50 Real")


def run_bp_sym_jiang(icflag):
    queryprofilesize=int(sys.argv[1])
    granularity=sys.argv[2]
    decaytype=sys.argv[4]
    
    dbprofiles,annotationpool=load_randomprofiles(granularity)
    ancestors=ancestors=load_ancestors(granularity)



    if icflag=="PIC":
        icdict=json.load(open("../data/"+granularity+"_ProfileIC.txt"))
        outfile=experiment_bestpairs_symmetric(dbprofiles,ancestors,queryprofilesize,annotationpool,icdict,"PIC","Jiang",granularity,decaytype)
    else:
        icdict=json.load(open("../data/"+granularity+"_AnnotationIC.txt"))
        outfile=experiment_bestpairs_symmetric(dbprofiles,ancestors,queryprofilesize,annotationpool,icdict,"AIC","Jiang",granularity,decaytype)
    os.system("python analyze_decay.py "+ outfile+" 50 Real")



def run_bp_asym_jiang(icflag):
    queryprofilesize=int(sys.argv[1])
    granularity=sys.argv[2]
    decaytype=sys.argv[4]
    dbprofiles,annotationpool=load_randomprofiles(granularity)
    ancestors=ancestors=load_ancestors(granularity)



    if icflag=="PIC":
        icdict=json.load(open("../data/"+granularity+"_ProfileIC.txt"))
        outfile=experiment_bestpairs_asymmetric(dbprofiles,ancestors,queryprofilesize,annotationpool,icdict,"PIC","Jiang",granularity,decaytype)
    else:
        icdict=json.load(open("../data/"+granularity+"_AnnotationIC.txt"))
        outfile=experiment_bestpairs_asymmetric(dbprofiles,ancestors,queryprofilesize,annotationpool,icdict,"AIC","Jiang",granularity,decaytype)
    os.system("python analyze_decay.py "+ outfile+" 50 Real")






def run_bp_asym_lin(icflag):
    queryprofilesize=int(sys.argv[1])
    granularity=sys.argv[2]
    decaytype=sys.argv[4]
    dbprofiles,annotationpool=load_randomprofiles(granularity)
    ancestors=ancestors=load_ancestors(granularity)



    if icflag=="PIC":
        icdict=json.load(open("../data/"+granularity+"_ProfileIC.txt"))
        outfile=experiment_bestpairs_asymmetric(dbprofiles,ancestors,queryprofilesize,annotationpool,icdict,"PIC","Lin",granularity,decaytype)
    else:
        icdict=json.load(open("../data/"+granularity+"_AnnotationIC.txt"))
        outfile=experiment_bestpairs_asymmetric(dbprofiles,ancestors,queryprofilesize,annotationpool,icdict,"AIC","Lin",granularity,decaytype)
    os.system("python analyze_decay.py "+ outfile+" 50 Real")


def run_bp_asym_resnik(icflag):
    queryprofilesize=int(sys.argv[1])
    granularity=sys.argv[2]
    decaytype=sys.argv[4]
    dbprofiles,annotationpool=load_randomprofiles(granularity)
    ancestors=ancestors=load_ancestors(granularity)



    if icflag=="PIC":
        icdict=json.load(open("../data/"+granularity+"_ProfileIC.txt"))
        outfile=experiment_bestpairs_asymmetric(dbprofiles,ancestors,queryprofilesize,annotationpool,icdict,"PIC","Resnik",granularity,decaytype)
    else:
        icdict=json.load(open("../data/"+granularity+"_AnnotationIC.txt"))
        outfile=experiment_bestpairs_asymmetric(dbprofiles,ancestors,queryprofilesize,annotationpool,icdict,"AIC","Resnik",granularity,decaytype)
    os.system("python analyze_decay.py "+ outfile+" 50 Real")

def run_bp_sym_resnik(icflag):
    
    queryprofilesize=int(sys.argv[1])
    granularity=sys.argv[2]
    decaytype=sys.argv[4]
    dbprofiles,annotationpool=load_randomprofiles(granularity)
    ancestors=ancestors=load_ancestors(granularity)
    if icflag=="PIC":
        icdict=json.load(open("../data/"+granularity+"_ProfileIC.txt"))
        outfile=experiment_bestpairs_symmetric(dbprofiles,ancestors,queryprofilesize,annotationpool,icdict,"PIC","Resnik",granularity,decaytype)
    else:

        icdict=json.load(open("../data/"+granularity+"_AnnotationIC.txt"))
        outfile=experiment_bestpairs_symmetric(dbprofiles,ancestors,queryprofilesize,annotationpool,icdict,"AIC","Resnik",granularity,decaytype)
    os.system("python analyze_decay.py "+ outfile+" 50 Real")
    sys.exit()


def run_hrss(icflag):
    queryprofilesize=int(sys.argv[1])
    granularity=sys.argv[2]
    decaytype=sys.argv[4]
    dbprofiles,annotationpool=load_randomprofiles(granularity)
    
    ancestors=ancestors=load_ancestors(granularity)
    if icflag=="PIC":
        icdict=json.load(open("../data/"+granularity+"_ProfileIC.txt"))
        subclasses=load_subclasses(granularity,icdict)
        apoutfile,bpsymoutfile,bpasymoutfile=experiment_hrss(dbprofiles,queryprofilesize,ancestors,"PIC",icdict,subclasses,granularity,decaytype)



    else:
        icdict=json.load(open("../data/"+granularity+"_AnnotationIC.txt"))
        subclasses=load_subclasses(granularity,icdict)
        apoutfile,bpsymoutfile,bpasymoutfile=experiment_hrss(dbprofiles,queryprofilesize,ancestors,"AIC",icdict,subclasses,granularity,decaytype)
    os.system("python analyze_decay.py "+ apoutfile+" 50 Real")
    os.system("python analyze_decay.py "+ bpsymoutfile+" 50 Real")
    os.system("python analyze_decay.py "+ bpasymoutfile+" 50 Real")

    sys.exit()



def run_bp_sym_jaccard():
    queryprofilesize=int(sys.argv[1])
    granularity=sys.argv[2]
    decaytype=sys.argv[4]
    dbprofiles,annotationpool=load_randomprofiles(granularity)
    ancestors=ancestors=load_ancestors(granularity)
    
    outfile=experiment_bestpairs_symmetric(dbprofiles,ancestors,queryprofilesize,annotationpool,dict(),"","Jaccard",granularity,decaytype)
    os.system("python analyze_decay.py "+ outfile+" 50 Real")
    sys.exit()

def run_bp_asym_jaccard():    
    queryprofilesize=int(sys.argv[1])
    granularity=sys.argv[2]
    decaytype=sys.argv[4]
    dbprofiles,annotationpool=load_randomprofiles(granularity)
    ancestors=ancestors=load_ancestors(granularity)
    
    outfile=experiment_bestpairs_asymmetric(dbprofiles,ancestors,queryprofilesize,annotationpool,dict(),"","Jaccard",granularity,decaytype)
    os.system("python analyze_decay.py "+ outfile+" 50 Real")
    
    sys.exit()

def run_groupwise_jaccard():
    queryprofilesize=int(sys.argv[1])
    granularity=sys.argv[2]
    decaytype=sys.argv[4]
    dbprofiles,annotationpool=load_randomprofiles(granularity)
    ancestors=ancestors=load_ancestors(granularity)
    
    outfile=experiment_groupwise_jaccard(dbprofiles,ancestors,queryprofilesize,annotationpool,'Groupwise_Jaccard',granularity,decaytype)
    os.system("python analyze_decay.py "+ outfile+" 50 Real")
    sys.exit()


def run_simGIC(icflag):
    
        

    queryprofilesize=int(sys.argv[1])
    granularity=sys.argv[2]
    decaytype=sys.argv[4]

    if icflag=="PIC":
        icdict=json.load(open("../data/"+granularity+"_ProfileIC.txt"))
        
    else:
        icdict=json.load(open("../data/"+granularity+"_AnnotationIC.txt"))
    dbprofiles,annotationpool=load_randomprofiles(granularity)
    ancestors=ancestors=load_ancestors(granularity)
    
    outfile=experiment_simgic(dbprofiles,ancestors,queryprofilesize,annotationpool,'simGIC',icdict,icflag,granularity,decaytype)
    os.system("python analyze_decay.py "+ outfile+" 50 Real")

    sys.exit()



def main():
    metriccombo=sys.argv[3] 
    if metriccombo=="bp_asym_pic_resnik":
        run_bp_asym_resnik("PIC")  
    elif metriccombo=="bp_sym_pic_resnik":
        run_bp_sym_resnik("PIC")
    elif metriccombo=="bp_asym_aic_resnik":
        run_bp_asym_resnik("AIC")  
    elif metriccombo=="bp_sym_aic_resnik":
        run_bp_sym_resnik("AIC")
    elif metriccombo=="bp_sym_pic_lin":
        run_bp_sym_lin("PIC")
    elif metriccombo=="bp_sym_aic_lin":
        run_bp_sym_lin("AIC")
    elif metriccombo=="bp_asym_pic_lin":
        run_bp_asym_lin("PIC")
    elif metriccombo=="bp_asym_aic_lin":
        run_bp_asym_lin("AIC")


    if metriccombo=="ap_pic_resnik":
        run_ap("PIC","Resnik")  
    if metriccombo=="ap_aic_resnik":
        run_ap("AIC","Resnik")
    if metriccombo=="ap_pic_lin":
        run_ap("PIC","Lin")  
    if metriccombo=="ap_aic_lin":
        run_ap("AIC","Lin")
    if metriccombo=="ap_pic_jiang":
        run_ap("PIC","Jiang")  
    if metriccombo=="ap_aic_jiang":
        run_ap("AIC","Jiang")
    if metriccombo=="ap_jaccard":
        run_ap("","Jaccard")  
   





    elif metriccombo=="bp_sym_pic_jiang":
        run_bp_sym_jiang("PIC")
    elif metriccombo=="bp_sym_aic_jiang":
        run_bp_sym_jiang("AIC")
    elif metriccombo=="bp_asym_pic_jiang":
        run_bp_asym_jiang("PIC")
    elif metriccombo=="bp_asym_aic_jiang":
        run_bp_asym_jiang("AIC")
    
    elif metriccombo=="aic_simGIC":
        run_simGIC("AIC")    
    elif metriccombo=="pic_simGIC":
        run_simGIC("PIC")
    

    elif metriccombo=="bp_sym_jaccard":
        run_bp_sym_jaccard()
    elif metriccombo=="bp_asym_jaccard":    
        run_bp_asym_jaccard()
    elif metriccombo=="groupwise_jaccard":
        run_groupwise_jaccard()
    elif metriccombo=="pic_hrss":
        run_hrss("PIC")
    elif metriccombo=="aic_hrss":
        run_hrss("AIC")
    
    


if __name__ == "__main__":
    import os
    import json
    import pandas as pd
    from copy import deepcopy
    from operator import itemgetter
    from scipy.stats import rankdata
    import numpy as np
    import math
    from scipy import stats
    import sys
    main()    