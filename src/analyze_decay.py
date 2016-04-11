from __future__ import division


def compute_decay(filename,quartile):
    infile=open(filename)
    print filename
    if "simGIC" in filename or "Groupwise" in filename:
        outfile=open(filename.replace("E","Decay/E").replace("ProfileSize","Decay_ProfileSize"),'w')
    else:
        outfile=open(filename.replace("E","Decay/E").replace("ProfileSize","Decay_Quartile"+str(quartile)+"_ProfileSize"),'w')
    queryids=set()
    

    outfile.write("QueryID\tNumber of Annotations Replaced\tBest Match\tSimilarity Score\n")
    matchdata=dict()
    infile.next()
    for line in infile:
        data=line.strip().split("\t")
        if len(data)==4:
            numreplaced,queryid,matchid,scorelist=int(data[0]),data[1],data[2],[float(x) for x in data[3].split(",")]
            
            score = np.percentile(scorelist, quartile)
            queryids.add(queryid)
        else:
            numreplaced,queryid,matchid,scorelist1,scorelist2=int(data[0]),data[1],data[2],[float(x) for x in data[3].split(",")],[float(x) for x in data[4].split(",")]
            queryids.add(queryid)
            score1 = np.percentile(scorelist1, quartile)
            score2 = np.percentile(scorelist2, quartile)
            score=np.mean([score1,score2])
        if numreplaced not in matchdata:
            matchdata[numreplaced]=dict()
        if queryid not in matchdata[numreplaced]:
            matchdata[numreplaced][queryid]=dict()
        matchdata[numreplaced][queryid][matchid]=score
    
    for queryid in queryids:
        print queryid
        for numreplaced in sorted(matchdata):
            if queryid in matchdata[numreplaced]:
                
                bestmatch=max(matchdata[numreplaced][queryid].iteritems(), key=operator.itemgetter(1))[0]
                bestmatchscore=max(matchdata[numreplaced][queryid].iteritems(), key=operator.itemgetter(1))[1]
                outfile.write(queryid+"\t"+str(numreplaced)+"\t"+bestmatch+"\t"+str(bestmatchscore)+"\n")
    infile.close()
    outfile.close()


def compute_noise_decay(filename,quartile):
    print filename
    infile=open(filename)
    outfile=open(filename.replace("_Noise_","_NoiseDist"+str(quartile)+"_").replace("/Noise/","/Noise/Distributions/"),'w')


    matchdata=dict()
    noisedist=[]
    for line in infile:
        data=line.strip().split("\t")
        if len(data)==2:
            scorelist1,scorelist2=[float(x) for x in data[0].split(",")],[float(x) for x in data[1].split(",")]
            
            score1 = np.percentile(scorelist1, quartile)
            score2 = np.percentile(scorelist2, quartile)
            score=np.mean([score1,score2])
            noisedist.append(score)
        else:
            scorelist=[float(x) for x in data[0].split(",")]
            score = np.percentile(scorelist, quartile)
            noisedist.append(score)
    json.dump(noisedist,outfile)    
    infile.close()
    outfile.close()



def error(scorelist):
    return 2*(np.std(scorelist)/math.sqrt(len(scorelist)))


def load_results(infile,quartile,scores,metric,granularity):
    infile.next()
    for line in infile:
        queryid,numreplaced,match,score=line.strip().split()
        numreplaced=int(numreplaced)
        if metric not in scores:
            scores[metric]=dict()
        if quartile not in scores[metric]:
            scores[metric][quartile]=dict()
        if granularity not in scores[metric][quartile]:
            scores[metric][quartile][granularity]=dict()
        if numreplaced not in scores[metric][quartile][granularity]:
            scores[metric][quartile][granularity][numreplaced]=[]
        scores[metric][quartile][granularity][numreplaced].append(float(score))

    infile.close()
    return scores



def main(): 
    
    filename=sys.argv[1]
    quartile=int(sys.argv[2])
    noiseflag=sys.argv[3]

    
    

   
    if noiseflag=='Noise':
        compute_noise_decay(filename,quartile)
    else:
        compute_decay(filename,quartile)
    


if __name__ == "__main__":
    import os
    import matplotlib.lines as mlines
    import matplotlib.pyplot as plt

    import json
    import pandas as pd
    from copy import deepcopy
    from operator import itemgetter
    import operator
    from scipy.stats import rankdata
    import numpy as np
    import math
    from scipy import stats
    import sys
    main()    