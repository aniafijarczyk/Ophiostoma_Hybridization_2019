'''
The script reads ouput of calcDxy.R script (ngsTools) with dxy estimated per site.
As a second input it reads ANGSD output with Fst estimated in windows.
For a given size of windows it sums dxy for positions within a window
and divides by the number of positions in the window estimated by ANGSD in Fst output.
If number of sites is less than 50% of a window length, window is not printed.
The script will read Fst file only if it was calculated for the same window size starting at position 1.


USAGE:
    
python combineDxy.py
requires:   <ame1_ame2>.Dxy_persite.txt
            <ame1.ame2>.fst.window
output:     <ame1.ame2>_window_Dxy.txt
'''

chroms = {'OphioH327chr_1':6937932,'OphioH327chr_2':6817711,'OphioH327chr_3':3669772,'OphioH327chr_4':3419703,\
          'OphioH327chr_5':2848703,'OphioH327chr_6':2801594,'OphioH327chr_7':2758224,'OphioH327chr_8':2531247}
#chroms = {'OphioH327chr_1':43100,'OphioH327chr_2':43100}

winSize=50000
stepSize=50000


def readTab(filename):
    fh = open(filename,'r')
    linie = fh.readlines()
    k = [ele.split() for ele in linie][1:]
    D = {a[0]+"."+a[1]:float(a[-1]) for a in k}
    return D

def readFst(fname):
    fh =open(fname,'r')
    linie = fh.readlines()
    k = [ele.split() for ele in linie][1:]
    D = {ele[1]+"."+str(int(ele[2])-1) : ele[3] for ele in k}
    return D

def readGlobalDxy(fname):
    fh = open(fname,'r')
    linie = fh.readlines()
    k = [ele.split() for ele in linie]
    D = {a:b for a,b in k}
    return D



#for para in ["ulm_nov","ulm_ame1","ulm_ame2","nov_ame1","nov_ame2","ame1_ame2"]:
for para in ["ame1_ame2"]:
    pop1=para.split("_")[0]
    pop2=para.split("_")[1]

    DxyResultsFile = pop1+"_"+pop2+".Dxy_persite.txt"
    fstFile = pop1+"."+pop2+".fst.window"

    dxyDict = readTab(DxyResultsFile)
    fstDict = readFst(fstFile)
    wh = open(pop1+"."+pop2+"_window_Dxy.txt",'w')
    D = {}
    i=1
    for chr in sorted(chroms.keys()):
        x=1
        y=winSize
        step = stepSize
        regionsList = {}
        while y < chroms[chr]:
            region_start = x
            region_stop = y
            win = range(x,y+1)
            midPos = int(x-1+step/2)
            Dxy = []
            if chr+"."+str(midPos) in fstDict.keys():
                nsites = int(fstDict[chr+"."+str(midPos)])
            else:
                x+=step
                y+=step
                i+=1
                continue
            for pos in win:
                #try:
                #print(chr+"."+str(pos))
                if chr+"."+str(pos) in dxyDict.keys():
                    dxy = dxyDict[chr+"."+str(pos)]
                    Dxy.append(dxy)
                #except KeyError:
                else:
                    continue

            regionsList[i]= x, y, midPos, sum(Dxy), nsites

            x+=step
            y+=step
            i+=1
            print(chr+"\twindow:"+str(i))
            if nsites < winSize/2.:
                continue
            if len(Dxy):
                wh.write(chr+"\t"+str(midPos)+"\t"+str(midPos)+"\t"+str(sum(Dxy)/float(nsites))+"\n")
        D[chr] = regionsList

    wh.flush()
    wh.close()


