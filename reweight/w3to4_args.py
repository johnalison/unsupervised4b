import os 
import time
start = time.time()

calc3to4Path = '/uscms/home/smurthy/nobackup/CMSSW_11_1_0_pre5/src/unsupervised4b/reweight/calculate3to4WeightsIn4b.py'
yearArr = ['2016','2017','2018']
vx = [i for i in range(0,15)]

for vn in vx:
    for year in (yearArr):    
            os.system('python3 ' + calc3to4Path + ' -v ' + str(vn) + ' -y ' + year) 

end = time.time()
print(end-start)