import numpy as np
import wquantiles as wq
import matplotlib.pyplot as plt
import matplotlib.colors as cls
from matplotlib.backends.backend_pdf import PdfPages
import uproot

def getUprootBranch(filename, vn = 0):
    f = uproot.open(filename)['Events']
    mcPsedoTagWeight = "mcPseudoTagWeight_3bDvTMix4bDvT_v" + str(vn)
    arrayNames = [mcPsedoTagWeight,"passHLT","leadStM","sublStM","m4j"]
    data = f.arrays(arrayNames)
    return data

def get3to4Reweight(filename):
    return np.array((uproot.open(filename)['Events']).arrays('w3to4'))['w3to4']

def getBGWeight(w3to4, ttData, vn=0):
    mcPsedoTagWeight = "mcPseudoTagWeight_3bDvTMix4bDvT_v" + str(vn)
    return np.append(w3to4, np.array(ttData[mcPsedoTagWeight]))

def getDataForPlot(filename, dataFilename, w3to4Filename, ttFilename = None, vn = 0):
    data3b = getUprootBranch(filename, vn = vn)
    dataMixed = getUprootBranch(dataFilename, vn = vn)
    bgW = get3to4Reweight(w3to4Filename)
    bgW[np.where(np.isnan(bgW)==True)] = 0
    # print(len(bgW), len(dataMixed))
    if ttFilename is not None:
        for ttFile in ttFilename:
            tt = getUprootBranch(ttFile, vn = vn)
            data3b = np.append(data3b, tt)
            bgW = getBGWeight(bgW, tt, vn=vn)
        
    # print(len(bgW), len(dataMixed))
    return data3b, bgW, dataMixed

def m4jPlot(m4jBinEdges, data3b, bgW, dataMixed, figName):
    m4jPass = getQuantile.dataSel(data3b, [m4jBinEdges[1][0],m4jBinEdges[-2][1]])
    m4j3b = data3b['m4j'][m4jPass]
    w = bgW[m4jPass]
    mixPass = getQuantile.dataSel(dataMixed, [0,1200])
    m4j = dataMixed['m4j'][mixPass]

    plt.figure(figsize=(10,5))
    plt.hist(m4j,bins=np.arange(0,1200,25), label="'Mixed' data weights", color = 'C0', alpha = 0.7)
    plt.hist(m4j3b,bins=np.arange(0,1200,25),weights=w, label="Background model weights", color = 'C1', alpha = 0.7)
    plt.title('Event Weights (rewighted: '+ str(m4jBinEdges[1][0])+'-'+str(m4jBinEdges[-2][1])+' GeV)')
    plt.xlabel("m4j (GeV)")#, fontsize=15)
    plt.ylabel("Event weights")#, fontsize=15)
    plt.vlines(m4jBinEdges[1][0],0,5000,color='k',linestyle='--')
    plt.vlines(m4jBinEdges[-2][1],0,5000,color='k',linestyle='--')
    plt.vlines(m4jBinEdges[-1][1],0,5000,color='k',linestyle='--', alpha=0.3)
    for i in range(len(m4jBinEdges)):
        plt.vlines(m4jBinEdges[i][0],0,5000,color='k',linestyle='--', alpha=0.3)

    plt.legend(); 
    plt.savefig(figName); plt.show() ; plt.close()



class getQuantile:
    def __init__(self, bgDataX, bgDataY, threshold, binNo=1, m4jBinEdge=[0,1200], weights = None):
        self.bkg = np.array([np.asarray(bgDataX), np.asarray(bgDataY)]).transpose()
        self.weights = weights
        self.bins = getQuantile.recursive_split(self.bkg, bins =[], threshold=threshold, weights=self.weights)
        self.binNo = binNo
        self.m4jBinEdge = m4jBinEdge
        self.bkgBinCont = self.getQuantBinContent(bgDataX, bgDataY, weights=self.weights)
        self.data = None
        self.dataBinCont = None
        
    def getBins(self):
        return self.bins
    
    def defData(self, dataX, dataY, dataW = None):
        self.data = np.array([np.asarray(dataX), np.asarray(dataY)]).transpose()
        self.dataBinCont = self.getQuantBinContent(dataX, dataY, weights = dataW )
        
    @staticmethod    
    def dataSel(data, binEdge, binEdgeSec = False):
        sel = np.asarray(np.where((data['m4j']>binEdge[0]) & (data['m4j']<=binEdge[1]), True, False))
        if binEdgeSec is not False:
            sel |= np.asarray(np.where((data['m4j']>binEdgeSec[0]) & (data['m4j']<=binEdgeSec[1]), True, False))
        sel &= np.asarray(np.where((data['leadStM']>0) & (data['leadStM']<=250), True,False))
        sel &= np.asarray(np.where((data['sublStM']>0) & (data['sublStM']<=250), True,False))
        sel &= np.asarray(np.where(data['passHLT'], True, False))
        return sel

    @staticmethod
    def get_median(arr, weights = None):
        if weights is None: return np.median(arr)
        else: return wq.median(arr, weights)

    @staticmethod
    def split(a, axis=0, box=np.array([[0.,250.],[0.,250.]]), weights = None):
        s = getQuantile.get_median(a[:,axis], weights) # finds the 'a' value that occurs at 50%
        # make boundary array with x,y points to define the line 
        boundary = np.array([[0.,0.],[0.,0.]]) # [[xstart,xend],[ystart,yend]]
        boundary[axis] = np.array([s,s]) # the 50% value along that axis for start & end
        other_axis = (axis+1)%2 # returns 0 if 1, and 1 if 0
        boundary[other_axis] = box[other_axis] # set default start & end values - basically full range of that axis from 'box'
        mask = a[:,axis]<s # returns an array of boolean based on each data point that lies below the cut
        lower, upper = a[mask], a[~mask] # get lower and upper halves of events in this plane
        if weights is not None:
            lw, uw = weights[mask], weights[~mask]
        else: lw, uw = None, None
        return lower, boundary, upper, lw, uw

    @staticmethod
    def recursive_split(a, bins, threshold, axis=0, box=np.array([[0.,250.],[0.,250.]]), weights = None):
        if (weights is not None and np.sum(weights)>threshold) or (weights is None and a.shape[0]>threshold):
            lower, boundary, upper, lw, uw = getQuantile.split(a, axis=axis, box=box, weights=weights)
            split_value = boundary[axis,axis] # returns s
            next_axis = (axis+1)%2
            box_l = box.copy()
            box_u = box.copy()
            box_l[axis, 1] = split_value # end of line at split for next divide; box_l ranges from lower bound of box up to split
            box_u[axis, 0] = split_value # start of line at split; box_u ranges from split up to upper bound of box
            
            bins = getQuantile.recursive_split(lower, bins, threshold=threshold, axis=next_axis, box=box_l, weights = lw)
            bins = getQuantile.recursive_split(upper, bins, threshold=threshold, axis=next_axis, box=box_u, weights = uw)
        else:
            # when there are events less than the threshold: bins contains [[xstart, xend], [ystart, yend]] and makes the box
            bins.append(box)
        return bins

    @staticmethod
    def get_bin_content(data,binbound,weights=None):
        mask  = (data[:,0]>binbound[0,0]) # x > xstart
        mask &= (data[:,0]<binbound[0,1]) # x < xend
        mask &= (data[:,1]>binbound[1,0]) # y > ystart
        mask &= (data[:,1]<binbound[1,1]) # y < yend
        if weights is None:
            return len(data[mask])
        else:
            return np.sum(weights[mask])

    def getQuantBinContent(self, dataX, dataY, weights=None):
        data = np.array([np.asarray(dataX), np.asarray(dataY)]).transpose()
        temp = []
        if weights is not None:
            for binbound in self.bins:
                temp.append(self.get_bin_content(data,binbound,weights))
            return temp
        for binbound in self.bins:
            temp.append(self.get_bin_content(data,binbound))
        return temp
    
    def calculatePull(self):
        return (np.array(self.dataBinCont) - np.array(self.bkgBinCont))/np.sqrt(2.5+np.array(self.bkgBinCont))
        
    @staticmethod
    def make_color_plot(content, bins, figf, axs, i=0, j=0, title = "m4j", xlabel = "lead_St M (GeV)", ylabel = "subl_St M (GeV)", rig = 0):      
        cmap = plt.cm.jet
        plt.rc('font', size=25)
    #     cmap = plt.cm.cividis
    #     cmap = plt.cm.viridis
        norm = cls.Normalize(vmin=min(content)-rig, vmax=max(content)+rig)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        for binN, bound in enumerate(bins):
            axs[i, j].fill_between([bound[0,0],bound[0,1]],bound[1,0],bound[1,1], color=cmap(norm(content[binN])))
        if i==0 and j ==1:
            for binN, bound in enumerate(bins):
                axs[i, j].hlines(bound[1,0],bound[0,0],bound[0,1], linewidth = .2, color = 'k')
                axs[i, j].hlines(bound[1,1],bound[0,0],bound[0,1], linewidth = .2, color = 'k')
                axs[i, j].vlines(bound[0,0],bound[1,0],bound[1,1], linewidth = .2, color = 'k')
                axs[i, j].vlines(bound[0,1],bound[1,0],bound[1,1], linewidth = .2, color = 'k')
        axs[i, j].set_title(title)
        axs[i, j].set_xlabel(xlabel);axs[i, j].set_ylabel(ylabel)
        axs[i, j].set_xlim(left=0, right=250); axs[i, j].set_ylim(bottom=0, top=250)
        plt.colorbar(sm, ax= axs[i, j])
    
    @staticmethod
    def make_color_hist(data, bins, figf, axs, i=0, j=0, title = "m4j", xlabel = "lead_St M (GeV)", ylabel = "subl_St M (GeV)", weights=None):      
        plt.rc('font', size=25) 
        a,b,c,d = axs[i,j].hist2d(data[:,0], data[:,1] , range=[[0,250],[0,250]],bins = bins, weights=weights);
        norm = cls.Normalize(vmin=np.amin(a), vmax=np.amax(a))
        sm = plt.cm.ScalarMappable(norm=norm)
        axs[i, j].set_title(title)
        axs[i, j].set_xlabel(xlabel);axs[i, j].set_ylabel(ylabel)
        axs[i, j].set_xlim(left=0, right=250); axs[i, j].set_ylim(bottom=0, top=250)
        plt.colorbar(sm, ax= axs[i, j])
    
    def getHistPlots(self):
        # figf, axs = plt.subplots(2, 3, figsize=(40,20.5), layout="constrained") # which version?
        plt.rc('font', size=25) 
        figf, axs = plt.subplots(2, 3, figsize=(40,20.5), constrained_layout=True)
        figf.suptitle('m4j: '+str(self.m4jBinEdge[0])+'-'+str(self.m4jBinEdge[1])+' GeV');
        # figf.suptitle('m4j: '+str(self.m4jBinEdge[0])+'-'+str(self.m4jBinEdge[1])+' GeV', y=0.96); # figf.tight_layout()
        pull = self.calculatePull()
        getQuantile.make_color_plot(self.bkgBinCont, self.bins, figf, axs,0,1, title = "Background, bin "+str(self.binNo), rig=1)
        getQuantile.make_color_plot(self.dataBinCont, self.bins, figf, axs,1,1, title = "Data, bin "+str(self.binNo), rig=1)
        getQuantile.make_color_plot(pull,self.bins, figf, axs,0,2, title = "Pull, bin "+str(self.binNo), rig=0)
        getQuantile.make_color_hist(self.bkg, 25, figf, axs,0,0,title = "Background", weights=self.weights)
        getQuantile.make_color_hist(self.data, 25, figf, axs,1,0,title = "Data")
        axs[1,2].hist(pull, bins = 15); axs[1, 2].set_xlim([-4,4]); axs[1, 2].set_xlabel("Pull value"); axs[1, 2].set_ylabel("No. of quantile bins"); axs[1, 2].set_title("Pull Distribution")
        axs[1,2].set_xticks(ticks=np.linspace(-4,4,9) , minor=True)
        # plt.tight_layout()
        return figf
        