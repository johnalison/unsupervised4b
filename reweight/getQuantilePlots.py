import numpy as np
import wquantiles as wq
import matplotlib.pyplot as plt
import matplotlib.colors as cls
from matplotlib.backends.backend_pdf import PdfPages
import uproot

def getUprootBranch(filename, vn = 0, closure = None):
    mcPsedoTagWeight = "mcPseudoTagWeight_3bDvTMix4bDvT_v" + str(vn)
    arrayNames = [mcPsedoTagWeight,"passHLT","leadStM","sublStM","m4j"]
    if closure is not None:
        arrayNames = arrayNames.append(closure)
        # arrayNames = [mcPsedoTagWeight,"nSelJets","m4j","passHLT","leadStM","sublStM"]
    if type(filename) is str: filename+=':Events'
    else: filename = [fn+":Events" for fn in filename]
    data = uproot.concatenate(filename, filter_name=arrayNames)
    return data

def get3to4Reweight(filename):
    if type(filename) is str: filename+=':Events'
    else: filename = [fn+":Events" for fn in filename]
    return np.asarray(uproot.concatenate(filename, filter_name='w3to4'))['w3to4']

def getDtoMReweight(filename):
    if type(filename) is str: filename+=':Events'
    else: filename = [fn+":Events" for fn in filename]
    return np.asarray(uproot.concatenate(filename, filter_name='wDtoM'))['wDtoM']

def getMCPTweight(ttData, vn=0):
    mcPsedoTagWeight = "mcPseudoTagWeight_3bDvTMix4bDvT_v" + str(vn)
    return np.array(ttData[mcPsedoTagWeight])

def getDataForQuantPlot(filename, dataFilename, w3to4Filename, ttFilename = None, vn = 0, mcpt = False, wDtoMFilename = None):
    data3b = getUprootBranch(filename, vn = vn)     # get relevant branches
    bgW = get3to4Reweight(w3to4Filename)  
    bgW[np.where(np.isnan(bgW)==True)] = 0
    if mcpt:
        bgW = getMCPTweight(data3b, vn = vn)
    if wDtoMFilename is not None:
        bgW = getDtoMReweight(wDtoMFilename)
    
    bgW[np.where(np.isnan(bgW)==True)] = 0
    bgW[np.where(np.isinf(bgW)==True)] = 0
    dataMixed = getUprootBranch(dataFilename, vn = vn) # get relevant branches
    dataW = getMCPTweight(dataMixed, vn = vn)
    if ttFilename is not None:
        dataTT = getUprootBranch(ttFilename, vn = vn)   # get relevant branches
        data3b = np.append(data3b, dataTT)  # append data3b with ttBar data for all branches
        bgW = np.append(bgW, getMCPTweight(dataTT, vn=vn)) # append bgW with ttBar-MCPtag weights

    return data3b, bgW, dataMixed, dataW

# wDtoM in the doesn't have mcp multiplied
def getDataForClosure(filename, wDtoMFilename, ttFilename, vn = 0, vars = None):
    print(filename)
    data = getUprootBranch(filename, vn = vn, closure=vars)     # get relevant branches
    dataDtoM = getDtoMReweight(wDtoMFilename)
    dataDtoM[np.where(np.isnan(dataDtoM)==True)] = 0
    dataDtoM[np.where(np.isinf(dataDtoM)==True)] = 0
    dataTT = getUprootBranch(ttFilename, vn = vn, closure=vars)
    return data, dataTT, dataDtoM

def getDataForPlot(filename, dataFilename, w3to4Filename, ttFilename = None, vn = 0, mcpt = False, wDtoMFilename = None):
    data3b = getUprootBranch(filename, vn = vn)     # get relevant branches
    bgW = get3to4Reweight(w3to4Filename)
    
    if mcpt:
        bgW = getMCPTweight(data3b, vn = vn)
    if wDtoMFilename is not None:
        bgW = getDtoMReweight(wDtoMFilename)
        
    bgW[np.where(np.isnan(bgW)==True)] = 0
    bgW[np.where(np.isinf(bgW)==True)] = 0
    dataMixed = getUprootBranch(dataFilename, vn = vn)  # get relevant branches
    dataW = getMCPTweight(dataMixed, vn = vn)

    dataTT, ttW = None, None
    if ttFilename is not None:
        dataTT = getUprootBranch(ttFilename, vn = vn)
        ttW = getMCPTweight(dataTT, vn=vn)
    return data3b, bgW, dataMixed, dataW, dataTT, ttW

def m4jPlot(m4jBinEdges, data3b, bgW, dataMixed, dataW, dataTT=None, ttW=None, figName="test.png", plotAll=False):
    datLabel = 'Multijet Model'
    # datLabel = '3b'
    reweightSel = [m4jBinEdges[1][0],m4jBinEdges[-2][1]]
    if plotAll: reweightSel = [0,1200]
    m4jPass = getQuantile.dataSel(data3b, reweightSel)
    m4j3b = data3b['m4j'][m4jPass]; w3b = bgW[m4jPass]
    mixPass = getQuantile.dataSel(dataMixed, [0,1200])
    m4j = dataMixed['m4j'][mixPass]; wData = dataW[mixPass]
    
    bins = np.arange(0,1200,25)
    bincount, binbound = np.histogram(m4j, bins=bins, weights=wData)
    bincenter = 0.5*(binbound[:-1] + binbound[1:])
    
    plt.rcParams.update({'font.size': 22})
    figf, axs = plt.subplots(2, 1, sharex = True, figsize=(12.5,10), gridspec_kw={'height_ratios': [0.8,0.2]})
    plt.subplots_adjust(hspace=0.1)
    if dataTT is not None:    
        TTm4jPass = getQuantile.dataSel(dataTT, reweightSel)
        m4jTT = dataTT['m4j'][TTm4jPass]; wTT = ttW[TTm4jPass]
        m4jHist = axs[0].hist([m4jTT,m4j3b], bins=bins,weights=[wTT, w3b], label=['$\mathrm{t\overline{t}}$', datLabel], color=['tab:blue', 'gold'], stacked=True, fill=True, alpha = 0.9, histtype='step', edgecolor='k')
        bkgd = m4jHist[0][1]
    else:
        m4jHist = axs[0].hist(m4j3b, bins=bins,weights=w3b, label=datLabel, color='gold', stacked=True, fill=True, alpha = 0.9, histtype='step', edgecolor='k')
        bkgd = m4jHist[0]
    
    axs[0].errorbar(bincenter, bincount, yerr=np.sqrt(bincount), marker='o', color='k', label='Data', ls='none')
    axs[0].errorbar(bincenter, bkgd, yerr = np.sqrt(bkgd), ls='none', color='k')
    # plt.title('Events (rewighted: '+ str(m4jBinEdges[1][0])+'-'+str(m4jBinEdges[-2][1])+' GeV)')
    plt.xlabel("m4j (GeV)")#, fontsize=15)
    axs[0].set_ylabel("Events")#, fontsize=15)
    axs[0].set_title('Events (rewighted: '+ str(m4jBinEdges[1][0])+'-'+str(m4jBinEdges[-2][1])+' GeV)')
    # axs[0].set_title('Events')

    axs[0].axvline(m4jBinEdges[1][0],color='k',linestyle='--')
    axs[0].axvline(m4jBinEdges[-2][1],color='k',linestyle='--')
    axs[0].axvline(m4jBinEdges[-1][1],color='k',linestyle='--', alpha=0.3)
    for i in range(len(m4jBinEdges)):
        axs[0].axvline(m4jBinEdges[i][0],color='k',linestyle='--', alpha=0.3)

    axs[0].minorticks_on()
    axs[0].set_ylim(bottom=0.0, top=20100)
    axs[0].legend(); 

    axs[1].plot(bincenter, bincount/bkgd, 'o', color='k')
    axs[1].minorticks_on()
    axs[1].axhline(1,color='k',linestyle='-', alpha=0.3)
    axs[1].set_ylabel("Data/Bkgd.")
    axs[1].set_ylim(bottom=0.1, top=1.8)
    
    # return figf
    plt.savefig(figName); 
    plt.show() ; plt.close()

def plotClosure(m4jBinEdges, data,  dataTT=None, dataDtoM = None, figName="test.png", vn='0', vars = ['m4j']):
    reweightSel = [0,1200]
    mcPsedoTagWeight = "mcPseudoTagWeight_3bDvTMix4bDvT_v" + str(vn)
    m4jPass = getQuantile.dataSel(data, reweightSel) # filter data
    wData = data[mcPsedoTagWeight][m4jPass]; 
    wDtoM = dataDtoM[m4jPass]*wData;      # 3b data and ttbar 
    
    for var in vars:
        print("Now plotting",  var)
        branchData = data[var][m4jPass]; 
        if var == 'm4j':          
            bins = np.arange(0,1200,25)
        if var == 'nSelJets':
            bins = np.arange(0,15,1)
        if var == 'leadStM' or var == 'sublStM':
            bins = np.arange(0,250,10)
        
        bincount, binbound = np.histogram(branchData, bins=bins, weights=wDtoM)
        bincenter = 0.5*(binbound[:-1] + binbound[1:])
        
        plt.rcParams.update({'font.size': 22})
        # figf = plt.figure(figsize=(11,8))
        figf, axs = plt.subplots(2, 1, sharex = True, figsize=(12.5,10), gridspec_kw={'height_ratios': [0.8,0.2]})
        plt.subplots_adjust(hspace=0.1)

        if dataTT is not None:    
            TTm4jPass = getQuantile.dataSel(dataTT, reweightSel)
            branchTT = dataTT[var][TTm4jPass]; wTT = dataTT[mcPsedoTagWeight][TTm4jPass]
            varHist = axs[0].hist(np.append(branchTT,branchData), bins=bins,weights=np.append(-wTT, wData), label="data - ttbar", color='gold',  fill=True, alpha = 0.9, histtype='step', edgecolor='k')
            axs[0].hist(branchData, bins=bins,weights= wData, label="data", hatch='/', fill=False,  alpha = 0.9, histtype='step', edgecolor='k')
        
        # print(varHist)
        axs[0].errorbar(bincenter, bincount, yerr=np.sqrt(bincount), marker='o', color='k', label='wDtoM', ls='none')
        # axs[0].errorbar(bincenter, varHist[0], yerr=np.sqrt(bincount), marker='o', color='r', label='check', ls='none')
        plt.xlabel(var)#"m4j (GeV)")#, fontsize=15)
        axs[0].set_ylabel("Events")#, fontsize=15)
        axs[0].set_title('Events')# (rewighted: '+ str(m4jBinEdges[1][0])+'-'+str(m4jBinEdges[-2][1])+' GeV)')
        axs[0].minorticks_on()
        axs[0].legend(); 

        # axs[1].plot(bincenter, 100*np.divide(bincount-varHist[0],varHist[0], where=bincount!=0), 'o', color='k')
        # axs[1].minorticks_on()
        # axs[1].axhline(0,color='k',linestyle='-', alpha=0.3)
        # axs[1].set_ylabel("percent error")
        # # axs[1].set_ylim(bottom=0.1, top=1.8)

        axs[1].plot(bincenter, np.divide(bincount,varHist[0], where=varHist[0]!=0), 'o', color='k')
        axs[1].minorticks_on()
        axs[1].axhline(1,color='k',linestyle='-', alpha=0.3)
        axs[1].set_ylabel("wDtoM/(data-ttbar)")
        axs[1].set_ylim(bottom=0.1)#, top=1.8)
        # axs[1].set_yscale("log")

        plt.savefig(figName+"_"+var+".png"); 
        # plt.show() ; 
        plt.close()


def plotClosureDijet(m4jBinEdges, data,  dataTT=None, dataDtoM = None, figName="test.png", vn='0', vars = ['m4j']):
    reweightSel = [0,1200]
    mcPsedoTagWeight = "mcPseudoTagWeight_3bDvTMix4bDvT_v" + str(vn)
    m4jPass = getQuantile.dataSel(data, reweightSel) # filter data
    wData = data[mcPsedoTagWeight][m4jPass]; 
    wDtoM = dataDtoM[m4jPass]*wData;      # 3b data and ttbar 
    leadData = data["leadStM"][m4jPass]; 
    sublData = data["sublStM"][m4jPass]; 

    TTm4jPass = getQuantile.dataSel(dataTT, reweightSel)
    leadTT = dataTT["leadStM"][TTm4jPass]; sublTT = dataTT["sublStM"][TTm4jPass];
    wTT = dataTT[mcPsedoTagWeight][TTm4jPass]

    lead = np.append(leadData, leadTT)
    subl = np.append(sublData, sublTT)
    weight = np.append(wData, wTT*-1)

    plt.rcParams.update({'font.size': 22})
    figf, axs = plt.subplots(1, 2, figsize=(10.5,4.2), gridspec_kw={'width_ratios': [1/2,1/2]})
    # plt.subplots_adjust(hspace=0.1)
    bincont0,xedge,yedge,im = axs[0].hist2d(lead, subl, range=[[0,250],[0,250]],bins = 25, weights=weight);
    bincont1,xedge,yedge,im = axs[1].hist2d(leadData, sublData,range=[[0,250],[0,250]],bins = 25, weights=wDtoM);
    # bincont2 = np.divide(bincont1, bincont0, where=bincont0!=0)
    # axs[2].imshow(bincont2, origin='lower',extent=np.array([0,250,0,250])); 
    
    axs[0].set_xlim(left=0, right=250); axs[0].set_ylim(bottom=0, top=250)
    axs[1].set_xlim(left=0, right=250); axs[1].set_ylim(bottom=0, top=250)

    norm = cls.Normalize(vmin=np.min([np.min(bincont0), np.min(bincont1)]), vmax=np.max([np.max(bincont0), np.max(bincont1)]))
    sm = plt.cm.ScalarMappable(norm=norm)

    for i in range(2):
        # axs[i].colorbar();
        axs[i].set_xlabel("leadStM (GeV)");
        axs[i].set_ylabel("sublStM (GeV)");
        plt.colorbar(sm, ax= axs[i])
    
    axs[0].set_title("data - ttbar");
    axs[1].set_title("wDtoM");
    # axs[2].set_title("wDtoM/(data-ttbar)");
        
    plt.savefig(figName+"_dijet.png"); 
    plt.show() ; 
    plt.close()


def make_color_plot_single(content, bins, figf, axs, i=0, title = "m4j", xlabel = "leadStM (GeV)", ylabel = "sublStM (GeV)", rig = 0):      
        cmap = plt.cm.jet
        plt.rc('font', size=25)
    #     cmap = plt.cm.cividis
    #     cmap = plt.cm.viridis
        norm = cls.Normalize(vmin=min(content)-rig, vmax=max(content)+rig)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        for binN, bound in enumerate(bins):
            axs[i].fill_between([bound[0,0],bound[0,1]],bound[1,0],bound[1,1], color=cmap(norm(content[binN])))
        if i==0:
            for binN, bound in enumerate(bins):
                axs[i].hlines(bound[1,0],bound[0,0],bound[0,1], linewidth = .2, color = 'k')
                axs[i].hlines(bound[1,1],bound[0,0],bound[0,1], linewidth = .2, color = 'k')
                axs[i].vlines(bound[0,0],bound[1,0],bound[1,1], linewidth = .2, color = 'k')
                axs[i].vlines(bound[0,1],bound[1,0],bound[1,1], linewidth = .2, color = 'k')
        axs[i].set_title(title)
        axs[i].set_xlabel(xlabel);axs[i].set_ylabel(ylabel)
        axs[i].set_xlim(left=0, right=250); axs[i].set_ylim(bottom=0, top=250)
        plt.colorbar(sm, ax= axs[i])

def plotClosureDijetQuant(m4jBinEdges, data,  dataTT=None, dataDtoM = None, figName="test.png", vn='0', vars = ['m4j']):
    reweightSel = [0,1200]
    mcPsedoTagWeight = "mcPseudoTagWeight_3bDvTMix4bDvT_v" + str(vn)
    m4jPass = getQuantile.dataSel(data, reweightSel) # filter data
    wData = data[mcPsedoTagWeight][m4jPass]; 
    wDtoM = dataDtoM[m4jPass]*wData;      # 3b data and ttbar 
    leadData = data["leadStM"][m4jPass]; 
    sublData = data["sublStM"][m4jPass]; 

    TTm4jPass = getQuantile.dataSel(dataTT, reweightSel)
    leadTT = dataTT["leadStM"][TTm4jPass]; sublTT = dataTT["sublStM"][TTm4jPass];
    wTT = dataTT[mcPsedoTagWeight][TTm4jPass]

    lead = np.append(leadData, leadTT)
    subl = np.append(sublData, sublTT)
    weightCombined = np.append(wData, wTT*-1)

    priData = np.array([np.array(lead), np.asarray(subl)]).transpose()
    bins = getQuantile.recursive_split(priData, bins =[], threshold=60, weights=weightCombined)

    priCount = []
    for binbound in bins:
        priCount.append(getQuantile.get_bin_content(priData,binbound,weightCombined))

    secCount = []
    secData = np.array([np.array(leadData), np.asarray(sublData)]).transpose()    
    for binbound in bins:
        secCount.append(getQuantile.get_bin_content(secData,binbound,wDtoM))
    
    plt.rcParams.update({'font.size': 22})
    # figf, axs = plt.subplots(1, 2, figsize=(10.2,4.5), gridspec_kw={'width_ratios': [1/2,1/2]}, constrained_layout=True)
    figf, axs = plt.subplots(1, 3, figsize=(13.2,4.5), gridspec_kw={'width_ratios': [1/3,1/3, 1/3]}, constrained_layout=True)
    make_color_plot_single(priCount,  bins, figf, axs,i=0, title = "data - ttbar", rig=1)
    make_color_plot_single(secCount, bins, figf, axs,i=1, title = "WDtoM", rig=1)
    ratioCount = np.divide(secCount, priCount, where=priCount!=0)
    make_color_plot_single(ratioCount, bins, figf, axs,i=2, title = "WDtoM/(data-ttbar)", rig=1)

    plt.savefig(figName+"_dijet_quant.png"); 
    plt.show() ; 
    plt.close()



def plotClosureDijetOld(m4jBinEdges, data,  dataTT=None, dataDtoM = None, figName="test.png", vn='0', vars = ['m4j']):
    reweightSel = [0,1200]
    mcPsedoTagWeight = "mcPseudoTagWeight_3bDvTMix4bDvT_v" + str(vn)
    m4jPass = getQuantile.dataSel(data, reweightSel) # filter data
    wData = data[mcPsedoTagWeight][m4jPass]; 
    wDtoM = dataDtoM[m4jPass]*wData;      # 3b data and ttbar 
    leadData = data["leadStM"][m4jPass]; 
    sublData = data["sublStM"][m4jPass]; 

    TTm4jPass = getQuantile.dataSel(dataTT, reweightSel)
    leadTT = dataTT["leadStM"][TTm4jPass]; sublTT = dataTT["sublStM"][TTm4jPass];
    wTT = dataTT[mcPsedoTagWeight][TTm4jPass]

    lead = np.append(leadData, leadTT)
    subl = np.append(sublData, sublTT)
    weight = np.append(wData, wTT*-1)

    plt.rcParams.update({'font.size': 22})
    figf, axs = plt.subplots(1, 3, figsize=(12.5,4.2), gridspec_kw={'width_ratios': [1/3,1/3,1/3]})
    # plt.subplots_adjust(hspace=0.1)
    bincont0,xedge,yedge,im = axs[0].hist2d(lead, subl, range=[[0,250],[0,250]],bins = 25, weights=weight);
    bincont1,xedge,yedge,im = axs[1].hist2d(leadData, sublData,range=[[0,250],[0,250]],bins = 25, weights=wDtoM);
    bincont2 = np.divide(bincont1, bincont0, where=bincont0!=0)
    axs[2].imshow(bincont2, origin='lower',extent=np.array([0,250,0,250])); 
    
    axs[0].set_xlim(left=0, right=250); axs[0].set_ylim(bottom=0, top=250)
    axs[1].set_xlim(left=0, right=250); axs[1].set_ylim(bottom=0, top=250)

    norm = cls.Normalize(vmin=np.min([np.min(bincont0), np.min(bincont1), np.min(bincont2)]), vmax=np.max([np.max(bincont0), np.max(bincont1), np.max(bincont2)]))
    sm = plt.cm.ScalarMappable(norm=norm)

    for i in range(3):
        # axs[i].colorbar();
        axs[i].set_xlabel("leadStM (GeV)");
        axs[i].set_ylabel("sublStM (GeV)");
        plt.colorbar(sm, ax= axs[i])
    
    axs[0].set_title("data - ttbar");
    axs[1].set_title("wDtoM");
    axs[2].set_title("wDtoM/(data-ttbar)");
        
    plt.savefig(figName+"_dijet.png"); 
    plt.show() ; 
    plt.close()


def m4jPlotClosureOld(m4jBinEdges, data,  dataTT=None, dataDtoM = None, figName="test.png", vn='0', vars = ['m4j']):
    reweightSel = [0,1200]
    mcPsedoTagWeight = "mcPseudoTagWeight_3bDvTMix4bDvT_v" + str(vn)
    m4jPass = getQuantile.dataSel(data, reweightSel) # filter data
    wData = data[mcPsedoTagWeight][m4jPass]; 
    wDtoM = dataDtoM[m4jPass]*wData;      # 3b data and ttbar 
    
    for var in vars:
        print("Now plotting",  var)
        branchData = data[var][m4jPass]; 
        if var == 'm4j':          
            bins = np.arange(0,1200,25)
        if var == 'nSelJets':
            bins = np.arange(0,15,1)
        bincount, binbound = np.histogram(branchData, bins=bins, weights=wDtoM)
        bincenter = 0.5*(binbound[:-1] + binbound[1:])
        
        plt.rcParams.update({'font.size': 22})
        figf = plt.figure(figsize=(11,8))

        if dataTT is not None:    
            TTm4jPass = getQuantile.dataSel(dataTT, reweightSel)
            branchTT = dataTT[var][TTm4jPass]; wTT = dataTT[mcPsedoTagWeight][TTm4jPass]
            m4jHist = plt.hist(np.append(branchTT,branchData), bins=bins,weights=np.append(-wTT, wData), label="data - ttbar", color='tab:blue',  fill=True, alpha = 0.9, histtype='step', edgecolor='k')
            plt.hist(branchData, bins=bins,weights= wData, label="data", hatch='/', fill=False, color='red', alpha = 0.9, histtype='step', edgecolor='k')
        
        plt.errorbar(bincenter, bincount, yerr=np.sqrt(bincount), marker='o', color='k', label='wDtoM', ls='none')
        plt.xlabel(var)#"m4j (GeV)")#, fontsize=15)
        plt.ylabel("Events")#, fontsize=15)
        plt.title(var)
        plt.title('Events')# (rewighted: '+ str(m4jBinEdges[1][0])+'-'+str(m4jBinEdges[-2][1])+' GeV)')
        plt.legend(); 
        plt.savefig(figName+"_"+var+".png"); 
        plt.show() ; plt.close()

def saveDijetPlots(m4jBinEdges,data3b, bgW, dataMixed, dataW, pdfFilename,threshold=100):
    print(pdfFilename)
    pp = PdfPages(pdfFilename)
    
    pullFull = []
    meanStdArr = []
    for binNo, m4jBinEdge in enumerate(m4jBinEdges):
        if binNo == 0 or binNo == len(m4jBinEdges)-1:
            continue
        if binNo == 3: break
        m4jPass_SR3b = getQuantile.dataSel(data3b, m4jBinEdges[binNo])
        leadStM_SR3b = data3b['leadStM'][m4jPass_SR3b]
        sublStM_SR3b = data3b['sublStM'][m4jPass_SR3b]
        w3b = bgW[m4jPass_SR3b]
        
        m4jPass_SRMixed = getQuantile.dataSel(dataMixed, m4jBinEdges[binNo])
        leadStM_mixed = dataMixed['leadStM'][m4jPass_SRMixed]
        sublStM_mixed = dataMixed['sublStM'][m4jPass_SRMixed]
        # wMixed = dataMixed['mcPseudoTagWeight_3bDvTMix4bDvT_v0'][m4jPass_SRMixed]
        wMixed = dataW[m4jPass_SRMixed]

        qData = getQuantile(leadStM_SR3b, sublStM_SR3b, threshold = threshold, weights = w3b, binNo=binNo, m4jBinEdge=m4jBinEdge)
        qData.defData(leadStM_mixed, sublStM_mixed, dataW = wMixed)
        figf = qData.getHistPlots()
        pullStat = qData.get_pull_weights(getData=True)
        meanStdArr.append(pullStat[0:2])
        pullFull.append(pullStat[2])
        plt.savefig("dijet_quals.png"); 
        pp.savefig(figf); # plt.show();
        plt.close()

    pp.close()
    np.save(pdfFilename[:-4]+'_pullFull.npy',pullFull)
    np.savetxt(pdfFilename[:-4]+"_MeanStd.txt", meanStdArr)
    print("pdf saved\n")


class getQuantile:
    def __init__(self, bgDataX, bgDataY, threshold, binNo=1, m4jBinEdge=[0,1200], weights = None):
        self.bkg = np.array([np.array(bgDataX), np.array(bgDataY)]).transpose()
        self.weights = weights
        self.bins = getQuantile.recursive_split(self.bkg, bins =[], threshold=threshold, weights=self.weights)
        self.binNo = binNo
        self.m4jBinEdge = m4jBinEdge
        self.bkgBinCont = self.getQuantBinContent(bgDataX, bgDataY, weights=self.weights)
        self.data = None
        self.dataBinCont = None
        self.pull_weights = np.zeros(len(self.bkg))
        
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

    @staticmethod
    def get_data_pull_val(data,binbound,pull_weights,pullVal):
        mask  = (data[:,0]>binbound[0,0]) # x > xstart
        mask &= (data[:,0]<binbound[0,1]) # x < xend
        mask &= (data[:,1]>binbound[1,0]) # y > ystart
        mask &= (data[:,1]<binbound[1,1]) # y < yend
        pull_weights[mask] = pullVal
        return pull_weights

    def get_pull_weights(self, getData = False):
        pull = self.calculatePull()
        for bInd, binbound in enumerate(self.bins):
            self.pull_weights = self.get_data_pull_val(self.bkg,binbound,self.pull_weights,pull[bInd])
        if getData:
           return np.mean(self.pull_weights), np.std(self.pull_weights), self.pull_weights 
        return np.mean(self.pull_weights), np.std(self.pull_weights)

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
        return (np.array(self.dataBinCont) - np.array(self.bkgBinCont))/(2.5+np.sqrt(np.array(self.bkgBinCont)))
    
    # makes quantile plots
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
    
    # makes 2d histograms
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
        getQuantile.make_color_plot(self.bkgBinCont, self.bins, figf, axs,0,1, title = "Background, SR "+str(self.binNo), rig=1)
        getQuantile.make_color_plot(self.dataBinCont, self.bins, figf, axs,1,1, title = "Data, SR "+str(self.binNo), rig=1)
        getQuantile.make_color_plot(pull,self.bins, figf, axs,0,2, title = "Pull, SR "+str(self.binNo), rig=0)
        getQuantile.make_color_hist(self.bkg, 25, figf, axs,0,0,title = "Background", weights=self.weights)
        getQuantile.make_color_hist(self.data, 25, figf, axs,1,0,title = "Data")
        axs[1,2].hist(pull, bins = 15); axs[1, 2].set_xlim([-5.2,6.5]); axs[1, 2].set_xlabel("Pull value"); axs[1, 2].set_ylabel("No. of quantile bins"); axs[1, 2].set_title("Pull Distribution")
        axs[1,2].set_xticks(ticks=np.linspace(-4,4,9) , minor=True)
        # plt.tight_layout()
        return figf
        