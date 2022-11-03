import numpy as np; import uproot; import matplotlib.pyplot as plt;

class file:
    def __init__(self, name, type = None):
        self.fileName = name
        self.type = type
        self.tree = uproot.open(self.fileName)['Events']
        self.m4j = self.get_branch('m4j')
        print(self.fileName)
        print("no of events (preselection):", len(self.m4j))
        self.SB = self.get_branch('SB')
        self.SR = self.get_branch('SR')
        self.fourTag = self.get_branch('fourTag')

        self.HLT = self.get_branch('passHLT')
        self.weights = self.get_branch('mcPseudoTagWeight_3bDvTMix4bDvT_v0')
        self.trigWeight_Data = self.get_branch('trigWeight_Data')

    def get_branch(self, branch_name):
        branch = []
        branch.extend(np.ndarray.tolist((self.tree).arrays([branch_name], namedecode='utf-8')[branch_name]))
        return branch

    # selection = '(df.SB|df.SR) & df.passHLT & ~(df.SR & df.fourTag)' for data
    # selection = '(df.SB|df.SR) & df.passHLT & (df.trigWeight_Data!=0)' for MC

    def getTTSelections(self):
        self.selHLT()
        self.selFourTag()
        self.selSB()
        self.selSR()
        self.selSBorSR()
        self.selSRandFourTag()
        self.selnotSRandFourTag()
        self.selSBandHLT()
        self.selSRandHLT()
        self.selTTAll()

    def getSelections(self):
        self.selHLT()
        self.selFourTag()
        self.selSB()
        self.selSR()
        self.selSBorSR()
        self.selSRandFourTag()
        self.selnotSRandFourTag()
        self.selAll()

    def selSB(self):
        print(len(self.m4j), '->',  len(np.array(self.m4j)[self.SB]), "for sel = SB")

    def selSR(self):
        print(len(self.m4j), '->', len(np.array(self.m4j)[self.SR]), "for sel = SR")

    def selSBorSR(self, getback=False):
        if getback: return np.logical_or(self.SB, self.SR); 
        print(len(self.m4j), '->', len(np.array(self.m4j)[np.logical_or(self.SB, self.SR)]), "for sel = SB|SR")

    def selHLT(self):
        print(len(self.m4j), '->',  len(np.array(self.m4j)[self.HLT]), "for sel = HLT")

    def selFourTag(self):
        print(len(self.m4j), '->',  len(np.array(self.m4j)[self.fourTag]), "for sel = fourTag")

    def selSRandFourTag(self):
        print(len(self.m4j), '->',  len(np.array(self.m4j)[np.logical_and(self.SR,self.fourTag)]), "for sel = SR&fourTag")

    def selnotSRandFourTag(self, getback = False):
        temp = np.logical_not(np.logical_and(self.SR,self.fourTag))
        if getback: return temp; 
        print(len(self.m4j), '->',  len(np.array(self.m4j)[temp]), "for sel = ~(SR&fourTag)")

    def selNonZeroTrigWeight(self, getback = False):
        temp = np.array(self.trigWeight_Data) != 0
        print(len(self.m4j), '->',  len(np.array(self.m4j)[temp]), "for sel = trigWeight_Data != 0")
        if getback: return temp; 

    def selSRandHLT(self):
        print(len(self.m4j), '->',  len(np.array(self.m4j)[np.logical_and(self.SR,self.HLT)]), "for sel = SR&HLT")

    def selSBandHLT(self):
        print(len(self.m4j), '->',  len(np.array(self.m4j)[np.logical_and(self.SB,self.HLT)]), "for sel = SB&HLT")

    def selAll(self):
        temp = np.logical_and(self.selSBorSR(True), self.HLT)
        temp2 = np.logical_and(temp, self.selnotSRandFourTag(True))
        print(len(self.m4j), '->',  len(np.array(self.m4j)[temp]), "for sel = (SB|SR) & HLT")
        print(len(self.m4j), '->',  len(np.array(self.m4j)[temp2]), "for sel = (SB|SR) & HLT & ~(SR&fourTag)")
        print("weights =", sum(np.array(self.weights)[temp2]))

    def selTTAll(self, getback=False):
        temp = np.logical_and(self.selSBorSR(True), self.HLT)
        temp2 = np.logical_and(temp, self.selNonZeroTrigWeight(True))
        if getback:
            return len(np.array(self.m4j)[temp2]), sum(np.array(self.weights)[temp2])
        print(len(self.m4j), '->',  len(np.array(self.m4j)[temp]), "for sel = (SB|SR) & HLT")
        print(len(self.m4j), '->',  len(np.array(self.m4j)[temp2]), "for sel = (SB|SR) & HLT & (weights!=0)")
        print("weights =", sum(np.array(self.weights)[temp2]))
  
