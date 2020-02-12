from ROOT import *
from config import config
import time

def getKey(phi,eta): return str( 100*phi + eta )
def getPhiEta(key): return int(key)/100,int(key)%100
def deltaR(iphi1,ieta1,iphi2,ieta2):
    return iphi2-iphi1,ieta2-ieta1
def irange(lo,hi,step=1): return range(lo,hi+1,step)
def BuildCalo(calo):
    calomap = {}
    for ieta in irange(1,calo.GetNbinsY()):
        for iphi in irange(1,calo.GetNbinsX()):
            key = getKey(iphi,ieta)
            calomap[key] = Tower(iphi,ieta,calo.GetBinContent(iphi,ieta))
    return calomap
def keyiter(keylist):
    for key in keylist: yield getPhiEta(key)
    
    
class Tower:
    def __init__(self,iphi,ieta,et):
        self.iphi = iphi
        self.ieta = ieta
        self.et = et
    def __str__(self):
        return 'iPhi: %i iEta: %i Et: %f' % (self.iphi,self.ieta,self.et)
    def __eq__(self,other):
        return self.iphi == other.iphi and self.ieta == other.ieta and self.et == other.et
    def __gt__(self,other):
        if other is None: return True
        dphi,deta = deltaR(self.iphi,self.ieta,other.iphi,other.ieta)
        if self.et == other.et:
            if deta == 0 and dphi == 0: return False
            if deta >= 0 and dphi >= 0: return False
            if deta <= 0 and dphi <= 0: return True
            if deta >  0 and dphi <  0: return abs(dphi) >= abs(deta)
            if deta <  0 and dphi >  0: return abs(dphi) < abs(deta)
        return self.et > other.et
class Jet:
    def __init__(self):
        self.iphi = -1
        self.ieta = -1
        self.et = -1
        self.towers = []
    def seed(self,seed):
        if config.debug: print '---Seeding Jet %s' % str(seed)
        self.iphi = seed.iphi
        self.ieta = seed.ieta
        self.seed_et = seed.et
        self.et = seed.et
        self.towers.append(seed)
    def add(self,tower):
        if config.debug and tower.et > 0: print '----Adding Tower %s' % str(tower)
        self.et += tower.et
        self.towers.append(tower)
    def __str__(self):
        return 'iPhi: %i iEta: %i Et: %f' % (self.iphi,self.ieta,self.et)

class Info:
    def __init__(self,fname):
        self.fname = fname
        self.tfile = TFile.Open(fname)
        self.nevents = sum([ key.GetName() == 'seed' for key in self.tfile.GetListOfKeys() ])
    def getCalo(self,i):
        if i < self.nevents:
            calo = self.tfile.Get('calo;%i'%(i+1))
            return BuildCalo(calo)
    def __len__(self): return self.nevents
