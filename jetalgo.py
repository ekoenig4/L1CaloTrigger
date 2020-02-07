from ROOT import *
from sys import argv
def getKey(phi,eta): return str( 100*phi + eta )
def getPhiEta(key): return int(key)/100,int(key)%100
def irange(lo,hi,step=1): return range(lo,hi,step)
class Tower:
    def __init__(self,iphi,ieta,et):
        self.iphi = iphi
        self.ieta = ieta
        self.et = et
class Jet:
    def __init__(self):
        self.iphi = -1
        self.ieta = -1
        self.et = -1
        self.towers = []
    def seed(self,seed):
        self.iphi = seed.iphi
        self.ieta = seed.et
        self.seed_et = seed.et
        self.et = seed.et
        self.towers.append(seed)
    def add(self,tower):
        self.et += tower.et
def BuildCalo(calo):
    calomap = {}
    for iphi in irange(1,calo.GetNbinsX()):
        for ieta in irange(1,calo.GetNbinsY()):
            calomap[getKey(iphi,ieta)] = Tower(iphi,ieta,calo.GetBinContent(iphi,ieta))
    return calomap
class Info:
    def __init__(self,fname):
        self.fname = fname
        self.tfile = TFile.Open(fname)
        self.nevents = sum([ key.GetName() == 'seed' for key in self.tfile.GetListOfKeys() ])
    def __iter__(self):
        self.ievent = 0
        return self
    def next(self):
        if self.ievent < self.nevents:
            calo = self.tfile.Get('calo;%i'%(self.ievent+1))
            self.ievent += 1
            return BuildCalo(calo)
        return StopIteration
def getGCTTowers(center_key,towers):
    subset = {}
    cphi,ceta = getPhiEta(center_key)
    dphi,deta = cphi - 6,ceta - 1
    for ieta in irange(1,34):
        for iphi in irange(1,32):
            rphi,reta = (iphi + dphi)%72+1,(ieta + deta)%34
            rkey = getKey(rphi,reta)
            nkey = getKey(iphi,ieta)
            if rkey not in towers: continue
            subset[nkey] = towers[rkey]
    return subset
def get25x25Towers(center_key,towers):
    subset = {}
    cphi,ceta = getPhiEta(center_key)
    dphi,deta = cphi - 13,ceta - 13
    for ieta in irange(1,25):
        for iphi in irange(1,25):
            rphi,reta = (iphi + dphi),(ieta + deta)
            rkey = getKey(rphi,reta)
            nkey = getKey(iphi,ieta)
            if rkey not in towers: continue
            subset[nkey] = towers[rkey]
    return subset
def get21x21Towers(center_key,towers):
    subset = {}
    cphi,ceta = getPhiEta(center_key)
    dphi,deta = cphi - 11,ceta - 11
    for ieta in irange(1,21):
        for iphi in irange(1,21):
            rphi,reta = (iphi + dphi),(ieta + deta)
            rkey = getKey(rphi,reta)
            nkey = getKey(iphi,ieta)
            if rkey not in towers: continue
            subset[nkey] = towers[rkey]
    return subset
def get7x7Towers(center_key,towers):
    subset = {}
    cphi,ceta = getPhiEta(center_key)
    dphi,deta = cphi - 4,ceta - 4
    for ieta in irange(1,7):
        for iphi in irange(1,7):
            rphi,reta = (iphi + dphi),(ieta + deta)
            rkey = getKey(rphi,reta)
            nkey = getKey(iphi,ieta)
            if rkey not in towers: continue
            subset[nkey] = towers[rkey]
    return subset
def getMaxTowerIn4x4(center,towers):
    cphi,ceta = getPhiEta(center)
    dphi,deta = cphi - 2,ceta - 3
    max_tower = {'key':-1,'et':-1}
    for iphi in irange(1,4):
        for ieta in irange(1,4):
            rphi,reta = iphi + dphi,ieta + deta
            key = getKey(rphi,reta)
            if key not in towers: continue
            if max_tower['et'] < towers[key].et:
                max_tower['key'] = key
                max_tower['et'] = towers.et
    return max_tower
def getMaxTowerIn7x7(center,towers):
    cphi,ceta = getPhiEta(center)
    dphi,deta = cphi - 4,ceta - 4
    max_tower = {'key':-1,'et':-1}
    for iphi in irange(1,7):
        for ieta in irange(1,7):
            rphi,reta = iphi + dphi,ieta + deta
            key = getKey(rphi,reta)
            if key not in towers: continue
            if max_tower['et'] < towers[key].et:
                max_tower['key'] = key
                max_tower['et'] = towers.et
    return max_tower
def Cluster7x7(jet,towers):
    seed = towers['404']
    jet.seed(seed)
    for key,tower in towers.iteritems():
        if key == '404': continue
        jet.add(tower)
def Cluster21x21(jet,towers):
    t7 = get7x7Towers('1111',towers)
    Cluster7x7(jet,t7)
def Cluster25x25(jet,towers):
    max_tower = getMaxTowerIn4x4('1313',towers)
    if max_tower['key'] == -1: return
    seed = towers[max_tower['key']]
    if max_tower['et'] < 2.5: return
    t21 = get21x21Towers(max_tower['key'])

    check = getMaxTowerIn7x7('1111',t21)
    check_seed = t21[check.key]
    if (check_seed != seed): return
    Cluster21x21(jet,t21)
def ClusterGCT(jetlist,towers):
    for ieta in irange(2,34,4):
        for iphi in irange(6,26,4):
            center = getKey(iphi,ieta)
            jet = Jet()
            t25 = get25x25Towers(center,towers)
            Cluster25x25(jet,t25)
            if jet.et > 0:
                jetlist.append(jet)
def ClusterCalorimeter(jetlist,towers):
    for iphi in irange(1,72,24):
        center = getKey(iphi,1)
        gct = getGCTTowers(center,towers)
        ClusterGCT(jetlist,gct)
def produce(calo):
    jetlist = []
    ClusterCalorimeter(jetlist,calo)
    jetlist.sort(key=lambda jet:jet.et,reverse=True)

    output = TFile("jet_output.root","RECREATE")
    h_seed = TH2F("seed","Clustered Jets;iPhi;iEta",72,0,72,34,0,34);
    h_tower = TH2F("tower","Jet Towers;iPhi;iEta",72,0,72,34,0,34);
    h_calo = TH2F("calo","Calo Towers;iPhi;iEta",72,0,72,34,0,34);
    for iphi in irange(1,72):
        h_seed.GetXaxis().SetBinLabel(iphi,str(iphi));
        h_tower.GetXaxis().SetBinLabel(iphi,str(iphi));
        h_calo.GetXaxis().SetBinLabel(iphi,str(iphi));
    for ieta in irange(1,34):
        h_seed.GetYaxis().SetBinLabel(ieta,str(ieta));
        h_tower.GetYaxis().SetBinLabel(ieta,str(ieta));
        h_calo.GetYaxis().SetBinLabel(ieta,str(ieta));
    for jet in jetlist:
        # if (debug1) printf("Writing Jet Seed iPhi: %i iEta: %i Et: %f\n",jet.iphi(),jet.ieta(),jet.jetClusterET);
        h_seed.SetBinContent(jet.iphi,jet.ieta,jet.et);
        for tower in jet.towers:
            if (tower.et == 0): continue
            # if (debug1) printf("--Writing Jet Tower iPhi: %i iEta: %i Et: %f\n",tower.iphi(),tower.ieta(),tower.total_tower_et);
            h_tower.SetBinContent(tower.iphi,tower.ieta,tower.et);
    for key,tower in calo.iteritems():
        h_calo.SetBinContent(tower.iphi,tower.ieta,tower.et);
    h_seed.Write();
    h_tower.Write();
    h_calo.Write();
    output.Close()
    return

events = Info('l1caloJetTest.root')
for event in events:
    print 'Analyzing Event: %i' % events.ievent
    produce(event)
