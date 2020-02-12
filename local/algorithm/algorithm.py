from ROOT import *
from sys import argv
from objects import *
from config import config
from debug import debugger

#----Algorithm----#
def ClusterCalorimeter(jetlist,towers):
    debugger.add('calo',towers)
    # debugger.draw('calo')
    for iphi in irange(1,72,24):
        center = getKey(iphi,1)
        gct = getGCTTowers(center,towers)
        debugger.add('gct',gct,center)
        if config.debug: print "Clustering GCT %s" % center
        ClusterGCT(jetlist,gct)
def ClusterGCT(jetlist,towers):
    # debugger.draw('gct')
    #-- Standard 34x32 --#
    # for ieta in irange(2,34,4):
        # for iphi in irange(6,26,4):
    
    #-- Larger 34x42 --#
    for ieta in irange(2,34,4):
        for iphi in irange(11,31,4):
            
            center = getKey(iphi,ieta)
            jet = Jet()
            t22 = get22x22Towers(center,towers)
            debugger.add('t22',t22,center)
            if config.debug: print "-Clustering 22x22 Jet %s" % center
            Cluster22x22(jet,t22)
            if jet.et > 0:
                jetlist.append(jet)
def Cluster22x22(jet,towers):
    max_tower = getMaxTowerIn4x4('1112',towers)
    if max_tower['key'] == -1: return
    seed = max_tower['tower']
    if config.debug: print "--Max Tower %s" % str(seed)
    # debugger.draw('t22')
    if seed.et < 2.5: return
    t19 = get19x19Towers(max_tower['key'],towers)

    check = getMaxTowerIn7x7('1010',t19)
    check_seed = t19[check['key']]
    if (check_seed != seed): return
    if config.debug: print "--Found Jet Seed %s" % str(seed)
    debugger.add('t19',t19,max_tower['key'])
    Cluster19x19(jet,t19)
def Cluster19x19(jet,towers):
    #--Overlap Handling needs to be done here--#
    overlaps = getOverlapTowers('1010',towers)
    # debugger.draw('t19',pause=True)
    t7 = get7x7Towers('1010',towers,overlaps)
    debugger.add('t7',t7,'1010')
    Cluster7x7(jet,t7)
def Cluster7x7(jet,towers):
    debugger.draw('t7',pause=True)
    seed = towers['404']
    jet.seed(seed)
    for key,tower in towers.iteritems():
        if key == '404': continue
        jet.add(tower)
    if config.debug: print '---Clustered Jet %s\n' % str(jet)
   
#----Helper Methods----# 
def getGCTTowers(center_key,towers):
    #-- Standard 34x32 (iEta x iPhi) --#
    # subset = {}
    # cphi,ceta = getPhiEta(center_key)
    # dphi,deta = cphi - 5,ceta - 1
    # for ieta in irange(1,34):
    #     for iphi in irange(1,32):
    #         rphi,reta = (iphi + dphi - 1)%72+1,(ieta + deta - 1)%34+1
    #         rkey = getKey(rphi,reta)
    #         nkey = getKey(iphi,ieta)
    #         if rkey not in towers: continue
    #         subset[nkey] = towers[rkey]
    # return subset

    #-- Larger 34x42 (iEta x iPhi) --#
    subset = {}
    cphi,ceta = getPhiEta(center_key)
    dphi,deta = cphi - 10,ceta - 1
    for ieta in irange(1,34):
        for iphi in irange(1,42):
            rphi,reta = (iphi + dphi - 1)%72+1,(ieta + deta - 1)%34+1
            rkey = getKey(rphi,reta)
            nkey = getKey(iphi,ieta)
            if rkey not in towers: continue
            subset[nkey] = towers[rkey]
    return subset
def get22x22Towers(center_key,towers):
    subset = {}
    cphi,ceta = getPhiEta(center_key)
    dphi,deta = cphi - 11,ceta - 12
    for ieta in irange(1,22):
        for iphi in irange(1,22):
            rphi,reta = (iphi + dphi),(ieta + deta)
            rkey = getKey(rphi,reta)
            nkey = getKey(iphi,ieta)
            if rkey not in towers: continue
            subset[nkey] = towers[rkey]
    return subset
def get19x19Towers(center_key,towers):
    subset = {}
    cphi,ceta = getPhiEta(center_key)
    dphi,deta = cphi - 10,ceta - 10
    for ieta in irange(1,19):
        for iphi in irange(1,19):
            rphi,reta = (iphi + dphi),(ieta + deta)
            rkey = getKey(rphi,reta)
            nkey = getKey(iphi,ieta)
            if rkey not in towers: continue
            subset[nkey] = towers[rkey]
    return subset
def get7x7Towers(center_key,towers,overlaps=[]):
    subset = {}
    cphi,ceta = getPhiEta(center_key)
    dphi,deta = cphi - 4,ceta - 4
    for ieta in irange(1,7):
        for iphi in irange(1,7):
            rphi,reta = (iphi + dphi),(ieta + deta)
            rkey = getKey(rphi,reta)
            nkey = getKey(iphi,ieta)
            if rkey not in towers or towers[rkey].et == 0: continue
            taken = False
            for okey in overlaps:
                ophi,oeta = getPhiEta(okey)
                d_iPhi,d_iEta = deltaR(rphi,reta,ophi,oeta)
                taken = abs(d_iPhi) <= 3 and abs(d_iEta) <= 3
                if taken: break
            if not taken:
                subset[nkey] = towers[rkey]
    return subset
def getOverlapTowers(center_key,towers):
    # Get List of all potential seed towers that are large than the center seed
    # There should only be on average 7.5 jets in the 3 ring around the 7x7 jet
    max_towers = []
    cphi,ceta = getPhiEta(center_key)
    dphi,deta = cphi - 7,ceta - 7

    ignore_iphi = irange(4,10)
    ignore_ieta = irange(4,10)
    
    for ieta in irange(1,13):
        for iphi in irange(1,13):
            if 4 <= ieta and ieta <= 10 and 4 <= iphi and iphi <= 10: continue

            rphi,reta = (iphi + dphi),(ieta + deta)
            rkey = getKey(rphi,reta)
            if rkey not in towers: continue
            if towers[rkey] > towers[center_key]:
                check = getMaxTowerIn7x7(rkey,towers)
                if check['key'] == rkey:
                    if config.debug: print '--Overlaping tower %s'  % str(towers[rkey])
                    debugger.addBox('t19',box=(rphi,reta,7,7))
                    max_towers.append(rkey)
    return max_towers
def getMaxTowerIn4x4(center,towers):
    cphi,ceta = getPhiEta(center)
    dphi,deta = cphi - 2,ceta - 3
    max_tower = {'key':-1,'tower':None}
    for iphi in irange(1,4):
        for ieta in irange(1,4):
            rphi,reta = iphi + dphi,ieta + deta
            key = getKey(rphi,reta)
            if key not in towers: continue
            if max_tower['tower'] < towers[key]:
                max_tower['key'] = key
                max_tower['tower'] = towers[key]
    return max_tower
def getMaxTowerIn7x7(center,towers):
    cphi,ceta = getPhiEta(center)
    dphi,deta = cphi - 4,ceta - 4
    max_tower = {'key':-1,'tower':None}
    for iphi in irange(1,7):
        for ieta in irange(1,7):
            rphi,reta = iphi + dphi,ieta + deta
            key = getKey(rphi,reta)
            if key not in towers: continue
            if max_tower['tower'] < towers[key]:
                max_tower['key'] = key
                max_tower['tower'] = towers[key]
    return max_tower
