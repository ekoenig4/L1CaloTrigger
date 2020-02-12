#!/usr/bin/env python
from ROOT import *
from algorithm import *

output = TFile("jet_output.root","RECREATE")
def produce(calo):
    jetlist = []
    ClusterCalorimeter(jetlist,calo)
    jetlist.sort(key=lambda jet:jet.et,reverse=True)

    output.cd()
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

    if config.debug: print "Clustered %i Jets" % len(jetlist)
    for jet in jetlist:
        if (config.debug): print "Writing Jet Seed iPhi: %i iEta: %i Et: %f" % (jet.iphi,jet.ieta,jet.et);
        h_seed.SetBinContent(jet.iphi,jet.ieta,jet.et);
        for tower in jet.towers:
            if (tower.et == 0): continue
            if (config.debug): print "--Writing Jet Tower iPhi: %i iEta: %i Et: %f" % (tower.iphi,tower.ieta,tower.et);
            if h_tower.GetBinContent(tower.iphi,tower.ieta) > 0: print 'WARNING. Overlapped Tower %s'%tower
            h_tower.SetBinContent(tower.iphi,tower.ieta,tower.et);
    for key,tower in calo.iteritems():
        h_calo.SetBinContent(tower.iphi,tower.ieta,tower.et);
    h_seed.Write();
    h_tower.Write();
    h_calo.Write();
    return

events = Info('data/tyler_output.root')
for ievent in range(events.nevents):
    calo = events.getCalo(ievent)
    print 'Analyzing Event: %i' % ievent
    produce(calo)
output.Close()
