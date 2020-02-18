#ifndef EVENTINFO_H
#define EVENTINFO_H
#include "Jet.h"

struct EventInfo {
  TFile* tfile;
  int nevents;
  EventInfo(const char* input) {
    tfile = TFile::Open(input);
  }
  map<int,SimpleCaloHit> getCalo(int ievent) {
    TH2F* calo = (TH2F*)tfile->Get( ("calo;"+to_string(ievent)).c_str() );
    map<int,SimpleCaloHit> calomap;
    for (int ieta = 1; ieta <= 34; ieta++) {
      for (int iphi = 1; iphi <= 72; iphi++) {
	int key = getKey(iphi,ieta);
	calomap[key].init(iphi,ieta,calo->GetBinContent(iphi,ieta));
      }
    }
    return calomap;
  }
};
#endif
