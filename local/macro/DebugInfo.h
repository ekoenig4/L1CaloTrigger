#include "Jet.h"
#include <unistd.h>

struct DebugInfo {
  TFile* output;
  TH2F *h_calo,*h_seed,*h_tower;
  
  DebugInfo() {
    output = new TFile("output_macro.root","recreate");
    h_seed = new TH2F("seed","Clustered Jets;iPhi;iEta",72,0,72,34,0,34);
    h_tower = new TH2F("tower","Jet Towers;iPhi;iEta",72,0,72,34,0,34);
    h_calo = new TH2F("calo","Calo Towers;iPhi;iEta",72,0,72,34,0,34);
    for (int iphi = 1; iphi <= 72; iphi++) {
      h_seed->GetXaxis()->SetBinLabel(iphi,std::to_string(iphi).c_str());
      h_tower->GetXaxis()->SetBinLabel(iphi,std::to_string(iphi).c_str());
      h_calo->GetXaxis()->SetBinLabel(iphi,std::to_string(iphi).c_str());
    }
    for (int ieta = 1; ieta <= 34; ieta++) {
      h_seed->GetYaxis()->SetBinLabel(ieta,std::to_string(ieta).c_str());
      h_tower->GetYaxis()->SetBinLabel(ieta,std::to_string(ieta).c_str());
      h_calo->GetYaxis()->SetBinLabel(ieta,std::to_string(ieta).c_str());
    }
  }
  void save(vector<l1CaloJetObj>& caloJetObjs,std::map<int,SimpleCaloHit>& l1CaloTowers) {
    TH2F* h_seed = (TH2F*)this->h_seed->Clone(); h_seed->Reset();
    TH2F* h_tower = (TH2F*)this->h_tower->Clone(); h_tower->Reset();
    TH2F* h_calo = (TH2F*)this->h_calo->Clone(); h_calo->Reset();
    for (auto jet : caloJetObjs) {
      // if (debug1) printf("Writing Jet Seed iPhi: %i iEta: %i Et: %f\n",jet.iphi(),jet.ieta(),jet.jetClusterET);
      h_seed->SetBinContent(jet.iphi(),jet.ieta(),jet.jetClusterET);
      for (auto tower : jet.towers) {
	if (tower.total_tower_et == 0) continue;
	// if (debug1) printf("--Writing Jet Tower iPhi: %i iEta: %i Et: %f\n",tower.iphi(),tower.ieta(),tower.total_tower_et);
	if (h_tower->GetBinContent(tower.iphi(),tower.ieta()) > 0)
	  cout << "WARNING Jet Tower " << tower.toString() << endl;
	h_tower->SetBinContent(tower.iphi(),tower.ieta(),tower.total_tower_et);
      }
    }
    for (auto key : l1CaloTowers) {
      auto tower = key.second;
      h_calo->SetBinContent(tower.iphi(),tower.ieta(),tower.total_tower_et);
    }
    auto& cwd = gDirectory;
    output->cd();
    h_seed->Write();
    h_tower->Write();
    h_calo->Write();
    cwd->cd();
  }
};
