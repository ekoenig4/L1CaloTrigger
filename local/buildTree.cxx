
struct TreeInfo {
  int nJet;
  vector<float> jetPhi;
  vector<float> jetEta;
  vector<float> jetEt;
  void init() {
    nJet = 0;
    jetPhi.clear();
    jetEta.clear();
    jetEt.clear();
  }
};

void buildTree() {
  const char* input = "data/tyler_output.root";
  TFile* tfile = TFile::Open(input);
  int nevents = 100;

  TreeInfo info;
  TFile* output = new TFile("tyler_tree.root","recreate");
  TTree* tree = new TTree("tree","tree");
  tree->Branch("nJet",&info.nJet);
  tree->Branch("jetPhi",&info.jetPhi);
  tree->Branch("jetEta",&info.jetEta);
  tree->Branch("jetEt",&info.jetEt);

  for (int ievent = 1; ievent <= nevents; ievent++) {
    info.init();
    TH2F* h_seed = (TH2F*)tfile->Get( ("seed;"+to_string(ievent)).c_str() );

    for (int iphi = 1; iphi <= h_seed->GetNbinsX(); iphi++) {
      for (int ieta = 1; ieta <= h_seed->GetNbinsY(); ieta++) {
	if ( h_seed->GetBinContent(iphi,ieta) > 0 ) {
	  info.nJet++;
	  info.jetPhi.push_back(iphi);
	  info.jetEta.push_back(ieta);
	  info.jetEt.push_back( h_seed->GetBinContent(iphi,ieta) );
	}
      }
    }
    tree->Fill();
  }
  tree->Write();
  output->Close();
  tfile->Close();
}

  
