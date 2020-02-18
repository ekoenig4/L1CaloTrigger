#include "Jet.h"
#include <unistd.h>

TBox* getBox(int iphi,int ieta,pair<float,float> geometry={7,7},int color = kBlack) {
  float dphi = geometry.first/2; float deta = geometry.second/2;
  TBox* box = new TBox(iphi-dphi-1,ieta-deta-1,iphi+dphi,ieta+deta);
  box->SetFillStyle(0);
  box->SetLineColor(color);
  box->SetLineWidth(2);
  return box;
}

struct Grid {
  float cphi = 0;
  float ceta = 0;
  int xmax = 0;
  int ymax = 0;
  TH2F* grid;
  vector<TBox*> boxes;
  Grid() {}
  Grid(string name,map<int,SimpleCaloHit> gridmap,int center=0) {
    auto coordinates = getPhiEta(center);
    cphi = coordinates.first; ceta = coordinates.second;
    if (name == "t22"){
      xmax = 22; ymax = 22;
    } else if (name == "t19"){
      xmax = 19; ymax = 19;
    } else if (name == "t7"){
      xmax = 7; ymax = 7;
    } else {
      for (auto key : gridmap) {
	coordinates = getPhiEta(key.first);
	if ( xmax < coordinates.first ) xmax = coordinates.first;
	if ( ymax < coordinates.second) ymax = coordinates.second;
      }
    }
    grid = new TH2F( ("grid_"+name).c_str(),";iPhi;iEta",xmax,0,xmax,ymax,0,ymax);
    for (int i = 1; i < xmax; i++) grid->GetXaxis()->SetBinLabel(i,to_string(i).c_str());
    for (int i = 1; i < ymax; i++) grid->GetYaxis()->SetBinLabel(i,to_string(i).c_str());

    for (int y = 1; y < ymax; y++) {
      for (int x = 1; x < xmax; x++) {
	int key = getKey(x,y);
	if ( gridmap.find(key) != gridmap.end() ) {
	  float et = gridmap[key].total_tower_et;
	  grid->SetBinContent(x,y,et);
	}
      }
    }
    boxes.clear();
  }
};
vector<string> heirarchy = {"t7","t19","t22","gct","calo"};
pair<float,float> t19_to_t22(pair<float,float> start,pair<float,float> center) {
  return { start.first - 10 + center.first,start.second - 10 + center.second };
}

pair<float,float> t22_to_gct(pair<float,float> start,pair<float,float> center) {
  return { start.first - 11 + center.first,start.second - 12 + center.second };
}

pair<float,float> gct_to_calo(pair<float,float> start,pair<float,float> center) {
  return { start.first - 10 + center.first,start.second - 1 + center.second };
}

struct DebugInfo {
  TFile* output;
  TH2F *h_calo,*h_seed,*h_tower;

  map<string,Grid> gridlist;
  map<string,TCanvas*> canvaslist;
  
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
  void add(string name,map<int,SimpleCaloHit> gridmap,int center=0) {
    if ( gridlist.find(name) != gridlist.end() ) {
      auto it = find(heirarchy.begin(),heirarchy.end(),name);
      ++it;
      if (it != heirarchy.end()) {
	for (auto box : gridlist[name].boxes) {
	  for (it; it != heirarchy.end(); ++it) {
	    if ( gridlist[*it].boxes.size() > 0 ) {
	      gridlist[*it].boxes.pop_back();
	    }
	  }
	}
      }
    }
    gridlist[name] = Grid(name,gridmap,center);
    if (name == "t22") addBox(name,11.5,11.5,3.0,3.0);
    if (name == "t19") addBox(name,10,10,7,7);
    if (canvaslist.find(name) != canvaslist.end())
      canvaslist[name]->Delete();
  }
  void addBox(string name,float cx,float cy,float wx,float wy) {
    auto grid = gridlist[name];
    grid.boxes.push_back( getBox(cx,cy,make_pair(wx,wy)) );

    
    pair<float,float> nr;
    auto it = find(heirarchy.begin(),heirarchy.end(),name);
    ++it;
    if (it != heirarchy.end()) {
      if (name == "t19")
	nr = t19_to_t22(make_pair(cx,cy),make_pair(wx,wy));
      else if (name == "t22")
	nr = t22_to_gct(make_pair(cx,cy),make_pair(wx,wy));
      else if (name == "gct")
	nr = gct_to_calo(make_pair(cx,cy),make_pair(wx,wy));
      addBox(*it,nr.first,nr.second,wx,wy);
    }
  }
  void draw(string name,map<int,SimpleCaloHit> gridmap=map<int,SimpleCaloHit>(),bool pause=false) {
    if ( gridlist.find(name) == gridlist.end() && gridmap.size() > 0 ) add(name,gridmap);

    auto grid = gridlist[name];

    if (grid.boxes.size() > 0 || name == "t7") {
      auto it = find(heirarchy.begin(),heirarchy.end(),name);
      ++it;
      if (it != heirarchy.end()) draw(*it);
    }
    auto hs = grid.grid;
    TCanvas* c = new TCanvas( ("canvas_"+name).c_str(),name.c_str() );
    c->cd();
    c->SetGrid();
    hs->Draw("COLZ");
    for (auto box : grid.boxes) box->Draw("same");
    canvaslist[name] = c;
    char h;
    if (pause) cin >> h;
    if (pause) usleep( 1000000 );
  }
};
