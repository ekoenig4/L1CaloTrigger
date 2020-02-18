#include "Jet.h"
#include "DebugInfo.h"
#include "EventInfo.h"

struct Algorithm {
  float EtMinForSeedHit = 2.5;
  bool debug = true;
  bool useLargerGCT = true;
  bool doOverlap = true;

  DebugInfo debuginfo;

  void produce(EventInfo info,int ievent) {
    auto l1CaloTowers = info.getCalo(ievent);
    vector<l1CaloJetObj> l1CaloJetObjs;
    ClusterCalorimeter(l1CaloJetObjs, l1CaloTowers);
    std::sort(begin(l1CaloJetObjs), end(l1CaloJetObjs), [](const l1CaloJetObj& a,
							   const l1CaloJetObj& b){return a.jetClusterET > b.jetClusterET;});
    debuginfo.save(l1CaloJetObjs,l1CaloTowers);
  }
  void SetJetSeed(l1CaloJetObj& caloJetObj,SimpleCaloHit& l1CaloTower) {
    if (debug) cout << "---Seeding Jet " << l1CaloTower.toString() << endl;
    caloJetObj.seed_iPhi = l1CaloTower.tower_iPhi;
    caloJetObj.seed_iEta = l1CaloTower.tower_iEta;
    caloJetObj.jetClusterET += l1CaloTower.total_tower_et;
    caloJetObj.towers.push_back(l1CaloTower);
  }
  void AddTower(l1CaloJetObj& caloJetObj,SimpleCaloHit& l1CaloTower) {
    if (l1CaloTower.total_tower_et > 0) {
      if (debug) cout << "----Adding Tower " << l1CaloTower.toString() << endl;
      caloJetObj.jetClusterET += l1CaloTower.total_tower_et;
      caloJetObj.towers.push_back(l1CaloTower);
    }
  }

  void ClusterCalorimeter(vector<l1CaloJetObj>& caloJetObjs,map<int,SimpleCaloHit>& l1CaloTowers) {
    auto get_iphi = [](int x) { return 24*x + 2; };
    // for (int iphi = 1; iphi <= 72; iphi += 24) {
    for (int nGCT = 0; nGCT < 3; nGCT++) {
      int iphi = get_iphi(nGCT);
      int center = getKey(iphi,1);
      auto gct = getGCTTowers(center,l1CaloTowers);
      if (debug) {
	printf("Clustering GCT %i\n",center);
      }
      ClusterGCT(caloJetObjs,gct);
    }
  }
  void ClusterGCT(vector<l1CaloJetObj>& caloJetObjs,map<int,SimpleCaloHit>& l1CaloTowers) {
    int START_IPHI,GCT_IPHI;
    if ( !useLargerGCT ){
      // Standard 32x34 (iPhi x iEta)
      // Map center tower to (5,1) of a 32x34 box
      START_IPHI = 6; GCT_IPHI = 32;
    } else {
      // Larger 42x34 (iPhi x iEta)
      // Map center tower to (10,1) of a 42x34 box
      START_IPHI = 11; GCT_IPHI = 42;
    }

    auto get_ieta = [](int y) { return 4*y + 2; };
    auto get_iphi = [START_IPHI](int x) { return 4*x + START_IPHI; };
    // for (int ieta = 2; ieta <= 34; ieta += 4) {
      // for (int iphi = START_IPHI; iphi <= GCT_IPHI; iphi += 4) {
    for (int y = 0; y < 9; y++) {
      for (int x = 0; x < 6; x++) {
	int ieta = get_ieta(y);
	int iphi = get_iphi(x);
	
	l1CaloJetObj jet;
	int center = getKey(iphi,ieta);
	auto towers22x22 = get22x22Towers(center,l1CaloTowers);
	if (debug) {
	  printf("-Clustering 22x22 Jet %i\n",center);
	}
	Cluster22x22(jet,towers22x22);
	if ( jet.jetClusterET > 0 ) caloJetObjs.push_back(jet);
      }
    }
  }
  void Cluster22x22(l1CaloJetObj& caloJetObj,map<int,SimpleCaloHit>& l1CaloTowers) {
    int max_key = getMaxTowerIn4x4(1112,l1CaloTowers);
    auto& seed = l1CaloTowers[max_key];
    if (debug) cout << "--Max Tower " << seed.toString() << endl;
    if (seed.total_tower_et < EtMinForSeedHit) return;
    auto towers19x19 = get19x19Towers(max_key,l1CaloTowers);

    int check_key = getMaxTowerIn7x7(1010,towers19x19);
    auto& check_seed = towers19x19[check_key];
    if (check_seed != seed) return;
    if (debug){
      cout << "--Found Jet Seed " << seed.toString() << endl;
    }

    Cluster19x19(caloJetObj,towers19x19);
  }
  void Cluster19x19(l1CaloJetObj& caloJetObj,map<int,SimpleCaloHit>& l1CaloTowers) {
    auto towers7x7 = get7x7Towers(1010,l1CaloTowers);
    if ( doOverlap ){
      auto overlaps = getOverlapTowers(1010,l1CaloTowers);
      Cluster7x7(caloJetObj,towers7x7,overlaps);
    } else {
      Cluster7x7(caloJetObj,towers7x7);
    }
  }
  void Cluster7x7(l1CaloJetObj& caloJetObj,map<int,SimpleCaloHit>& l1CaloTowers,vector<SimpleCaloHit> overlaps=vector<SimpleCaloHit>()) {
    auto& seed = l1CaloTowers[404];
    SetJetSeed(caloJetObj,seed);
    for (auto& pair : l1CaloTowers) {
      int key = pair.first;
      auto& tower = pair.second;
      // Don't add the seed tower twice and make sure the tower isn't in the overlap list
      if ( key != 404 ) {
	bool taken = false;
	for (auto& overlap : overlaps) {
	  int dphi = tower_diPhi(overlap.tower_iPhi,tower.tower_iPhi);
	  int deta = tower_diEta(overlap.tower_iEta,tower.tower_iEta);
	  taken = abs(dphi) <= 3 && abs(deta) <= 3;
	  if (taken) break;
	}
	if (!taken) AddTower(caloJetObj,tower);
      }
    }
    if (debug) cout << "---Clustered Jet " << caloJetObj.toString() << endl << endl;
  }

  map<int,SimpleCaloHit> getGCTTowers(int center,map<int,SimpleCaloHit>& l1CaloTowers) {
    auto coordinates = getPhiEta(center);
    int cphi = coordinates.first; int ceta = coordinates.second;
    int dphi,deta,GCT_IPHI;
    if ( !useLargerGCT) {
      // Standard 32x34 (iPhi x iEta)
      // Map center tower to (5,1) of a 32x34 box
      dphi = cphi - 5; deta = ceta - 1; GCT_IPHI = 32;
    } else {
      // Larger 42x34 (iPhi x iEta)
      // Map center tower to (10,1) of a 42x34 box
      dphi = cphi - 10; deta = ceta - 1; GCT_IPHI = 42;
    }
  
    map<int,SimpleCaloHit> towers;
    for (int ieta = 1; ieta <= 34; ieta++) {
      for (int iphi = 1; iphi <= GCT_IPHI; iphi++) {
	int rphi = module(iphi + dphi - 1,72)+1; int reta = module(ieta + deta - 1,34)+1;
	int rkey = getKey(rphi,reta);
	int nkey = getKey(iphi,ieta);
	if ( l1CaloTowers.find(rkey) != l1CaloTowers.end() ) {
	  towers[nkey] = l1CaloTowers[rkey];
	} 
      }
    }
    return towers;
  }
  map<int,SimpleCaloHit> get22x22Towers(int center,map<int,SimpleCaloHit>& l1CaloTowers) {
    // Map center tower to (11,12) of a 22x22 box
    auto coordinates = getPhiEta(center);
    int cphi = coordinates.first; int ceta = coordinates.second;
    int dphi = cphi - 11; int deta = ceta - 12;

    map<int,SimpleCaloHit> towers;
    for (int ieta = 1; ieta <= 22; ieta++) {
      for (int iphi = 1; iphi <= 22; iphi++) {
	int rphi = iphi + dphi; int reta = ieta + deta;
	int rkey = getKey(rphi,reta);
	if ( l1CaloTowers.find(rkey) != l1CaloTowers.end() ) {
	  int nkey = getKey(iphi,ieta);
	  towers[nkey] = l1CaloTowers[rkey];
	}
      }
    }
    return towers;
  }
  map<int,SimpleCaloHit> get19x19Towers(int center,map<int,SimpleCaloHit>& l1CaloTowers) {
    // Map center tower to (10,10) of a 19x19 box
    auto coordinates = getPhiEta(center);
    int cphi = coordinates.first; int ceta = coordinates.second;
    int dphi = cphi - 10; int deta = ceta - 10;

    map<int,SimpleCaloHit> towers;
    for (int ieta = 1; ieta <= 19; ieta++) {
      for (int iphi = 1; iphi <= 19; iphi++) {
	int rphi = iphi + dphi; int reta = ieta + deta;
	int rkey = getKey(rphi,reta);
	if ( l1CaloTowers.find(rkey) != l1CaloTowers.end() ) {
	  int nkey = getKey(iphi,ieta);
	  towers[nkey] = l1CaloTowers[rkey];
	}
      }
    }
    return towers;
  }
  map<int,SimpleCaloHit> get7x7Towers(int center,map<int,SimpleCaloHit>& l1CaloTowers) {
    // Map center tower to (4,4) of a 7x7 box
    auto coordinates = getPhiEta(center);
    int cphi = coordinates.first; int ceta = coordinates.second;
    int dphi = cphi - 4; int deta = ceta - 4;

    map<int,SimpleCaloHit> towers;
    for (int ieta = 1; ieta <= 7; ieta++) {
      for (int iphi = 1; iphi <= 7; iphi++) {
	int rphi = iphi + dphi; int reta = ieta + deta;
	int rkey = getKey(rphi,reta);
	if ( l1CaloTowers.find(rkey) != l1CaloTowers.end() ) {
	  int nkey = getKey(iphi,ieta);
	  towers[nkey] = l1CaloTowers[rkey];
	}
      }
    }
    return towers;
  }
  vector<SimpleCaloHit> getOverlapTowers(int center,map<int,SimpleCaloHit>& l1CaloTowers) {
    // Get list of all potential seed towers that are larger than the center seed
    // There should only be on average 7.5 jets in the 3 ring around the 7x7 jet
    auto& center_tower = l1CaloTowers[center];
    auto coordinates = getPhiEta(center);
    int cphi = coordinates.first; int ceta = coordinates.second;
    int dphi = cphi - 7; int deta = ceta - 7;

    vector<SimpleCaloHit> overlap;
    for (int ieta = 1; ieta <= 13; ieta++) {
      for (int iphi = 1; iphi <= 13; iphi++) {
	// Ignore center 7x7
	if ( 4 <= ieta && ieta <=10 && 4 <= iphi && iphi <= 10 ) continue;

	int rphi = iphi + dphi; int reta = ieta + deta;
	int key = getKey(rphi,reta);
	if ( l1CaloTowers.find(key) != l1CaloTowers.end() ) {
	  auto tower = l1CaloTowers[key];
	  if ( center_tower < tower ) {
	    int check_key = getMaxTowerIn7x7(key,l1CaloTowers);
	    if ( key == check_key ) {
	      if (debug) cout << "--Overlaping with Tower " << tower.toString() << endl;
	      overlap.push_back(tower);
	    }
	  }
	}
      }
    }
    return overlap;
  }
  int getMaxTowerIn4x4(int center,map<int,SimpleCaloHit>& l1CaloTowers) {
    // Map center tower to (2,3) of a 4x4 box
    auto coordinates = getPhiEta(center);
    int cphi = coordinates.first; int ceta = coordinates.second;
    int dphi = cphi - 2; int deta = ceta - 3;
  
    pair<int,SimpleCaloHit> max_tower;
    for (int ieta = 1; ieta <= 4; ieta++) {
      for (int iphi = 1; iphi <= 4; iphi++) {
	int rphi = iphi + dphi; int reta = ieta + deta;
	int key = getKey(rphi,reta);
	if ( l1CaloTowers.find(key) != l1CaloTowers.end() ) {
	  auto& tower = l1CaloTowers[key];
	  if ( max_tower.second < tower ) {
	    max_tower = make_pair(key,tower);
	  }
	}
      }
    }
    return max_tower.first;
  }

  int getMaxTowerIn7x7(int center,map<int,SimpleCaloHit>& l1CaloTowers) {
    // Map center tower to (4,4) of a 7x7 box
    auto coordinates = getPhiEta(center);
    int cphi = coordinates.first; int ceta = coordinates.second;
    int dphi = cphi - 4; int deta = ceta - 4;
  
    pair<int,SimpleCaloHit> max_tower;
    for (int ieta = 1; ieta <= 7; ieta++) {
      for (int iphi = 1; iphi <= 7; iphi++) {
	int rphi = iphi + dphi; int reta = ieta + deta;
	int key = getKey(rphi,reta);
	if ( l1CaloTowers.find(key) != l1CaloTowers.end() ) {
	  auto& tower = l1CaloTowers[key];
	  if ( max_tower.second < tower ) {
	    max_tower = make_pair(key,tower);
	  }
	}
      }
    }
    return max_tower.first;
  }
};
