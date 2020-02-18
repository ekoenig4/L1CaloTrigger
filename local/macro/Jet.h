#ifndef JET_H
#define JET_H
inline int getKey(int iphi,int ieta) { return 100*iphi + ieta; }
inline std::pair<int,int> getPhiEta(int tower_key) { return std::make_pair((int)tower_key/100,(int)tower_key%100); }
int tower_diPhi( int iPhi_1, int iPhi_2 ) {
  // 360 Crystals in full, 72 towers, half way is 36
  int PI = 36;
  int result = iPhi_1 - iPhi_2;
  while (result > PI) result -= 2*PI;
  while (result <= -PI) result += 2*PI;
  return result;
}


// Added b/c of the iEta jump from +1 to -1 across the barrel mid point
int tower_diEta( int iEta_1, int iEta_2 ) {
  // On same side of barrel
  if (iEta_1 * iEta_2 > 0) return iEta_1 - iEta_2;
  else return iEta_1 - iEta_2 - 1;
}

int module(int n,int m) {
  int result = n%m;
  if (result < 0) result += m;
  return result;
}

struct SimpleCaloHit {
  int tower_iPhi;
  int tower_iEta;
  float total_tower_et;
  SimpleCaloHit() {
    tower_iPhi = -99;
    tower_iEta = -99;
    total_tower_et = -99;
  }
  void init(int iphi,int ieta,float et) {
    tower_iPhi = iphi;
    tower_iEta = ieta;
    total_tower_et = et;
  }
  inline int iphi() { return tower_iPhi; }
  inline int ieta() { return tower_iEta; }
  inline bool operator ==(const SimpleCaloHit& rhs) {
    return tower_iEta == rhs.tower_iEta && tower_iPhi == rhs.tower_iPhi && total_tower_et == rhs.total_tower_et;
  }
  inline bool operator !=(const SimpleCaloHit& rhs) {
    return !(*this == rhs);
  }
  inline bool operator <(const SimpleCaloHit& rhs) {
    if ( total_tower_et == rhs.total_tower_et ) {
      int dphi = tower_diPhi(rhs.tower_iPhi,tower_iPhi);
      int deta = tower_diEta(rhs.tower_iEta,tower_iEta);
      int mag = dphi + deta;
      if (mag == 0) return deta < 0;
      return mag > 0;
    }
    return total_tower_et < rhs.total_tower_et;
  }
  string toString() const {
    return "iPhi: "+to_string(tower_iPhi)+" iEta: "+to_string(tower_iEta)+" Et: "+to_string(total_tower_et); 
  }
};

struct l1CaloJetObj {
  int seed_iPhi = -99;
  int seed_iEta = -99;
  float jetClusterET = 0;
  vector<SimpleCaloHit> towers;
  l1CaloJetObj() {
    seed_iPhi = -99;
    seed_iEta = -99;
    jetClusterET = 0;
    towers.clear();
  }
  inline bool equals(const l1CaloJetObj& rhs) { return this->seed_iPhi == rhs.seed_iPhi && this->seed_iEta == rhs.seed_iEta; }
  inline int iphi() { return seed_iPhi; }
  inline int ieta() { return seed_iEta; }
  string toString() {
    return ( "iPhi: "+to_string(seed_iPhi)+" iEta: "+to_string(seed_iEta)+" Et: "+to_string(jetClusterET) ); 
  }
};
  
#endif
