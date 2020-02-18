#include "EventInfo.h"
#include "Algorithm.h"

void run() {
  EventInfo info = EventInfo("data/tyler_output.root");
  Algorithm algorithm;
  algorithm.debug = false;
  for (int ievent = 1; ievent <= 100; ievent++) {
    cout << "Analyzing Event: " << ievent << endl;
    algorithm.produce(info,ievent);
  }
}
