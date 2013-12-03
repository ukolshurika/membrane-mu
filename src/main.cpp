#include <iostream>
#include <fstream>
#include <cmath>

#include "membrane.h"
#include "bound.h"
#include "gaus.h"

using namespace std;


int main(){
  double h0 = 2/100.0;
  double q = 2.65*1000/88.3/1000000;
  double n = 3.4;

  Membrane m(q, h0, n);

  m.free(999);
  m.constrained(999);

  return 0;
}