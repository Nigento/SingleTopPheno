#include  "../include/root_simu.hpp"
#include  "../include/matrixSingleTop.hpp"
#include <iostream>
#include <vector>

using namespace std;

int main (){

  TLorentzVector tw (-4.81086e-06, 3.67358e-06, -0.0330773, -0.0215645);
  TLorentzVector tbarw ( -6.42135e-06, 3.8503e-06, 0.0332121, 0.0217606);
  TLorentzVector tq (2.88377e-12, 1.0515e-12, 3.73404e-08, 4.39321e-08);
  TLorentzVector tbarq (1.41135e-11, 4.08558e-12, 5.47861e-08, 5.01681e-08);


  TLorentzVector bmu (100,0.,0.,0.);

  double stw = bmu*tw;
  double stbarw = bmu*tbarw;
  double stq = bmu*tq;
  double stbarq = bmu*tbarq;

  double asymtw = (stw*2.981 - stbarw*2.981)/(stw*2.981 + stbarw*2.981);
  cout<<"Asymmetry for tw-chan is:"<<asymtw<<endl;
  double asymtq = (stq*12.26 - stbarq*7.004)/(stq*12.26 + 7.004*stbarq);
  cout<<"Asymmetry for t-chan is:"<<asymtq<<endl;
  return 0;
}
