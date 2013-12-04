#include "membrane.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <algorithm> 

#include "utils.h"
#include "simpson.h"

using namespace std;

namespace{

const double kSqrt3 = sqrt(3);
const int kSimpsonStep = 99;

double ValueAsLine(double time,
                   const pair<double, double>& point1,
                   const pair<double, double>& point2) {
  double k, b;
  k = (point2.second - point1.second)/(point2.first - point1.first);
  b = -k*point1.first+point1.second;
  return k*time + b;
}


void Kahan(double offset, const vector<pair<double, double>>& src, vector<pair<double, double>>* dst){
  // DCHECK(dst->empty());
  // for (auto& v : src)
  //   DCHECK(!utils::IsNaN(v.first));
  // cerr << 'v' <<endl;
  // DCHECK(!utils::IsNaN(offset))
  double s = offset, c = 0, t, y;
  cerr << " #### " << s << endl;
  for(auto j=src.begin(); j!=src.end(); j++){
    y = j->first - c;
    t = s + y;
    c = (t - s) - y;
    s = t;
    // DCHECK(kSqrt3*s/2.0 > offset);
    // (*dst).push_back(make_pair(kSqrt3*s/2.0, j->second));
    (*dst).push_back(make_pair(s, j->second));
  }
}

}  // namespace

struct Free {
  Free(double h0, double q, double n):q_(q), h0_(h0), n_(n){};
  double operator()(double alpha) const{
    // cerr << "FREE" << endl;
    // cerr<< (1/alpha-1/tan(alpha)) << endl;
    return (1/alpha-1/tan(alpha))*pow(((2*h0_*sin(alpha)*sin(alpha))/(kSqrt3*q_*alpha)-1), n_);
  }

  double h(double alpha){
    return sin(alpha)/alpha*h0_;
  }

  private:
  double q_;
  double h0_;
  double n_;
};

Membrane::Membrane(double q, double h0, double n, double mu):q_(q), n_(n), h0_(h0), mu_(mu) {
  alpha1_ = 0.41; // from Maxima flexible step(boolshit just get it from terraud)
  alpha2_ = 0.93;//atan(2*Bound::kB/(Bound::kB2));//M_PI/2;
  h1_ = sin(alpha2_)/alpha2_*h0_;
  cerr << h1_ << ' ' << h1_/h0_ << endl;
}

void Membrane::free(int steps){
  double dalpha = (alpha2_ - alpha1_)/steps;
  double t;
  Free f(h0_, q_, n_);

  vector<pair<double, double>> v;

  for(double a = alpha1_; a < alpha2_; a+=dalpha){
    t = Simpson::Integrate(a, a+dalpha, kSimpsonStep, f);
    v.push_back(make_pair(t, a));
  }

  t_free_.clear();
  double offset = 0;
  cerr << " ### "<< offset << endl;
  Kahan(offset, v, &t_free_);
  // for (auto it = v.begin(); it != v.end(); ++it) {
  //   t_free_.push_back(make_pair((it->first + offset), it->second));
  //   offset += it->first;
  // }
}

void Membrane::constrained(int steps){
  ofstream h_data("data/h.dat");
  ofstream sigma_data ("data/s.dat");

  double init_sigma_k = q_/h1_; // Q*RHO/H (from free stadia)
  double dt = 10000000;

  vector <double> sigma_k(steps, init_sigma_k), sigma_k1(steps, 0.0),
                  ds_k(steps, 0.0), ds_k1(steps, 0.0), 
                  delta_ds_k(steps, 0.0), delta_ds_k1(steps, 0.0),
                  h_k(steps, 0.63), h_k1(steps, 0.0);

  for(auto t = 1; t < steps ; ++t){

    for(auto i = 1; i < t; ++i)
      delta_ds_k1[i] = pow((sigma_k[i-1]+sigma_k[i])/(4/sqrt(3)-(sigma_k[i-1]+sigma_k[i])), n_)*ds_k[i]*dt;

    for(auto i = 0; i< t; ++i)
      ds_k1[i] = ds_k[i] + delta_ds_k1[i];
    
    for(auto i = 1; i<t; ++i) 
      h_k1[i] = h_k[i]*(1-pow(1/(1-sqrt(3)/2*q_/h_k[t-1])-1, n_));
    
    /** WANING!!! THINK ABOUT HK+1K+1 !!!*/

    sigma_k1[t] = sqrt(3)/2*q_/h_k1[t];
    for(auto i = t-1; i>0; --i)
      sigma_k1[i] = sigma_k1[i+1]*h_k1[i+1]/h_k1[i]-mu_*ds_k1[i+1]*q_/h_k1[i];
    
    cerr<<h_k[t-1]<< endl;
    ds_k1[t] = M_PI_2*pow(1/(1-sqrt(3)/2*q_/h_k[t-1])-1, n_);
    
    h_data << t_free_.back().first + t*dt << ' ' << exp(-2/M_PI*ds_k1[t]) << endl;
    sigma_data << t_free_.back().first + t*dt << ' ' << sqrt(3)/2*q_/exp(-2/M_PI*ds_k1[t])<<endl;
      
    /* --- SWAPPING --- */
    swap(sigma_k, sigma_k1);
    swap(ds_k, ds_k1);
    swap(delta_ds_k, delta_ds_k1);
    swap(h_k, h_k1);
  }

}
