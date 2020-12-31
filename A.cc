#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <functional>
using namespace std;
using ll = int64_t;
using ld = double;

//We want stats on some distributions over civs in C. For each distribution, can describe it in terms of a mean and percentile cuts: 1%, 5%, 10%, 20%, 50%, 80%, 90%, 95%, 99%. 
ld average(const vector<ld>& A) {
  ld sum = 0.0;
  for(auto& x : A) {
    sum += x;
  }
  return sum / static_cast<ld>(A.size());
}
template<typename U,typename T>
vector<U> map(const vector<T>& A, std::function<U(T)> fn) {
  vector<U> ans;
  for(auto& x : A) {
    ans.push_back(fn(x));
  }
  return ans;
}
ld percentile(const vector<ld>& C, ll p) {
  if(C.size() == 0) {
    return 0.0;
  } else if(C.size() == 1) {
    return C[0];
  } else {
    ll idx = static_cast<size_t>((C.size()-1)*(p/100.0));
    return C[idx];
  }
}

struct Civ {
  Civ(ll D) : V(D, 0.0), T(0.0) {}
  vector<ld> V; // position in space
  ld T; // origin time
};
struct Distro {
  Distro(const vector<ld>& A) :
    mean(average(A)),
    p1(percentile(A, 1)),
    p5(percentile(A, 5)),
    p10(percentile(A, 10)),
    p20(percentile(A, 20)),
    p50(percentile(A, 50)),
    p80(percentile(A, 80)),
    p90(percentile(A, 90)),
    p95(percentile(A, 95)),
    p99(percentile(A, 99)) {}
  Distro(ld m, ld q1, ld q5, ld q10, ld q20, ld q50, ld q80, ld q90, ld q95, ld q99) :
    mean(m), p1(q1), p5(q5), p10(q10), p20(q20), p50(q50), p80(q80), p90(q90), p95(q95), p99(q99) {}
  ld mean;
  ld p1;
  ld p5;
  ld p10;
  ld p20;
  ld p50;
  ld p80;
  ld p90;
  ld p95;
  ld p99;
};
ostream& operator<<(ostream& o, const Distro& d) {
  o << "mean=" << d.mean << " p1=" << d.p1 << " p5=" << d.p5 << " p10=" << d.p10 << " p20=" << d.p20 << " p50=" << d.p50 << " p80=" << d.p80 << " p90=" << d.p90 << " p95=" << d.p95 << " p99=" << d.p99;
  return o;
}
Distro combine_distros(const vector<Distro>& D) {
  return Distro {
    average(map<ld,Distro>(D, [](Distro d) { return d.mean; })),
    average(map<ld,Distro>(D, [](Distro d) { return d.p1; })),
    average(map<ld,Distro>(D, [](Distro d) { return d.p5; })),
    average(map<ld,Distro>(D, [](Distro d) { return d.p10; })),
    average(map<ld,Distro>(D, [](Distro d) { return d.p20; })),
    average(map<ld,Distro>(D, [](Distro d) { return d.p50; })),
    average(map<ld,Distro>(D, [](Distro d) { return d.p80; })),
    average(map<ld,Distro>(D, [](Distro d) { return d.p90; })),
    average(map<ld,Distro>(D, [](Distro d) { return d.p95; })),
    average(map<ld,Distro>(D, [](Distro d) { return d.p99; }))
  };
}
Distro distro(vector<ld>& T) {
  sort(T.begin(), T.end());
  if(T.size() == 0) {
    return Distro{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  } else {
    return Distro(T);
  }
}

struct Sample {
  ld C;
  Distro F;
  Distro G;
  ld percent_see;
  ld median_wait_see;
  ld median_wait_nosee;
  Distro H;
};
Sample combine_samples(const vector<Sample>& S) {
  return Sample{
    average(map<ld,Sample>(S, [](Sample s) { return s.C; })),
    combine_distros(map<Distro,Sample>(S, [](Sample s) { return s.F; })),
    combine_distros(map<Distro,Sample>(S, [](Sample s) { return s.G; })),
    average(map<ld,Sample>(S, [](Sample s) { return s.percent_see; })),
    average(map<ld,Sample>(S, [](Sample s) { return s.median_wait_see; })),
    average(map<ld,Sample>(S, [](Sample s) { return s.median_wait_nosee; })),
    combine_distros(map<Distro,Sample>(S, [](Sample s) { return s.H; }))
  };
}
ostream& operator<<(ostream& o, const Sample& S) {
  o << "|C|: " << S.C << endl;
  o << "F: " << S.F << endl;
  o << "G: " << S.G << endl;
  o << "percent_see: " << S.percent_see 
    << " median_wait_among_see: " << S.median_wait_see
    << " median_wait_nosee: " << S.median_wait_nosee
    << endl;
  o << "H: " << S.H << endl;
  return o;
}

ld r(ld lo, ld hi) {
  static default_random_engine RNG;
  return uniform_real_distribution<ld>(lo, hi)(RNG);
}

// Between every pair of points i,j is a Euclidean distance of di,j = ((xi-xj)2 + (yi-yj)2 + (zi-zj)2)1/2.
// Use “hypertorus” distance metric, that identifies opposite sides of box [0,1]d. (Can calc via, instead of dx = xi-xj, use (let a = abs(xi-xj), if a < ½, use a, else use 1-a).
ld distance(const vector<ld>& A, const vector<ld>& B) {
  ld d2 = 0.0;
  for(ll i=0; i<A.size(); i++) {
    ld dx = abs(A[i]-B[i]);
    dx = min(dx, 1.0-dx);
    d2 += dx*dx;
  }
  return sqrt(d2);
}

Sample simulate(ll D, ld speed, ld n, ll N) {
  vector<Civ> CIV(N, Civ{D});
	for(ll i=0; i<N; i++) {
    for(ll j=0; j<D; j++) {
      CIV[i].V[j] = r(0.0, 1.0);
    }
    CIV[i].T = pow(r(0.0, 1.0), 1.0/(1.0+n));
    // If save & reuse same r to generate vi,ti, can get less noisy comparisons between runs using different parameter values?
    // What does this mean, since there are multiple components to V[i]? Should V[i] have the same value in all dimensions?
	}

  //Colonies from civ at j first arrive at point i at time ai,j = tj + di,j/s. If ti > ai,j , then j colony prevents a civ origin at i. We seek subset C of points i that are not so prevented by any other point j.
  //(To find C fast, outer loop over all points i starting with highest ti, inner loop over all points j starting with lowest tj. If j prevents i then delete i from any further consideration in either loop.)
  sort(CIV.begin(), CIV.end(), [](const Civ& A, const Civ& B) { return A.T < B.T; });
  vector<Civ> ALIVE;
  for(ll i=N-1; i>=0; i--) {
    bool alive = true;
    for(ll j=0; j<N; j++) {
      ld aij = CIV[j].T + distance(CIV[j].V, CIV[i].V)/speed;
      if(i!=j && CIV[i].T > aij) {
        alive = false;
        break;
      }
    }
    if(alive) {
      ALIVE.push_back(CIV[i]);
    }
  }
  assert(ALIVE.size() > 0);

  vector<ld> T;
  vector<ld> MIN_WAIT;
  vector<ld> MIN_WAIT_SEE;
  vector<ld> MIN_WAIT_NOSEE;
  vector<ld> MAX_ANGLE;
  for(ll i=0; i<ALIVE.size(); i++) {
    T.push_back(ALIVE[i].T);

    bool see = false;
    ld min_wait = 1e9;
    ld max_angle = -1.0;
    for(ll j=0; j<ALIVE.size(); j++) {
      if(j!=i) {
        ld dij = distance(ALIVE[i].V, ALIVE[j].V);
        ld wij = (dij/speed - (ALIVE[i].T-ALIVE[j].T))/2.0;
        ld oij = ALIVE[j].T + dij;
        ld bij = speed*(ALIVE[i].T-oij)/dij;
        if(ALIVE[i].T > oij) {
          see = true;
          max_angle = max(max_angle, bij);
        }
        if(wij < min_wait) {
          min_wait = wij;
        }
      }
    }
    // max angle = max{j that i can see in C} bi,j.
    assert(min_wait < 1e6);
    MIN_WAIT.push_back(min_wait);
    if(see) {
      MIN_WAIT_SEE.push_back(min_wait);
      assert(max_angle >= 0.0);
      MAX_ANGLE.push_back(max_angle);
    } else {
      MIN_WAIT_NOSEE.push_back(min_wait);
    }
  }
  ld percent_see = static_cast<ld>(MIN_WAIT_SEE.size())/static_cast<ld>(ALIVE.size());
  return Sample{static_cast<ld>(ALIVE.size()), distro(T), distro(MIN_WAIT), percent_see, percentile(MIN_WAIT_SEE, 50), percentile(MIN_WAIT_NOSEE, 50), distro(MAX_ANGLE)};
}

// For same value of N/s3, should get about same |C|/s3, same distributions over C of ti , wi,j etc. 


int main() {
  for(auto speed : vector<ld>{0.5,0.25,0.125,1.0/16.0}) {
    for(auto N : vector<ll>{8000,4000,2000,1000}) {
      vector<Sample> S1;
      for(ll t=0; t<1000; t++) {
        Sample S = simulate(1, speed, 10.0, N);
        S1.push_back(S);
      }
      Sample s1 = combine_samples(S1);
      cout << "N=" << N << " s=" << speed << " |C|=" << s1.C
           << " |C|*s/log(N)=" << s1.C*speed/(log(N))
           << endl;
      cout << s1 << endl;
    }
  }
}
