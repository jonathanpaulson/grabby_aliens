#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <functional>
#include <unordered_set>
#include <cfloat>
#include <iomanip>
#include <fstream>
#include <queue>
using namespace std;
using ll = int64_t;
using ld = long double;

std::minstd_rand RNG(0);
uniform_real_distribution<ld> DIST(0.0, 1.0);
ld r01() {
  return DIST(RNG);
}

// https://apps.dtic.mil/dtic/tr/fulltext/u2/a066739.pdf
// Samples from a sorted list of N U(0,1) variables
struct SortedRNG {
  SortedRNG(ll N) : I(N), LnCurMax(0.0) {}
  ld next() {
    LnCurMax += log(r01())/I;
    I--;
    ld ans = exp(LnCurMax);
    //cerr << setprecision(30) << LnCurMax << " " << ans << endl;
    return ans;
  }
  ll I;
  ld LnCurMax = 0.0;
};

struct Civ {
  Civ(ll D) : V(D, 0.0), T(0.0) {}

  static Civ mk_random(ll D, ld power, ld L, SortedRNG& R) {
    Civ ret(D);
    for(ll j=0; j<D; j++) {
      ret.V[j] = r01()*L;
    }

    ld t = 1.0 - R.next();
    ret.T = pow(t, 1.0/(1.0+power));
    return ret;
  }
  vector<ld> V; // position in space
  ld T; // origin time
  ld min_wait = 1e9; // min "wait time" until we see another civ
  ll nsee = 0; // number of other civs whose signals we see at our origin time
  ld max_angle = 0.0; // max "angle" among civs we see at our origin time
};
ostream& operator<<(ostream& o, const Civ& C) {
  for(ll i=0; i<C.V.size(); i++) {
    o << C.V[i] << ",";
  }
  o << C.T << "," << C.min_wait << "," << C.nsee << "," << C.max_angle;
  return o;
}

// Between every pair of points i,j is a Euclidean distance of di,j = ((xi-xj)2 + (yi-yj)2 + (zi-zj)2)1/2.
// Use “hypertorus” distance metric, that identifies opposite sides of box [0,1]d. (Can calc via, instead of dx = xi-xj, use (let a = abs(xi-xj), if a < ½, use a, else use 1-a).
ld distance2(const vector<ld>& A, const vector<ld>& B, ld L) {
  ld d2 = 0.0;
  for(ll i=0; i<A.size(); i++) {
    ld dx = abs(A[i]-B[i]);
    dx = min(dx, L-dx);
    d2 += dx*dx;
  }
  return d2;
}
ld distance(const vector<ld>& A, const vector<ld>& B, ld L) {
  return sqrt(distance2(A,B,L));
}
ld sq(ld x) { return x*x; }

vector<ld> ratio_distribution(const vector<ld>& NUM, const vector<ld>& DEN) {
  assert(NUM.size() > 0);
  assert(DEN.size() > 0);
  for(ll i=0; i+1<DEN.size(); i++) {
    assert(DEN[i] < DEN[i+1]);
  }
  for(ll i=0; i+1<NUM.size(); i++) {
    assert(NUM[i] < NUM[i+1]);
  }

  using Element = pair<ld,pair<ll,ll>>;

  priority_queue<Element,vector<Element>,std::greater<Element>> Q;
  for(ll i=0; i<DEN.size(); i++) {
    ld value = 13.8/DEN[i] * NUM[0];
    Q.push(make_pair(value, make_pair(i, static_cast<ll>(0))));
  }
  ll k = Q.size();
  ll ai = 0;
  vector<ld> ANS;
  while(!Q.empty()) {
    Element x = Q.top(); Q.pop();
    ll di = x.second.first;
    ll ni = x.second.second;
    if(ai%k==0) {
      //cerr << "DEBUG: " << x.first << " " << NUM[ni] << " " << DEN[di] << " di=" << di << " DEN.size=" << DEN.size() << endl;
      ANS.push_back(x.first);
    }
    ai++;

    if(ni+1 < NUM.size()) {
      ld value = 13.8/DEN[di] * NUM[ni+1];
      Q.push(make_pair(value, make_pair(di, ni+1)));
    } else if(Q.empty() && ai%k!=0) { // always include the max
      ANS.push_back(x.first);
    }
  }

  assert(NUM.size()<=ANS.size() && ANS.size()<=NUM.size()+1);
  for(ll i=0; i+1<ANS.size(); i++) {
    assert(ANS[i]<ANS[i+1]);
  }
  return ANS;
}

vector<pair<ld,ld>> to_years(const vector<Civ>& C) {
  for(ll i=0; i+1<C.size(); i++) {
    assert(C[i].T < C[i+1].T);
  }
  vector<ld> N_ORIGIN;
  vector<ld> N_WAIT;
  vector<ld> DEN;
  for(ll i=0; i<C.size(); i++) {
    N_ORIGIN.push_back(C[i].T);
    N_WAIT.push_back(C[i].min_wait);
    if(C[i].nsee == 0) {
      DEN.push_back(C[i].T);
    }
  }
  
  sort(N_WAIT.begin(), N_WAIT.end());
  vector<ld> ORIGIN_YEARS = ratio_distribution(N_ORIGIN, DEN);
  vector<ld> WAIT_YEARS = ratio_distribution(N_WAIT, DEN);
  assert(ORIGIN_YEARS.size() == WAIT_YEARS.size());
  vector<pair<ld,ld>> ANS;
  for(ll i=0; i<ORIGIN_YEARS.size(); i++) {
    ANS.push_back(make_pair(ORIGIN_YEARS[i], WAIT_YEARS[i]));
  }
  return ANS;
}

vector<Civ> simulate(ll D, ld speed, ld n, ll N, ld c, ld L) {
  /*
  vector<Civ> C;
  C.reserve(N);
  for(ll i=0; i<N; i++) {
    C.push_back(Civ::mk_random(D, n, L));
  }
  sort(C.begin(), C.end(), [](Civ& A, Civ& B) { return A.T < B.T; });
  */

  SortedRNG R(N);
  ll last_alive = 0;
  vector<Civ> ALIVE;
  for(ll i=0; i<N; i++) {
    Civ cand = Civ::mk_random(D, n, L, R);
    bool is_alive = true;

    for(ll j=0; j<ALIVE.size(); j++) {
      Civ& alive = ALIVE[j];
      assert(cand.T > alive.T);
      // dead iff
      // i.T > j.T + dij/speed
      // i.T*speed > j.T*speed + dij
      // i.T*speed - j.T*speed > dij
      // dij < speed*(i.T*j.T)
      // dij^2 < (speed*(i.T*j.T))^2
      ld d2 = distance2(alive.V, cand.V,L);
      bool dead = d2 < sq(speed*(cand.T-alive.T));
      if(dead) {
        is_alive = false;
        break;
      }
    }
    if(is_alive) {
      cerr << "i=" << i << " |C|=" << ALIVE.size() << endl;
      ALIVE.push_back(cand);
      last_alive = i;
    }
    if(i > last_alive + 1000000) { break; } // probably no more survivors
  }
  assert(ALIVE.size() > 0);
  cerr << "last_alive=" << last_alive << endl;

  for(ll i=0; i<ALIVE.size(); i++) {
    auto& c1 = ALIVE[i];
    for(ll j=0; j<ALIVE.size(); j++) {
      auto c2 = ALIVE[j];
      if(i!=j) {
        ld dij = distance(c1.V, c2.V,L);
        ld wij = (dij/speed - (c1.T-c2.T))/2.0;
        ld oij = c2.T + dij/c;
        ld bij = 2*speed*(c1.T-oij)/dij;
        assert(bij < 2.0);
        if(c1.T > oij) {
          c1.nsee++;
          c1.max_angle = max(c1.max_angle, bij);
        }
        c1.min_wait = min(c1.min_wait, wij);
      }
    }
    // max angle = max{j that i can see in C} bi,j.
    assert(c1.min_wait < 1e6);
  }
  return ALIVE;
}

int main(int, char** argv) {
  ll D = stoll(argv[1]);
  ld n = atof(argv[2]);
  ll N = stoll(argv[3]);
  ld speed = atof(argv[4]);
  ld c = atof(argv[5]);
  ld L = atof(argv[6]);
  string fname = argv[7];
  ll seed = stoll(argv[8]);

  RNG.seed(seed);

  cerr << "D=" << D << " n=" << n << " N=" << N << " speed=" << speed << " c=" << c << " L=" << L << endl;

  vector<Civ> CIVS = simulate(D, speed, n, N, c, L);
  std::ofstream civ_out (fname+".csv", std::ofstream::out);
  for(ll i=0; i<D; i++) {
    civ_out << static_cast<char>('X'+i) << ",";
  }
  civ_out << "OriginTime,MinWait,NumberSeen,MaxAngle" << endl;
  for(auto& civ : CIVS) {
    civ_out << civ << endl;
  }
  civ_out.close();

  vector<pair<ld,ld>> years = to_years(CIVS);
  std::ofstream year_out (fname+"_years.txt", std::ofstream::out);
  year_out << "OriginTime,MinWait" << endl;
  for(auto& y : years) {
    year_out << y.first << "," << y.second << endl;
  }
  year_out.close();
}
