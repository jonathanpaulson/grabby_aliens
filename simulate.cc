#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <functional>
#include <unordered_set>
using namespace std;
using ll = int64_t;
using ld = double;

std::minstd_rand RNG(std::random_device{}());
uniform_real_distribution<ld> DIST(0.0, 1.0);
ld r01() {
  return DIST(RNG);
}

struct Civ {
  Civ(ll D) : V(D, 0.0), T(0.0) {}

  static Civ mk_random(ll D, ld power, ld L) {
    Civ ret(D);
    for(ll j=0; j<D; j++) {
      ret.V[j] = r01()*L;
    }
    ret.T = pow(r01(), 1.0/(1.0+power));
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
ld distance2(const vector<ld>& A, const vector<ld>& B) {
  ld d2 = 0.0;
  for(ll i=0; i<A.size(); i++) {
    ld dx = abs(A[i]-B[i]);
    dx = min(dx, 1.0-dx);
    d2 += dx*dx;
  }
  return d2;
}
ld distance(const vector<ld>& A, const vector<ld>& B) {
  return sqrt(distance2(A,B));
}
ld sq(ld x) { return x*x; }

void simulate(ll D, ld speed, ld n, ll N, ld c, ld L) {
  vector<Civ> C;
  C.reserve(N);
  for(ll i=0; i<N; i++) {
    C.push_back(Civ::mk_random(D, n, L));
    if(i%1000000==0) {
      cerr << "GENERATE i=" << i << endl;
    }
  }
  sort(C.begin(), C.end(), [](Civ& A, Civ& B) { return A.T < B.T; });

  vector<Civ> ALIVE;
  for(ll i=0; i<N; i++) {
    Civ cand = C[i];
    bool is_alive = true;

    for(ll j=0; j<ALIVE.size(); j++) {
      Civ& alive = ALIVE[j];
      // dead iff
      // i.T > j.T + dij/speed
      // i.T*speed > j.T*speed + dij
      // i.T*speed - j.T*speed > dij
      // dij < speed*(i.T*j.T)
      // dij^2 < (speed*(i.T*j.T))^2
      ld d2 = distance2(alive.V, cand.V);
      bool dead = d2 < sq(speed*(cand.T-alive.T));
      if(dead) {
        is_alive = false;
        break;
      }
    }
    if(is_alive) {
      ALIVE.push_back(cand);
    }
    if(i%1000000==0) {
      cerr << "i=" << i << " |C|=" << ALIVE.size() << endl;
    }
  }
  assert(ALIVE.size() > 0);

  for(ll i=0; i<D; i++) {
    cout << static_cast<char>('X'+i) << ",";
  }
  cout << "OriginTime,MinWait,NumberSeen,MaxAngle" << endl;

  for(ll i=0; i<ALIVE.size(); i++) {
    auto c1 = ALIVE[i];
    for(ll j=0; j<ALIVE.size(); j++) {
      auto c2 = ALIVE[j];
      if(i!=j) {
        ld dij = distance(c1.V, c2.V);
        ld wij = (dij/speed - (c1.T-c2.T))/2.0;
        ld oij = c2.T + dij/c;
        ld bij = speed*(c1.T-oij)/dij;
        if(c1.T > oij) {
          c1.nsee++;
          c1.max_angle = max(c1.max_angle, bij);
        }
        c1.min_wait = min(c1.min_wait, wij);
      }
    }
    // max angle = max{j that i can see in C} bi,j.
    assert(c1.min_wait < 1e6);
    cout << c1 << endl;
  }
}

int main(int, char** argv) {
  ll D = stoll(argv[1]);
  ld n = atof(argv[2]);
  ll N = stoll(argv[3]);
  ld speed = 1.0/atof(argv[4]);
  ld c = stoll(argv[5]);
  ld L = atof(argv[6]);
  cerr << "D=" << D << " n=" << n << " N=" << N << " speed=" << speed << " c=" << c << " L=" << L << endl;
  simulate(D, speed, n, N, c, L);
}
