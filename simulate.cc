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
#include <tuple>
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
    return ans;
  }
  ll I;
  ld LnCurMax = 0.0;
};

struct Civ {
  Civ(ll D) : V(D, 0.0), T(0.0) {}

  static vector<ld> random_point(ll D, ld L) {
    vector<ld> R(D, 0.0);
    for(ll i=0; i<D; i++) {
      R[i] = r01()*L;
    }
    return R;
  }
  static Civ mk_random(ll D, ld power, ld L, SortedRNG& R) {
    Civ ret(D);
    ret.V = random_point(D, L);

    ld t = 1.0 - R.next();
    ret.T = pow(t, 1.0/(1.0+power));
    return ret;
  }
  vector<ld> V; // position in space
  ld T; // origin time
  ld min_arrival = 1e9; // min time when another civ arrives at our origin
  ld min_see = 1e9; // min time when we see signals from another civ
  ll nsee = 0; // number of other civs whose signals we see at our origin time
  ld max_angle = 0.0; // max "angle" among civs we see at our origin time
  ld percent_empty = 0.0; // how much of the universe is empty at our origin time
  ld volume = 0.0; // Fraction of the universe controlled by this civ at the end of time
};
ostream& operator<<(ostream& o, const Civ& C) {
  for(ll i=0; i<C.V.size(); i++) {
    o << C.V[i] << ",";
  }
  o << C.T << "," << C.min_arrival << "," << C.min_see << "," << C.nsee << "," << C.max_angle << "," << C.percent_empty << "," << C.volume;
  return o;
}

// Between every pair of points i,j is a Euclidean distance of dij = sqrt((xi-xj)^2 + (yi-yj)^2 + (zi-zj)^2)
// Use “hypertorus” distance metric, that identifies opposite sides of box [0,L]^D.
// i.e. instead of "dx=xi-xj", use "dx=abs(xi-xj); dx=min(dx,L-dx)"
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
    assert(NUM[i] <= NUM[i+1]);
  }

  using Element = pair<ld,pair<ll,ll>>;

  priority_queue<Element,vector<Element>,std::greater<Element>> Q;
  for(ll i=0; i<DEN.size(); i++) {
    ld value = 13.787/DEN[i] * NUM[0];
    Q.push(make_pair(value, make_pair(i, static_cast<ll>(0))));
  }
  ll k = Q.size();
  ll ai = 0;
  vector<ld> ANS;
  while(!Q.empty()) {
    Element x = Q.top(); Q.pop();
    ll di = x.second.first;
    ll ni = x.second.second;
    // Only sample every kth element, to save memory
    if(ai%k==0) {
      ANS.push_back(x.first);
    }
    ai++;

    if(ni+1 < NUM.size()) {
      ld value = 13.787/DEN[di] * NUM[ni+1];
      Q.push(make_pair(value, make_pair(di, ni+1)));
    } else if(Q.empty() && ai%k!=0) { // always include the max
      ANS.push_back(x.first);
    }
  }

  assert(NUM.size()<=ANS.size() && ANS.size()<=NUM.size()+1);
  for(ll i=0; i+1<ANS.size(); i++) {
    assert(ANS[i]<=ANS[i+1]);
  }
  return ANS;
}

vector<tuple<ld,ld,ld>> to_years(const vector<Civ>& C, ld m) {
  for(ll i=0; i+1<C.size(); i++) {
    assert(C[i].T < C[i+1].T);
  }
  vector<ld> N_ORIGIN;
  vector<ld> N_WAIT;
  vector<ld> N_SEE;
  vector<ld> DEN;
  for(ll i=0; i<C.size(); i++) {
    ld time_power = 1.0 / (1.0 - m);
    ld real_t = pow(C[i].T, time_power);
    ld real_arrival = pow(C[i].min_arrival, time_power);
    ld real_see = pow(C[i].min_see, time_power);

    ld wait = (real_arrival - real_t) / 2.0;
    ld wait_see = max(static_cast<ld>(0.0), real_see - real_t);

    N_ORIGIN.push_back(real_t);
    N_WAIT.push_back(wait);
    N_SEE.push_back(wait_see);
    if(C[i].nsee == 0) {
      DEN.push_back(real_t);
    }
  }
  
  sort(N_WAIT.begin(), N_WAIT.end());
  sort(N_SEE.begin(), N_SEE.end());
  vector<ld> ORIGIN_YEARS = ratio_distribution(N_ORIGIN, DEN);
  vector<ld> WAIT_YEARS = ratio_distribution(N_WAIT, DEN);
  vector<ld> SEE_YEARS = ratio_distribution(N_SEE, DEN);
  assert(ORIGIN_YEARS.size() == WAIT_YEARS.size());
  assert(ORIGIN_YEARS.size() == SEE_YEARS.size());
  vector<tuple<ld,ld,ld>> ANS;
  for(ll i=0; i<ORIGIN_YEARS.size(); i++) {
    ANS.push_back(make_tuple(ORIGIN_YEARS[i], WAIT_YEARS[i], SEE_YEARS[i]));
  }
  return ANS;
}

vector<Civ> simulate(ll D, ld speed, ld n, ll N, ld c, ld L, ll empty_samples, ll volume_samples) {
  /* This is a simpler way of generating the candidate civs.
     Instead, I generate them one-by-one by generating the origin times already in sorted
     order (see SortedRNG for more details on this).

  vector<Civ> C;
  C.reserve(N);
  for(ll i=0; i<N; i++) {
    C.push_back(Civ(D));
    C[i].V = Civ::random_point(D, L);
    C[i].T = r01();
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
      if(empty_samples == 0) {
        cand.percent_empty = 0.0;
      } else {
        ll nalive = 0;
        for(ll k=0; k<empty_samples; k++) {
          vector<ld> PT = Civ::random_point(D, L);
          bool pt_alive = true;
          for(ll j=0; j<ALIVE.size(); j++) {
            auto& alive = ALIVE[j];
            ld d2 = distance2(alive.V, PT, L);
            bool dead = d2 < sq(speed*(cand.T-alive.T));
            if(dead) {
              pt_alive = false;
              break;
            }
          }
          if(pt_alive) { nalive++; }
        }
        cand.percent_empty = static_cast<ld>(nalive)/static_cast<ld>(empty_samples);
      }

      cerr << "i=" << i << " |C|=" << ALIVE.size() << " percent_empty=" << cand.percent_empty << " n=" << n << endl;
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
        ld arrival = c2.T + dij/speed;
        ld oij = c2.T + dij/c;

        // t*(s+c) = d + s*t0 + c*t1
        ld see_time = (dij + speed*c1.T + c*c2.T) / (speed + c);
        c1.min_see = min(c1.min_see, see_time);

        if(c1.T > oij) {
          c1.nsee++;
          ld dt = abs(c1.T - c2.T);
          assert(dt > 0);
          ld angle_b = 1 + sq(speed/c);
          ld angle_a = (1.0 - sqrt(1.0 - angle_b*(1.0 - sq(dij/(c*dt)))))/angle_b;
          ld angle = 2*atan((speed/c)*(angle_a/(1-angle_a)));
          assert(0.0 < angle && angle < M_PI);
          c1.max_angle = max(c1.max_angle, angle);
        }
        c1.min_arrival = min(c1.min_arrival, arrival);
      }
    }
    assert(c1.min_arrival < 1e6);
    assert(c1.min_see < 1e6);
  }

  for(ll t=0; t<volume_samples; t++) {
    vector<ld> PT = Civ::random_point(D, L);
    pair<ld,ld> best = make_pair(1e9, 0);
    for(ll i=0; i<ALIVE.size(); i++) {
      auto& c1 = ALIVE[i];
      ld di = distance(c1.V, PT, L);
      ld ti = (c1.T + di/speed);
      if(ti < best.first) {
        best = make_pair(ti, i);
      }
    }
    ALIVE[best.second].volume += 1.0/volume_samples;
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
  ll empty_samples = stoll(argv[9]);
  ld m = atof(argv[10]);
  ll volume_samples = stoll(argv[11]);

  RNG.seed(seed);

  cerr << "D=" << D << " n=" << n << " N=" << N << " speed=" << speed << " c=" << c << " L=" << L << endl;

  vector<Civ> CIVS = simulate(D, speed, n, N, c, L, empty_samples, volume_samples);
  std::ofstream civ_out (fname+"_civs.csv", std::ofstream::out);
  for(ll i=0; i<D; i++) {
    civ_out << static_cast<char>('X'+i) << ",";
  }
  civ_out << "OriginTime,MinArrival,MinSee,NumberSeen,MaxAngle,PctEmpty,Volume" << endl;
  for(auto& civ : CIVS) {
    civ_out << civ << endl;
  }
  civ_out.close();

  vector<tuple<ld,ld,ld>> years = to_years(CIVS, m);
  std::ofstream year_out (fname+"_years.csv", std::ofstream::out);
  year_out << "OriginTime,MinWait,MinSETI" << endl;
  for(auto& y : years) {
    year_out << get<0>(y) << "," << get<1>(y) << "," << get<2>(y) << endl;
  }
  year_out.close();
}
