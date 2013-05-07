
#include <string>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <limits>
#include <climits>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <cassert>

using namespace std;

#define rep(i,n) for(int i = 0; i < (n); i++)
#define rap(i,n,m) for(int i = (n); i < (m); i++)
#define iter(it, v) for(typeof((v).begin()) it = (v).begin(); it != (v).end(); ++it)
#define MP(x,y) make_pair(x,y)
#define F first
#define S second
typedef unsigned int uint;
typedef unsigned long long ull;
typedef long long ll;
typedef pair<int,int> pii;
typedef pair<double,double> pdd;

const double Inf = numeric_limits<double>::infinity();

int sumi(int *a, int n) { int sum = 0; rep(i,n) sum += a[i]; return sum; }
double sumd(double *a, int n) { double sum = 0; rep(i,n) sum += a[i]; return sum; }
int maxi(int *a, int n) { int m = INT_MIN; rep(i,n) if(a[i] > m) m = a[i]; return m; }
double maxd(double *a, int n) { double m = -1e300; rep(i,n) if(a[i] > m) m = a[i]; return m; }
int mini(int *a, int n) { int m = INT_MAX; rep(i,n) if(a[i] < m) m = a[i]; return m; }
double mind(double *a, int n) { double m = 1e300; rep(i,n) if(a[i] < m) m = a[i]; return m; }

//some settings
const unsigned batchsizemin = 4, //minimum nr of words shown at once
               batchsizemax = 50; //maximum nr of words shown at once

double target_dunno = 5; //there should be approximately how many forgotten words among the ones that are shown at once?

string trainfilename, auxfilename;
uint emptylines, lookahead;
double target;

struct stamp {
  stamp(int t, int k) : time(t), know(k) {}
  int time, know;
  };

struct data {
  data(const vector<stamp> &s, char l) : stamps(s), lang(l) {}
  data(char l) : lang(l) {}
  data() {}
  vector<stamp> stamps;
  char lang;
  };

map<string, data> db;
map<char, string> dbs;

pii get_lastkd(const vector<stamp> &v) {
  int lastknow = 0, lastdunno = 0;
  rep(i, v.size()) (v[i].know ? lastknow : lastdunno) = v[i].time;
  return pii(lastknow, lastdunno);
  }
pii get_nkd(const vector<stamp> &v) {
  int nknow = 0, ndunno = 0;
  rep(i, v.size()) (v[i].know ? nknow : ndunno)++;
  return pii(nknow, ndunno); 
  }
pii get_maxkt(const vector<stamp> &v) {
  int maxkt = 0, max2kt = 0; //longest and second longest time without forgetting
  rap(i, 1, v.size()) {
    int kt = v[i].know*(v[i].time - v[i-1].time);
    if(kt > maxkt) max2kt = maxkt, maxkt = kt;
    else if(kt > max2kt) max2kt = kt;
    }
  return pii(maxkt, max2kt);
  }
pii get_sumkdt(const vector<stamp> &v) {
  int sumk = 0, sumd = 0;
  rap(i, 1, v.size()) (v[i].know ? sumk : sumd) += v[i].time - v[i-1].time;
  return pii(sumk, sumd);
  }
pdd get_sumlkdt(const vector<stamp> &v) {
  double sumk = 0, sumd = 0;
  rap(i, 1, v.size()) (v[i].know ? sumk : sumd) += log(1 + v[i].time - v[i-1].time);
  return pdd(sumk, sumd);
  }

struct params {
  int firstseen, knew, lastknow, lastdunno, nknow, ndunno, maxkt, max2kt, sumkt, sumdt, k[10];
  double sumlkt, sumldt;
  time_t now;
  };

void getparams(const vector<stamp> &v, time_t now, params *p) {
  if(v.empty()) {
    cerr << "getparams called with empty v\n";
    exit(1);
    return;
    }
  p->firstseen = v[0].time;
  p->knew = v[v.size()-1].know;
  pii lastkd = get_lastkd(v);
  p->lastknow = lastkd.F, p->lastdunno = lastkd.S;
  pii nkd = get_nkd(v);
  p->nknow = nkd.F, p->ndunno = nkd.S;
  pii maxkt = get_maxkt(v);
  p->maxkt = maxkt.F, p->max2kt = maxkt.S;
  pii sumkdt = get_sumkdt(v);
  p->sumkt = sumkdt.F, p->sumdt = sumkdt.S;
  rep(i, 10) p->k[i] = int(v.size()) > i ? v[v.size()-1-i].know*1 : 1;
  pdd sumlkdt = get_sumlkdt(v);
  p->sumlkt = sumlkdt.F, p->sumldt = sumlkdt.S;
  p->now = now;
  }

double ilogit(double x) {
  double t = exp(x);
  return t/(1+t);
  }

const int npredstot = 30;
int npreds[4], predinds[4][npredstot];
const char *prednames[npredstot] = {
	"(Intercept)",		"lsl", 			"I(la*lsl)", 		"la",
	"I(la^2*lsl)",		"I(la^2)", 		"I(lmaxkt*lsl)", 	"lmaxkt",
	"I(lmax2kt*lsl)", 	"lmax2kt", 		"I(d*lsl)", 		"d",
	"I(lsumkt*lsl)", 	"lsumkt",		"I(lsumdt*lsl)",	"lsumdt",
	"I(sumlkt*lsl)", 	"sumlkt", 		"I(sumldt*lsl)",	"sumldt",
	"I(lsumldt*lsl)", 	"lsumldt", 		"I(c*lsl)",		"c",
	"I((k2+k3+k4)*lsl)", 	"I(k2+k3+k4)", 		"I((k5+k6+k7)*lsl)", 	"I(k5+k6+k7)",
	"I(lnseen*lsl)", 	"lnseen"};
double coefs[4][npredstot];

double rate(const vector<stamp> v) {
  if(v.empty()) return 0;
  params p;
  getparams(v, time(0)+lookahead, &p);
  int knew=p.knew,
      sincelast=p.now-max(p.lastknow,p.lastdunno),
      a=p.lastknow - (p.lastdunno ? p.lastdunno : p.firstseen-400),
      maxkt = p.maxkt,
      nseen = p.nknow + p.ndunno,
      nknow = p.nknow,
      ndunno = p.ndunno;
  double lsl = log(1+sincelast),
         la = log(1+a),
         lmaxkt = log(1+maxkt),
         lmax2kt = log(1+p.max2kt),
         c = double(nknow)/nseen,
	 d = nknow/(1.+ndunno),
         lsumkt = log(1+p.sumkt),
         lsumdt = log(1+p.sumdt),
         sumlkt = p.sumlkt,
         sumldt = p.sumldt,
         lsumldt = log(1+sumldt),
         lnseen = log(nseen),
         k234 = sumi(p.k+1, 3),
         k567 = sumi(p.k+4, 3);
  double preds[] = {
	1,			lsl,			la*lsl,			la,
	la*la*lsl,		la*la,			lmaxkt*lsl,		lmaxkt,
	lmax2kt*lsl,		lmax2kt,		d*lsl,			d,
	lsumkt*lsl,		lsumkt,			lsumdt*lsl,		lsumdt,
	sumlkt*lsl,		sumlkt,			sumldt*lsl,		sumldt,
	lsumldt*lsl,		lsumldt,		c*lsl,			c,
	k234*lsl,		k234,			k567*lsl,		k567,
	lnseen*lsl,		lnseen};
  int ind = 2*!knew + !!maxkt;
  double sum = 0;
  rep(i, npreds[ind])
    sum += coefs[ind][i]*preds[predinds[ind][i]];
  return ilogit(sum);
  }

istream &getline2(istream &is, string &s) {
  getline(is, s);
  while(s[s.size()-1] == '\\') {
    string ss;
    getline(is, ss);
    s += '\n' + ss;
    }
  return is;
  }

void start() {
  if(target_dunno <= 0) {
    ofstream out(trainfilename.c_str());
    return;
    }
  vector<pair<pair<double, double>, string> > cand;
  cand.reserve(db.size());
  iter(it, db) {
    double rating = rate(it->S.stamps);
    int lastupdate = it->S.stamps.size() ? it->S.stamps.back().time : 0;
    double sortby;
    if(rating > target)
      sortby = rating;
    else
      sortby = -lastupdate;
    cand.push_back(make_pair(make_pair(sortby, rating), it->F));
    }
  double totaldunno = 0;
  int undertarget = 0;
  iter(it, cand) {
    totaldunno += 1-fabs(it->F.S);
    if(it->F.S < target)
      undertarget++;
    }
  sort(cand.begin(), cand.end());
  uint numcand = 0;
  for(double e_dunno = 0;
      e_dunno < target_dunno && numcand < db.size();
      numcand++)
    e_dunno += 1 - fabs(cand[numcand].F.S);
  if(numcand < batchsizemin)
    numcand = min(batchsizemin, uint(db.size()));
  else if(numcand > batchsizemax)
    numcand = batchsizemax;
  if(!numcand || db.empty()) {
    cout << "database is empty\n";
    return;
    }
  cout << "count: " << numcand
    << ", worst: " << cand[0].F.S
    << ", best: " << cand[numcand-1].F.S
    << ", totdunno: " << totaldunno
    << ", <target: " << undertarget << endl;
  srand(time(0));

  for(unsigned i = 0; i < numcand; i++) {
    int j = i + rand()%(numcand-i);
    swap(cand[i], cand[j]);
    }
  ofstream out(trainfilename.c_str());
  if(!out.good()) {
    cout << "could not open file:\n" << trainfilename <<'\n';
    return;
    }
  for(unsigned i = 0; i < emptylines; i++)
    out << endl;
  for(map<char, string>::iterator it0 = dbs.begin(); it0 != dbs.end(); it0++) {
    char ch = it0->first;
    bool first = 1;
    for(unsigned i = 0; i < numcand; i++)
      if(db[cand[i].S].lang == ch) {
        if(first)
          out << it0->S << ":\n";
        first = 0;
        out << cand[i].S << endl;
        }
    }
  }

bool finish(string filename) {
  ifstream in(filename.c_str());
  if(!in.good()) {
    cout << "could not open file:\n" << filename <<'\n';
    return 0;
    }
  string s;
  vector<string> tochange;
  vector<pair<string, int> > input;
  bool hasoldornew = 0, isnew = 1, haspipe = 0;
  char ch;
  while(getline2(in, s))
    if(!s.empty()) {
      if(s[s.size()-1] == ':' || s[0] == '#')
        continue;
      if(s == "old" || (s.size() == 5 && s.substr(0, 4) == "new " && dbs.count(ch = s[4]))) {
        if(hasoldornew) {
          cout << "old or new can only be specified once!\n"
            << "nothing done.\n";
          return 0;
          }
        hasoldornew = 1;
        if(s == "old")
          isnew = 0;
        continue;
        }
      bool dunno = false;
      size_t pos = s.find_first_of('*');
      if(pos != s.npos) {
        s.erase(pos, 1);
        dunno = true;
        while((pos = s.find_first_of('*')) != s.npos)
          s.erase(pos, 1);
        }
      pos = s.find_first_of('#');
      if(pos != s.npos) {
        s.erase(pos, 1);
        while((pos = s.find_first_of('#')) != s.npos)
          s.erase(pos, 1);
        tochange.push_back(s);
        }
      pos = s.find_first_of('|');
      if(pos != s.npos) haspipe = 1;
      input.push_back(make_pair(s, !dunno));
      }
  if(haspipe) {
    cout << "entries cannot contain \'|\'.\n"
      << "nothing done.\n";
    return 0;
    }
  if(!input.empty() && !hasoldornew) {
    cout << "you need to let one line contain only \"old\" or \"new <lang>\"!\n"
      << "nothing done.\n";
    return 0;
    }
  for(unsigned i = 0; i < input.size(); i++) {
    if(isnew && db.count(input[i].first)) {
      cout << input[i].first << endl
        << "already existed, but you specified \"new\"!\n"
        << "nothing done.\n";
      return 0;
      }
    if(!isnew && !db.count(input[i].first)) {
      cout << input[i].first << endl
        << "didn't exist, but you specified \"old\"!\n"
        << "nothing done.\n";
      return 0;
      }
    }
  if(isnew && !tochange.empty()) {
    cout << "you used \"#\", but you specified \"new\"!\n"
      << "nothing done.\n";
    return 0;
    }
  for(unsigned i = 0; i < input.size(); i++) {
    if(isnew)
      db[input[i].first] = data(ch);
    db[input[i].first].stamps.push_back(stamp(time(0), input[i].second));
    }

  if(!tochange.empty()) {
    ofstream out(filename.c_str());
    for(unsigned i = 0; i < tochange.size(); i++)
      out << tochange[i] << endl;
    out.close();
    cout << "press enter when you're done changing\n";
    cin.ignore(1024, '\n');
    ifstream in(filename.c_str());
    for(unsigned i = 0; getline2(in, s) && i < tochange.size(); i++) {
      if(s != tochange[i]) {
        if(!s.empty())
          db[s] = db[tochange[i]], cout << "copied to " << s.substr(0, 5) << "...\n";
        db.erase(tochange[i]), cout << "deleting " << tochange[i].substr(0, 5) << "...\n";
        }
      }
    }
  ofstream out(filename.c_str());
  return 1;
  }

void stat() {
  vector<pair<char, string> > dbs2;
  dbs2.push_back(MP('*', string("all languages")));
  iter(it,dbs)
    dbs2.push_back(*it);
  iter(it0, dbs2) {
    int total = 0, last24h = 0, last48h = 0, last7d = 0, last30d = 0;
    int now = time(0);
    char ch = it0->first;
    iter(it, db) {
      if(it->second.lang != ch && ch != '*')
        continue;
      int time = (it->second.stamps.empty() ? now : it->second.stamps[0].time);
      total++;
      if(time > now - 24*60*60)
        last24h++;
      if(time > now - 48*60*60)
        last48h++;
      if(time > now - 7*24*60*60)
        last7d++;
      if(time > now - 30*24*60*60)
        last30d++;
      }
    if(it0 != dbs2.begin())
      cout << endl;
    cout << it0->second << ":\n"
      << "total:\t\t" << total << endl
      << "last 24h\t" << last24h << endl
      << "last 48h\t" << last48h << endl
      << "last 7d\t\t" << last7d << endl
      << "last 30d\t" << last30d << endl;
    }
  }

void undo() {
  time_t last = 0;
  for(map<string, data>::iterator it = db.begin(); it != db.end(); it++)
    if(it->second.stamps.back().time > last)
      last = it->second.stamps.back().time;
  int num = 0;
  for(map<string, data>::iterator it = db.begin(); it != db.end(); it++)
    if(it->second.stamps.back().time == last)
      num++, cout << it->first << endl;
  cout << "restore " << num << " entries, ok?\n";
  string s;
  getline(cin, s);
  if(s[0] == 'y') {
    vector<string> todelete;
    for(map<string, data>::iterator it = db.begin(); it != db.end(); it++)
      if(it->second.stamps.back().time == last) {
        it->second.stamps.pop_back();
        if(it->second.stamps.empty())
          todelete.push_back(it->first);
        }
    for(unsigned i = 0; i < todelete.size(); i++)
      db.erase(todelete[i]);
    }
  }

void loadsettings() {
  ifstream settingsfile("settings");
  if(!settingsfile.good()) {
    cout << "could not open file \"settings\"\n";
    return;
    }
  getline(settingsfile, trainfilename);
  getline(settingsfile, auxfilename);
  settingsfile >> emptylines >> target;
  }

void savesettings() {
  ofstream settingsfile("settings");
  if(!settingsfile.good()) {
    cout << "could not open file \"settings\"\n";
    return;
    }
  settingsfile << trainfilename << endl
    << auxfilename << endl
    << emptylines << endl
    << target << endl;
  }

void printraw() {
  ofstream out("raw");
  out << "knew know firstseen lastknow lastdunno nknow ndunno maxkt max2kt sumkt sumdt sumlkt sumldt now k1 k2 k3 k4 k5 k6 k7 k8 k9 k10\n";
  iter(it, db) {
    params p;
    vector<stamp> v = it->S.stamps;
    if(v.empty()) continue;
    time_t now = v.back().time;
    int know = v.back().know;
    v.pop_back();
    while(v.size() > 0) {
      getparams(v, now, &p);
      out<<p.knew<<' '<<know<<' '<<p.firstseen<<' '<<p.lastknow<<' '
         <<p.lastdunno<<' '<<p.nknow<<' '<<p.ndunno<<' '<<p.maxkt<<' '
         <<p.max2kt<<' '<<p.sumkt<<' '<<p.sumdt<<' '<<p.sumlkt<<' '
         <<p.sumldt<<' '<<p.now<<' ';
      rep(i,10) out<<p.k[i]<<' ';
      out<<'\n';
      now = v.back().time;
      know = v.back().know;
      v.pop_back();
      }
    }
  }


int getpredind(string s) {
  rep(i, npredstot) if(s == prednames[i]) return i;
  return -1;
  }

void loadcoefs() {
  ifstream in("coefs");
  string tmp;
  rep(i, 4) {
    in >> tmp;
    in >> npreds[i];
    if(npreds[i] > npredstot)
      cerr << "too many predictors\n", exit(1);
    rep(j, npreds[i]) {
      in >> tmp, predinds[i][j] = getpredind(tmp);
      if(predinds[i][j] < 0) cerr<<"unknown predictor: "<<tmp<<endl, exit(1);
      }
    rep(j, npreds[i])
      in >> coefs[i][j];
    }
  }

void calibrate() {
  printraw();
  system("./cal.R > coefs");
  loadcoefs();
  }

//database format:
//string|int:1/0|int:1/0|...|0\n

int main(int argc, char **argv) {
  if(argc < 2) {
    printf("Usage: %s <command>\n"
"available commands: stat[istics], c[ontinue], aux[iliary], cal[ibrate], undo\n", argv[0]);
    return 1;
    }
  string command = argv[1];
  loadsettings();
  loadcoefs();
  ifstream dbsfile("databases");
  char ch;
  string st;
  while(dbsfile >> ch >> st) {
    if(st.empty()) {
      cout << "syntax error in file \"databases\", aborting\n";
      return 1;
      }
    dbs[ch] = st;
    ifstream dbfile(st.c_str());
    if(!dbfile.good()) {
      cout << "could not open file:\n" << st <<'\n';
      return 1;
      }
    string s;
    while(getline(dbfile, s, '|')) {
      vector<stamp> vi;
      int t, b;
      dbfile >> t;
      while(t != 0) {
        dbfile.ignore();
        dbfile >> b;
        dbfile.ignore();
        vi.push_back(stamp(t, b));
        dbfile >> t;
        }
      dbfile.ignore();
      db[s] = data(vi, ch);
      }
    dbfile.close();
    }

//(does a start with b?)
#define starts(a, b) ((b) == string(a).substr(0, b.size()))
  if(starts("statistics", command))
    {stat(); return 0;}
  else if(starts("continue", command)) {
    if(argc >= 3)
      target_dunno = atof(argv[2]);
    if(argc >= 4)
      lookahead = 3600*atof(argv[3]);
    if(finish(trainfilename))
      start();
    }
  else if(starts("finish", command))
    finish(trainfilename);
  else if(starts("auxiliary", command))
    finish(auxfilename);
  else if(starts("calibrate", command))
    {calibrate(); return 0;}
  else if(starts("undo", command))
    undo();
  else if(starts("printraw", command))
    printraw(), cout << "printed to file \'raw\'\n";
  else {
    cout << "invalid command\n";
    return 1;
    }

  for(map<char, string>::iterator it0 = dbs.begin(); it0 != dbs.end(); it0++) {
    ofstream dbout((it0->second+".new").c_str(), ios_base::trunc);
    if(!dbout.good()) {
      cout << "could not open file:\n" << it0->second+".new" <<'\n';
      return 1;
      }
    for(map<string, data>::iterator it = db.begin(); it != db.end(); it++)
      if(it0->first == it->second.lang) {
        dbout << it->first << '|';
        for(unsigned i = 0; i < it->second.stamps.size(); i++)
          dbout << it->second.stamps[i].time << ':' << it->second.stamps[i].know << '|';
        dbout << "0\n";
        }
    dbout.close();
    system(("mv " + it0->second + ".new " + it0->second).c_str());
    }
  }
