
#include <string>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <limits>
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

const double Inf = numeric_limits<double>::infinity();

//some settings
const unsigned batchsizemin = 8, //minimum words shown at once
               batchsizemax = 50; //maximum words shown at once

double target_dunno = 5; //there should be approximately how many forgotten words among the ones that are shown at once?

const double ncoefs[4] = {2, 10, 3, 6};
double coefs[4][10];

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

struct params {
  int firstseen, knew, lastknow, lastdunno, b, nseen, nknow;
  time_t now;
  };

void getparams(const vector<stamp> &v, time_t now, params *p) {
  int firstseen;
  if(v.empty()) {
    firstseen = 0;
    return;
    }
  firstseen = v[0].time;
  int prev = 0, lastdunno = 0, lastknow = 0,
      b = 0, // longest time without forgetting, but is decreased once we forget twice in a row
      nseen = 0, nknow = 0;
  iter(it, v) {
    if(prev > 0 && it->know && it->time - prev > b)
      b = it->time - prev;
    else if(prev > 0 && !it->know && it->time - prev < b)
      b = (b + it->time - prev)/2;
    (it->know ? lastknow : lastdunno) = it->time;
    nseen++;
    if(it->know) nknow++;
    prev = it->time;
    }
  *p = {firstseen, v.back().know, lastknow, lastdunno, b, nseen, nknow, now};
  }

double ilogit(double x) {
  double t = exp(x);
  return t/(1+t);
  }

double rate(const vector<stamp> v) {
  params p;
  getparams(v, time(0)+lookahead, &p);
  
  int knew=p.knew,
      sincelast=p.now-max(p.lastknow,p.lastdunno),
      a=p.lastknow-p.lastdunno-(1-!!p.lastdunno)*(p.firstseen-400),
      b=p.b,
      nseen = p.nseen,
      nknow = p.nknow;
  double lsl = log(1+sincelast),
         la = log(1+a),
         lb = log(1+b),
         c = double(nknow)/nseen;
  double params[4][10] = {{1,lsl,0,0,0,0,0,0,0,0},
                          {1,lsl,la*lsl,la*la*lsl,lb*lsl,c*lsl,la,la*la,lb,c},
                          {1,lsl,c,0,0,0,0,0,0,0},
                          {1,lsl,lb*lsl,c*lsl,lb,c,0,0,0,0}};
  if(!p.firstseen) return 0;
  int ind = 2*!knew + !!b;
  double sum = 0;
  rep(i, ncoefs[ind])
    sum += coefs[ind][i]*params[ind][i];
  return ilogit(sum);
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
    int lastupdate = it->S.stamps.back().time;
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
  assert(out.good() && !db.empty());
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
  assert(in.good());
  string s;
  vector<string> tochange;
  vector<pair<string, int> > input;
  bool hasoldornew = 0, isnew = 1;
  char ch;
  while(getline(in, s))
    if(!s.empty()) {
      if(s[s.size()-1] == ':')
        continue;
      if(s == "old" || (s.size() == 5 && s.substr(0, 4) == "new " && dbs.count(ch = s[4]))) {
        if(hasoldornew) {
          cout << "old or new can only be specified once!\n"
            << "nothing doe.\n";
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
      input.push_back(make_pair(s, !dunno));
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
    for(unsigned i = 0; getline(in, s) && i < tochange.size(); i++) {
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
  assert(settingsfile.good());
  getline(settingsfile, trainfilename);
  getline(settingsfile, auxfilename);
  settingsfile >> emptylines >> target;
  }

void savesettings() {
  ofstream settingsfile("settings");
  assert(settingsfile.good());
  settingsfile << trainfilename << endl
    << auxfilename << endl
    << emptylines << endl
    << target << endl;
  }

void printraw() {
  ofstream out("raw");
  out << "knew know firstseen lastknow lastdunno b nseen nknow now\n";
  iter(it, db) {
    params p;
    vector<stamp> v = it->S.stamps;
    if(v.empty()) continue;
    time_t now = v.back().time;
    int know = v.back().know;
    v.pop_back();
    while(v.size() > 0) {
      getparams(v, now, &p);
      out<<p.knew<<' '<<know<<' '<<p.firstseen<<' '<<p.lastknow<<' '<<p.lastdunno<<' '<<p.b<<' '<<p.nseen<<' '<<p.nknow<<' '<<p.now<<'\n';
      now = v.back().time;
      know = v.back().know;
      v.pop_back();
      }
    }
  }

void loadcoefs() {
  ifstream in("coefs");
  rep(i, 4) {
    in.ignore(100000, '\n');
    in.ignore(100000, '\n');
    rep(j, ncoefs[i])
      in >> coefs[i][j];
    in.ignore(100000, '\n');
    }
  }

void calibrate() {
  printraw();
  system("./cal.R > coeefs");
  loadcoefs();
  }

//database format:
//string|int:1/0|int:1/0|...|0\n

int main(int argc, char **argv) {
  assert(argc >= 2);
  string command = argv[1];
  loadsettings();
  loadcoefs();
  ifstream dbsfile("databases");
  char ch;
  string st;
  while(dbsfile >> ch >> st) {
    assert(!st.empty());
    dbs[ch] = st;
    ifstream dbfile(st.c_str());
    assert(dbfile.good());
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
  //  {analyze(); return 0;}
  else if(starts("continue", command)) {
    if(argc >= 3)
      target_dunno = atof(argv[2]);
    if(argc >= 4)
      lookahead = 3600*atof(argv[3]);
    if(finish(trainfilename))
      start();
    }
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
    assert(dbout.good());
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
