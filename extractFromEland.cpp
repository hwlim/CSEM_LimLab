#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<fstream>
#include<iostream>
#include<string>
using namespace std;

string line;
ifstream fin;
ofstream fout;

long long num_filt;

int main(int argc, char* argv[]) {
  if (argc != 3) {
    printf("Usage : extractFromEland inpF outF\n");
    exit(-1);
  }

  size_t pos0, pos1, pos2;
  string tag, seq, tp;

  long long cnt = 0;

  fin.open(argv[1]);
  fout.open(argv[2]);

  num_filt = 0;

  while (getline(fin, line)) {
    ++cnt;
    pos0 = line.find_first_of('\t');
    tag = line.substr(0, pos0);
    pos1 = line.find_first_of('\t', pos0 + 1);
    seq = line.substr(pos0 + 1, pos1 - pos0 - 1);
    pos2 = line.find_first_of('\t', pos1 + 1);
    tp = line.substr(pos1 + 1, pos2 - pos1 - 1);

    if (cnt % 1000000 == 0) cout<< cnt<< " FIN"<< endl;

    if (tp.compare("QC") == 0 || tp.compare("RM") == 0) {
      ++num_filt; continue;
    }

    fout<<tag<<endl<<seq<<endl;
  }

  fin.close();
  fout.close();

  cout<< "Number of Filtered Reads: "<< num_filt<< endl;

  return 0;
}
