#include <iostream>

using namespace std;

class Assign1 {
public:
  Assign1();
  ~Assign1();
  void stats(std::ifstream& filename);
  void printStats(std::ofstream& filename);
  void gauss(std::ofstream& filename);
  static int countLines;
  static int countBigrams;
  static int sum;
  static float mean;
  static float variance;
  static float sumDiffs;
  static float sdev;
  static int countA;
  static int countC;
  static int countT;
  static int countG;
  static int countAA;
  static int countAC;
  static int countAT;
  static int countAG;
  static int countCA;
  static int countCC;
  static int countCT;
  static int countCG;
  static int countTA;
  static int countTC;
  static int countTT;
  static int countTG;
  static int countGA;
  static int countGC;
  static int countGT;
  static int countGG;
  static float probA;
  static float probC;
  static float probT;
  static float probG;
  static float probAA;
  static float probAC;
  static float probAT;
  static float probAG;
  static float probCA;
  static float probCC;
  static float probCT;
  static float probCG;
  static float probTA;
  static float probTC;
  static float probTT;
  static float probTG;
  static float probGA;
  static float probGC;
  static float probGT;
  static float probGG;
};
