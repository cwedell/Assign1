#include "Assign1.h"
#include <iostream>
#include <fstream>
#include <iterator>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <ctime>
#define PI 3.1415926535

using namespace std;

Assign1::Assign1() {
}

Assign1::~Assign1() {
}

int Assign1::countLines;
int Assign1::countBigrams;
int Assign1::sum;
float Assign1::mean;
float Assign1::variance;
float Assign1::sumDiffs;
float Assign1::sdev;
int Assign1::countA;
int Assign1::countC;
int Assign1::countT;
int Assign1::countG;
int Assign1::countAA;
int Assign1::countAC;
int Assign1::countAT;
int Assign1::countAG;
int Assign1::countCA;
int Assign1::countCC;
int Assign1::countCT;
int Assign1::countCG;
int Assign1::countTA;
int Assign1::countTC;
int Assign1::countTT;
int Assign1::countTG;
int Assign1::countGA;
int Assign1::countGC;
int Assign1::countGT;
int Assign1::countGG;
float Assign1::probA;
float Assign1::probC;
float Assign1::probT;
float Assign1::probG;
float Assign1::probAA;
float Assign1::probAC;
float Assign1::probAT;
float Assign1::probAG;
float Assign1::probCA;
float Assign1::probCC;
float Assign1::probCT;
float Assign1::probCG;
float Assign1::probTA;
float Assign1::probTC;
float Assign1::probTT;
float Assign1::probTG;
float Assign1::probGA;
float Assign1::probGC;
float Assign1::probGT;
float Assign1::probGG;

void Assign1::stats(std::ifstream& filename) { // to calculate statistics on file
  if(filename.is_open()) {
    string line;
    string cleanup;
    string compare;
    while(getline(filename, line)) {
      cleanup = line.c_str();
      compare = "";
      for(int i = 0; i < cleanup.length(); ++i) {
        compare = toupper(cleanup[i]);
        if(compare != "A" && compare != "C" && compare != "T" && compare != "G" && line.find_first_of(compare) != std::string::npos) {
          line.erase(line.find_first_of(compare), 1); // clean up original string
        }
      }
      for(int m = 0; m < line.length(); ++m) { // counting nucleotides
        string letter = "";
        letter = toupper(line[m]); // get current letter
        if(letter == "A") {
          ++Assign1::countA;
        } else if(letter == "C") {
          ++Assign1::countC;
        } else if(letter == "T") {
          ++Assign1::countT;
        } else {
          ++Assign1::countG;
        }
        if(m != 0) { // counting bigrams
          string lastLetter = "";
          lastLetter = toupper(line[m - 1]); // get previous letter
          if(lastLetter == "A") {
            if(letter == "A") {
              ++Assign1::countAA;
            } else if(letter == "C") {
              ++Assign1::countAC;
            } else if(letter == "T") {
              ++Assign1::countAT;
            } else {
              ++Assign1::countAG;
            }
          } else if(lastLetter == "C") {
            if(letter == "A") {
              ++Assign1::countCA;
            } else if(letter == "C") {
              ++Assign1::countCC;
            } else if(letter == "T") {
              ++Assign1::countCT;
            } else {
              ++Assign1::countCG;
            }
          } else if(lastLetter == "T") {
            if(letter == "A") {
              ++Assign1::countTA;
            } else if(letter == "C") {
              ++Assign1::countTC;
            } else if(letter == "T") {
              ++Assign1::countTT;
            } else {
              ++Assign1::countTG;
            }
          } else {
            if(letter == "A") {
              ++Assign1::countGA;
            } else if(letter == "C") {
              ++Assign1::countGC;
            } else if(letter == "T") {
              ++Assign1::countGT;
            } else {
              ++Assign1::countGG;
            }
          }
        }
      }
      Assign1::sum += line.length();
      ++Assign1::countLines;
      Assign1::countBigrams += line.length() - 1;
    }
    Assign1::mean = Assign1::sum / Assign1::countLines;
    filename.clear();
    filename.seekg(0, ios::beg);
    while(getline(filename, line)) {
      Assign1::sumDiffs += pow((line.length() - Assign1::mean), 2); // calculating differences for variance
    }
    Assign1::variance = sqrt(Assign1::sumDiffs) / (Assign1::countLines - 1);
    Assign1::sdev = sqrt(Assign1::variance);
    Assign1::probA = (float)Assign1::countA / (float)Assign1::sum;
    Assign1::probC = (float)Assign1::countC / (float)Assign1::sum;
    Assign1::probT = (float)Assign1::countT / (float)Assign1::sum;
    Assign1::probG = (float)Assign1::countG / (float)Assign1::sum;
    Assign1::probAA = (float)Assign1::countAA / (float)Assign1::countBigrams;
    Assign1::probAC = (float)Assign1::countAC / (float)Assign1::countBigrams;
    Assign1::probAT = (float)Assign1::countAT / (float)Assign1::countBigrams;
    Assign1::probAG = (float)Assign1::countAG / (float)Assign1::countBigrams;
    Assign1::probCA = (float)Assign1::countCA / (float)Assign1::countBigrams;
    Assign1::probCC = (float)Assign1::countCC / (float)Assign1::countBigrams;
    Assign1::probCT = (float)Assign1::countCT / (float)Assign1::countBigrams;
    Assign1::probCG = (float)Assign1::countCG / (float)Assign1::countBigrams;
    Assign1::probTA = (float)Assign1::countTA / (float)Assign1::countBigrams;
    Assign1::probTC = (float)Assign1::countTC / (float)Assign1::countBigrams;
    Assign1::probTT = (float)Assign1::countTT / (float)Assign1::countBigrams;
    Assign1::probTG = (float)Assign1::countTG / (float)Assign1::countBigrams;
    Assign1::probGA = (float)Assign1::countGA / (float)Assign1::countBigrams;
    Assign1::probGC = (float)Assign1::countGC / (float)Assign1::countBigrams;
    Assign1::probGT = (float)Assign1::countGT / (float)Assign1::countBigrams;
    Assign1::probGG = (float)Assign1::countGG / (float)Assign1::countBigrams;
  }
  filename.close();
}

void Assign1::printStats(std::ofstream& filename) { // print statistics to new output file
  if(filename.is_open()) {
    filename << "--------------------" << endl;
    filename << "Sum: " << Assign1::sum << endl;
    filename << "Mean: " << Assign1::mean << endl;
    filename << "Variance: " << Assign1::variance << endl;
    filename << "Standard deviation: " << Assign1::sdev << endl << endl;
    filename << "Probabilities:" << endl;
    filename << "\tA: " << Assign1::probA << endl;
    filename << "\tC: " << Assign1::probC << endl;
    filename << "\tT: " << Assign1::probT << endl;
    filename << "\tG: " << Assign1::probG << endl << endl;
    filename << "\tAA: " << Assign1::probAA << endl;
    filename << "\tAC: " << Assign1::probAC << endl;
    filename << "\tAT: " << Assign1::probAT << endl;
    filename << "\tAG: " << Assign1::probAG << endl;
    filename << "\tCA: " << Assign1::probCA << endl;
    filename << "\tCC: " << Assign1::probCC << endl;
    filename << "\tCT: " << Assign1::probCT << endl;
    filename << "\tCG: " << Assign1::probCG << endl;
    filename << "\tTA: " << Assign1::probTA << endl;
    filename << "\tTC: " << Assign1::probTC << endl;
    filename << "\tTT: " << Assign1::probTT << endl;
    filename << "\tTG: " << Assign1::probTG << endl;
    filename << "\tGA: " << Assign1::probGA << endl;
    filename << "\tGC: " << Assign1::probGC << endl;
    filename << "\tGT: " << Assign1::probGT << endl;
    filename << "\tGG: " << Assign1::probGG << endl << endl;
  }
}

void Assign1::gauss(std::ofstream& filename) { // create Gaussian distribution of new nucleotides
  double randA = 0.0;
  double randB = 0.0;
  double randC = 0.0;
  double randD = 0.0;
  string nukeString;
  float randChar = 0.0f;
  string toAdd;
  if(filename.is_open()) {
    filename << "Random DNA strings" << endl;
    filename << "--------------------" << endl;
    srand(time(NULL)); // set random seed
    for(int i = 0; i < 1000; ++i) {
      randA = (double)rand() / RAND_MAX;
      randB = (double)rand() / RAND_MAX;
      randC = sqrt(-2 * log(randA)) * cos(2 * PI * randB);
      randD = abs((Assign1::sdev * randC) + Assign1::mean);
      nukeString = ""; // initialize string for nucleotide
      for(int j = 0; j <= ceil(randD); ++j) {
        randChar = (double)rand() / RAND_MAX; // choose random character based on probabilities
        if(randChar <= Assign1::probA) {
          toAdd = "A";
        } else if(randChar <= Assign1::probA + Assign1::probC) {
          toAdd = "C";
        } else if(randChar <= Assign1::probA + Assign1::probC + Assign1::probT) {
          toAdd = "T";
        } else {
          toAdd = "G";
        }
        nukeString += toAdd; // add character to string
      }
      filename << nukeString << endl;
    }
  }
}

int main() {
  ifstream filein;
  bool fileValid = false;
  string fileToOpen = "nuketest.txt";
  while(!fileValid) { // check for valid filename
    try {
      filein.open(fileToOpen.c_str());
      fileValid = true;
    } catch(const std::exception& e) {
      cout << "File not found!" << endl;
      cout << "Enter the file name:" << endl;
      cin >> fileToOpen;
    }
  }
  ofstream fileout;
  fileout.open("ColtonWedell.out");
  fileout << "Colton Wedell" << endl;
  fileout << "Assignment 1" << endl;
  fileout << "9/14/2019" << endl;
  fileout << "CPSC350-03" << endl;
  fileout << "Rene German" << endl << endl;

  bool userFinished = false;

  Assign1 a;

  while(!userFinished) {
    a.stats(filein);
    cout << "Calculation complete." << endl;
    a.printStats(fileout);
    a.gauss(fileout);
    cout << "Gaussian distribution complete." << endl;
    string answer;
    cout << "Analyze another file? (y/n)" << endl;
    cin >> answer;
    bool doneAnalyzing = false;
    bool newFileValid = false;
    while(!doneAnalyzing) {
      if(answer == "y" || answer == "Y") {
        string newFileToOpen;
        while(!newFileValid) { // check for new valid filename
          cout << "Enter file name:" << endl;
          cin >> newFileToOpen;
          filein.open(newFileToOpen.c_str());
          if(filein.good()) {
            newFileValid = true;
            doneAnalyzing = true;
          } else {
            cout << "File not found." << endl;
            filein.close();
          }
        }
      } else if(answer == "n" || answer == "N") {
          cout << "Goodbye!" << endl;
          doneAnalyzing = true;
          userFinished = true;
      } else {
          cout << "Sorry, I didn't get that. Analyze another file? (y/n)" << endl;
          cin >> answer;
      }
    }
  }
  fileout.close();
  return 0;
}
