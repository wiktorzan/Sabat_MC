// #include "AnalysisApp.h"
#include "AnalysisApp.h"
#include <iostream>

int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cout << "Not enough arguments. Up to two needed: name of the file | name of the second file ..." << std::endl;
        return 0;
      }

      std::vector<std::string> filenames;
      for (int i = 1; i < argc; ++i) {
          filenames.push_back(argv[i]);
      }

 
      AnalysisApp app;
      app.RunAnalysis(filenames);
    
    return 0;
}


