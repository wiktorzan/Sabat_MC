#ifndef FILE_ANALYZER_H
#define FILE_ANALYZER_H

#include "HistCollection.h"
#include "EnergySmearer.h"
#include <string>

class FileAnalyzer {
    private:
        HistCollection& histo;
        EnergySmearer smearer;
        int GetMassNumberFromProducts(std::string label);
        std::vector<std::string> SplitString(std::string& s, const std::string& delimiter);
    public:
        FileAnalyzer(HistCollection& hist, EnergySmearer smearer);
        void Analyze(const std::string& filename);
};

#endif