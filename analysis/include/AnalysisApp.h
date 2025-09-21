#ifndef ANALYSIS_APP_H
#define ANALYSIS_APP_H

#include "FileAnalyzer.h"
#include "HistCollection.h"

class AnalysisApp {
    private:
    HistCollection hist;
    EnergySmearer smearer;
    public:
        AnalysisApp();
        void RunAnalysis(std::vector<std::string> filenames);
        std::vector<std::string> AnalyzePattern(std::string fileOrPattern);
};


#endif