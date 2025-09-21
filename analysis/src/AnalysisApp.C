#include "AnalysisApp.h"
#include <TSystemDirectory.h>
#include <TSystemFile.h>

AnalysisApp::AnalysisApp() : smearer(2.0*pow(10, -4), 2.22*pow(10, -2),  0.5)
{
    Double_t minEnergy = 0, maxEnergy = 12, binSize = 0.01; // in MeV
    Double_t minTime = 0, timeSeparator = 1, maxTime = 1500, binSizeSmall = 0.001, binSizeLarge = 1; // in us

    hist.CreateEnergyHistos(minEnergy, maxEnergy, binSize);
    hist.CreateTimingHistos(minTime, timeSeparator, maxTime, binSizeSmall, binSizeLarge);
    hist.CreateEnergyVsTimingHistos(minEnergy, maxEnergy, binSize, minTime, maxTime, binSizeLarge);  

    //Smearing parameters
    // double a = 2.0*pow(10, -4); // in MeV
    // double b = 2.22*pow(10, -2);
    // double c = 0.5;
    // this->smearer = EnergySmearer(a, b, c);
}

void AnalysisApp::RunAnalysis(std::vector<std::string> filenames) {

    FileAnalyzer analyzer(hist, smearer);

    if(filenames.size() > 1) 
    {
        for (const auto& file : filenames) 
        {
            analyzer.Analyze(file);
        }
    } else if (filenames.size() == 1) 
    {
        std::vector<std::string> filesFromPattern = AnalyzePattern(filenames[0]);
        for (const auto& file : filesFromPattern) 
        {
            analyzer.Analyze(file);
        }
    } else {
        std::cout << "No files to analyze." << std::endl;
        return;
    }



    hist.SaveHistos("Results.root");
}

std::vector<std::string> AnalysisApp::AnalyzePattern(std::string fileOrPattern)
{
    std::vector<std::string> filesToAnalyze;
    TString root_file = fileOrPattern;
    TString star = "*";
    TString slash = "/";
    std::size_t starPlace = fileOrPattern.rfind(star);
    std::size_t slashPlace = fileOrPattern.rfind(slash);

    if(starPlace == std::string::npos) {
        filesToAnalyze.push_back(fileOrPattern);
        return filesToAnalyze;
    } else if (root_file[strlen(root_file) - 1] == star) {
        std::string pattern;
        std::string dirName;
        if (slashPlace > 0) {
            dirName = fileOrPattern.substr(0, slashPlace+1);
            pattern = fileOrPattern.substr(slashPlace+1, starPlace-slashPlace-1);
        }
        else {
            dirName = "";
            pattern = fileOrPattern.substr(0, starPlace);
        }
        std::cout << "Getting files from directory: " << dirName << " and pattern " << pattern << std::endl;
        TString dirNameRoot = dirName;
        if (dirNameRoot == "")
        dirNameRoot = ".";
        TSystemDirectory dir(dirNameRoot, dirNameRoot);
        TList *files = dir.GetListOfFiles();
    
        if (files) {
        TSystemFile *file;
        TString fname;
        TIter next(files);
        while ((file = (TSystemFile*)next())) {
            fname = file->GetName();
    
            if (!file->IsDirectory() && fname.EndsWith(".root") && fname.BeginsWith(pattern)) {
            std::cout << "Adding file to the analysis queue " << fname.Data() << std::endl;
            filesToAnalyze.push_back(dirName + fname.Data());
            }
        }
        }
    }    
    return filesToAnalyze;
}
