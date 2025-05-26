// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "G4NistManager.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"
#include "G4PVPlacement.hh"
#include "G4Transform3D.hh"
#include "G4UnitsTable.hh"
#include "G4Material.hh"
#include "G4Colour.hh"

#include "G4SubtractionSolid.hh"
// Geometries
#include "G4Box.hh"
#include "G4Trap.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"

#include "G4MultiFunctionalDetector.hh"
#include <G4VSensitiveDetector.hh>
#include <G4SDManager.hh>

#include "SensitiveVetoSD.hh"
#include "SensitiveSD.hh"

#include <sstream>

DetectorConstruction::DetectorConstruction() : G4VUserDetectorConstruction()
{
    fDetMess = new DetectorMessenger(this);
}

DetectorConstruction::~DetectorConstruction()
{
    ;
}

G4VPhysicalVolume *DetectorConstruction::Construct()
{
    G4bool checkOverlaps = true;

    G4NistManager *man = G4NistManager::Instance();

    G4Element *Li = man->FindOrBuildElement("Li");
    G4Element *N = man->FindOrBuildElement("N");
    G4Element *La = man->FindOrBuildElement("La");
    G4Element *Br = man->FindOrBuildElement("Br");
    G4Element *Ce = man->FindOrBuildElement("Ce");

    G4Material *Silicon = new G4Material("Silicon", 14., 28.0855 * g / mole, 2.33 * g / cm3);
    G4Material *VetoMat = Silicon;

    G4Element *H = man->FindOrBuildElement("H");
    G4Element *B = man->FindOrBuildElement("B");
    G4Element *C = man->FindOrBuildElement("C");
    G4Element *O = man->FindOrBuildElement("O");
    G4Element *F = man->FindOrBuildElement("F");
    G4Element *Na = man->FindOrBuildElement("Na");
    G4Element *Mg = man->FindOrBuildElement("Mg");
    G4Element *Al = man->FindOrBuildElement("Al");
    G4Element *Si = man->FindOrBuildElement("Si");
    G4Element *S = man->FindOrBuildElement("S");
    G4Element *Cl = man->FindOrBuildElement("Cl");
    G4Element *K = man->FindOrBuildElement("K");
    G4Element *Ca = man->FindOrBuildElement("Ca");
    G4Element *V = man->FindOrBuildElement("V");
    G4Element *Cr = man->FindOrBuildElement("Cr");
    G4Element *Fe = man->FindOrBuildElement("Fe");
    G4Element *Co = man->FindOrBuildElement("Co");
    G4Element *Ni = man->FindOrBuildElement("Ni");
    G4Element *Cu = man->FindOrBuildElement("Cu");
    G4Element *Zn = man->FindOrBuildElement("Zn");
    G4Element *As = man->FindOrBuildElement("As");
    G4Element *Sr = man->FindOrBuildElement("Sr");
    G4Element *Cd = man->FindOrBuildElement("Cd");
    G4Element *Ba = man->FindOrBuildElement("Ba");
    G4Element *Hg = man->FindOrBuildElement("Hg");
    G4Element *Pb = man->FindOrBuildElement("Pb");

    // Sand Sediment - 27 elemental composition
    G4Material *SandSediment = new G4Material("SandSediment", 1.535 * g / cm3, 27);
    SandSediment->AddElement(H, 3.30144 * perCent);
    SandSediment->AddElement(B, 0.00003 * perCent);
    SandSediment->AddElement(C, 0.35344 * perCent);
    SandSediment->AddElement(O, 61.9999 * perCent);
    SandSediment->AddElement(F, 0.00001 * perCent);
    SandSediment->AddElement(Na, 1.34864 * perCent);
    SandSediment->AddElement(Mg, 0.00734 * perCent);
    SandSediment->AddElement(Al, 2.59482 * perCent);
    SandSediment->AddElement(Si, 27.5365 * perCent);
    SandSediment->AddElement(S, 0.00517 * perCent);
    SandSediment->AddElement(Cl, 0.11055 * perCent);
    SandSediment->AddElement(K, 0.87891 * perCent);
    SandSediment->AddElement(Ca, 0.84806 * perCent);
    SandSediment->AddElement(V, 0.00100 * perCent);
    SandSediment->AddElement(Cr, 0.00100 * perCent);
    SandSediment->AddElement(Fe, 1.00237 * perCent);
    SandSediment->AddElement(Co, 0.00050 * perCent);
    SandSediment->AddElement(Ni, 0.00100 * perCent);
    SandSediment->AddElement(Cu, 0.00100 * perCent);
    SandSediment->AddElement(Zn, 0.00200 * perCent);
    SandSediment->AddElement(As, 0.00080 * perCent);
    SandSediment->AddElement(Br, 0.00038 * perCent);
    SandSediment->AddElement(Sr, 0.00104 * perCent);
    SandSediment->AddElement(Cd, 0.00010 * perCent);
    SandSediment->AddElement(Ba, 0.00250 * perCent);
    SandSediment->AddElement(Hg, 0.000002 * perCent);
    SandSediment->AddElement(Pb, 0.00150 * perCent);

    G4Material *vacuum = man->FindOrBuildMaterial("G4_Galactic");

    G4int ncomponents;
    G4double density, temperature, pressure;

    density = 0.001225 * g / cm3;
    pressure = 98658.96 * pascal;
    temperature = 273 * kelvin;

    G4Element *Ar = man->FindOrBuildElement("Ar");

    G4Material *Air = new G4Material("Air", density, ncomponents = 4, kStateGas, temperature, pressure);
    Air->AddElement(N, 75.53 * perCent);
    Air->AddElement(O, 23.18 * perCent);
    Air->AddElement(Ar, 1.28 * perCent);
    Air->AddElement(C, 0.01 * perCent);

    G4NistManager *nist_man = G4NistManager::Instance();
    G4Material *AIR_mat = nist_man->FindOrBuildMaterial("G4_AIR");

    G4Material *H2O = new G4Material("Water", 1.000 * g / cm3, 2);
    H2O->AddElement(H, 2);
    H2O->AddElement(O, 1);
    H2O->GetIonisation()->SetMeanExcitationEnergy(78.0 * eV);

    G4Material *H3BO3 = new G4Material("Boric_Acid", 1.44 * g / cm3, 3);
    H3BO3->AddElement(H, 3);
    H3BO3->AddElement(B, 1);
    H3BO3->AddElement(O, 3);
    // GetIonisation()->SetMeanExcitationEnergy(78.0*eV);

    G4Material *NaCl = new G4Material("Salt", 2.165 * g / cm3, 2);
    NaCl->AddElement(Na, 1);
    NaCl->AddElement(Cl, 1);

    G4Material *LiF = new G4Material("Lithium_Fluoride", 2.64 * g / cm3, 2);
    LiF->AddElement(Li, 1);
    LiF->AddElement(F, 1);

    G4Material *SeaWater = new G4Material("SeaWater", 1.027 * g / cm3, 2);
    SeaWater->AddMaterial(H2O, 98.8 * perCent);
    SeaWater->AddMaterial(NaCl, 1.2 * perCent);

    G4Material *LaBr3 = new G4Material("LaBr3", 5.08 * g / cm3, 2);
    LaBr3->AddElement(La, 1);
    LaBr3->AddElement(Br, 3);

    G4Material *LaBr3_Ce = new G4Material("LaBr3_Ce", 5.29 * g / cm3, 2);
    LaBr3_Ce->AddMaterial(LaBr3, 95 * perCent);
    LaBr3_Ce->AddElement(Ce, 5 * perCent);

    G4Material *C4H8Cl2S = new G4Material("Mustard_Gas", 1.27 * g / cm3, 4);
    C4H8Cl2S->AddElement(C, 4);
    C4H8Cl2S->AddElement(H, 8);
    C4H8Cl2S->AddElement(Cl, 2);
    C4H8Cl2S->AddElement(S, 1);

    G4Material *TNT = new G4Material("TNT", 1.654 * g / cm3, 4);
    TNT->AddElement(O, 42.26 * perCent);
    TNT->AddElement(C, 37.02 * perCent);
    TNT->AddElement(N, 18.50 * perCent);
    TNT->AddElement(H, 2.22 * perCent);

    G4Material *Clark1 = new G4Material("Clark1", 1.422 * g / cm3, 4);
    Clark1->AddElement(C, 12);
    Clark1->AddElement(H, 10);
    Clark1->AddElement(As, 1);
    Clark1->AddElement(Cl, 1);

    G4Material *Clark2 = new G4Material("Clark2", 1.45 * g / cm3, 4);
    Clark2->AddElement(C, 13);
    Clark2->AddElement(H, 10);
    Clark2->AddElement(As, 1);
    Clark2->AddElement(N, 1);

    G4Material *SiO2 = new G4Material("Silicon_Dioxide", 2.196 * g / cm3, 2);
    SiO2->AddElement(Si, 1);
    SiO2->AddElement(O, 2);

    G4Material *Iron = new G4Material("Iron", 7.850 * g / cm3, 1);
    Iron->AddElement(Fe, 1);

    G4Material *Lead = new G4Material("Lead", 11.34 * g / cm3, 1);
    Lead->AddElement(Pb, 1);

    G4Material *Paraffin = man->FindOrBuildMaterial("G4_PARAFFIN");

    G4Material *Aluminum = new G4Material("Al", 2.7 * g / cm3, 1);
    Aluminum->AddElement(Al, 1);

    G4Material *Paraffin_Boric_Acid = new G4Material("Paraffin_Boric_Acid", 0.972 * g / cm3, 2);
    Paraffin_Boric_Acid->AddMaterial(Paraffin, 90 * perCent);
    Paraffin_Boric_Acid->AddMaterial(H3BO3, 10 * perCent);

    G4Material *Cadmium = new G4Material("Cadmium", 8.65 * g / cm3, 1);
    Cadmium->AddElement(Cd, 1);

    G4Material *Carbon = new G4Material("Carbon", 2.3 * g / cm3, 1);
    Carbon->AddElement(C, 1);

    G4Material *Low_Carbon_Steel = new G4Material("Low_Carbon_Steel", 7.8499 * g / cm3, 2);
    Low_Carbon_Steel->AddMaterial(Iron, 99.75 * perCent);
    Low_Carbon_Steel->AddMaterial(Carbon, 0.25 * perCent);

    G4Material *LCSt = man->FindOrBuildMaterial("Low_Carbon_Steel");

    G4Material *polypropylene = man->FindOrBuildMaterial("G4_POLYPROPYLENE");

    G4Material *TargetMat = SeaWater;

    switch (targetType)
    {
    case TargetVariables::fWater:
        TargetMat = SeaWater;
        break;
    case TargetVariables::fMustardGas:
        TargetMat = C4H8Cl2S;
        break;
    case TargetVariables::fTNT:
        TargetMat = TNT;
        break;
    case TargetVariables::fClark1:
        TargetMat = Clark1;
        break;
    case TargetVariables::fClark2:
        TargetMat = Clark2;
        break;
    }

    

    
 

//--------------------geometry--------------------//
//---------------volumes-definition---------------//
    //World
    G4double worldSizeDimX = 4.2 * m;
    G4double worldSizeDimY = 5.2 * m;
    G4double worldSizeDimZ = 4.2 * m;

    G4Box *solidWorld = new G4Box("World", 0.5 * worldSizeDimX, 0.5 * worldSizeDimY, 0.5 * worldSizeDimZ);
    G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld, SeaWater, "World");
    G4VPhysicalVolume *physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0);

    G4VisAttributes *worldVisAtt = new G4VisAttributes(G4Colour(0., 0., 0.));
    // worldVisAtt->SetForceWireframe(true);
    logicWorld->SetVisAttributes(worldVisAtt);

    // SAND
    G4double sandXLength = 4 * m;
    G4double sandYLength = 100 * cm;
    G4double sandZLength = 4 * m;
    G4Box *solidSANDVolume = new G4Box("SAND", 0.5 * sandXLength, 0.5 * sandYLength, 0.5 * sandZLength);
    G4LogicalVolume *logicSANDVolume = new G4LogicalVolume(solidSANDVolume, SandSediment, "SAND");
    
    G4VisAttributes *sandVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    sandVisAtt->SetForceSolid(true);
    logicSANDVolume->SetVisAttributes(sandVisAtt);

    // Target
    G4double targetCoverDimX = 194 * cm;
    G4double targetCoverDimY = 50 * cm;
    G4double targetCoverDimZ = 50 * cm;
    G4Box *solidTGVolumeCover = new G4Box("TargetVolumeCover", 0.5 * targetCoverDimX, 0.5 * targetCoverDimY, 0.5 * targetCoverDimZ);
    G4LogicalVolume *logicTGVolumeCover = new G4LogicalVolume(solidTGVolumeCover, LCSt, "TargetVolumeCover");
    
    G4double targetWallSize = 3 * mm;
    G4Box *solidTGVolume = new G4Box("TargetVolume", 0.5 * targetCoverDimX - targetWallSize, 0.5 * targetCoverDimY - targetWallSize, 0.5 * targetCoverDimZ - targetWallSize);
    G4LogicalVolume *logicTGVolume = new G4LogicalVolume(solidTGVolume, TargetMat, "TargetVolume");
    
    G4VisAttributes *mustardVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0));
    mustardVisAtt->SetForceSolid(true);
    logicTGVolume->SetVisAttributes(mustardVisAtt);


    // Detector Case
    G4double detectorCaseRmin = 0 * cm;
    G4double detectorCaseRmax = 15. * cm;
    G4double detectorCaseHalfLength = 15 * cm;
    G4double detectorCaseWallSize = 3 * mm;
    
    G4VSolid *solidDetectorCaseOuter = new G4Tubs("DetectorCase", detectorCaseRmin, detectorCaseRmax, detectorCaseHalfLength, 0., 2 * M_PI * rad);
    G4LogicalVolume *logicDetectorCaseOuter = new G4LogicalVolume(solidDetectorCaseOuter, LCSt, "DetectorCase");

    G4VSolid *solidDetectorCaseInner = new G4Tubs("DetectorCaseInner", detectorCaseRmin, detectorCaseRmax - detectorCaseWallSize, detectorCaseHalfLength - detectorCaseWallSize, 0., 2 * M_PI * rad);
    G4LogicalVolume *logicDetectorCaseInner = new G4LogicalVolume(solidDetectorCaseInner, Air, "DetectorCaseInner");

    G4VisAttributes *caseAtt = new G4VisAttributes(G4Color(0, 0., 0.));
    caseAtt->SetForceWireframe(true);
    logicDetectorCaseOuter->SetVisAttributes(caseAtt);

    
    // Detector:
    G4double detectorRadiusMin = 0;
    G4double detectorRadiusMax = 2.54 * cm;
    G4double detectorHalfLengthZ = 2.54 * cm;
    G4VSolid *solidAbsor = new G4Tubs("DetectorLaBr", detectorRadiusMin, detectorRadiusMax, detectorHalfLengthZ, 0., 2 * M_PI * rad);
    G4LogicalVolume *logicAbsor = new G4LogicalVolume(solidAbsor, LaBr3_Ce, "DetectorLaBr");
    
    G4VisAttributes *solidVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
    solidVisAtt->SetForceSolid(true);
    logicAbsor->SetVisAttributes(solidVisAtt);

    //Generator Case
    G4double generatorCaseRmin = 0 * cm;
    G4double generatorCaseRmax = 15. * cm;
    G4double generatorCaseHalfLength = 15 * cm;
    G4double generatorCaseWallSize = 3 * mm;
    
    G4VSolid *solidGeneratorCaseOuter = new G4Tubs("GeneratorCase", generatorCaseRmin, generatorCaseRmax, generatorCaseHalfLength, 0., 2 * M_PI * rad);
    G4LogicalVolume *logicGeneratorCaseOuter = new G4LogicalVolume(solidGeneratorCaseOuter, LCSt, "GeneratorCase");

    G4VSolid *solidGeneratorCaseInner = new G4Tubs("GeneratorCaseInner", generatorCaseRmin, generatorCaseRmax - generatorCaseWallSize, generatorCaseHalfLength - generatorCaseWallSize, 0., 2 * M_PI * rad);
    G4LogicalVolume *logicGeneratorCaseInner = new G4LogicalVolume(solidGeneratorCaseInner, Air, "GeneratorCaseInner");

    logicGeneratorCaseOuter->SetVisAttributes(caseAtt);

    //Generator
    G4double generatorDimX = 10.8 * cm;
    G4double generatorDimY = 20.3 * cm;
    G4double generatorDimZ = 25.4 * cm;
    G4double generatorWallSize = 3 * mm;
    G4VSolid *solidGeneratorOuter = new G4Box("Generator", 0.5 * generatorDimX, 0.5 * generatorDimY, 0.5 * generatorDimZ);
    G4LogicalVolume *logicGeneratorOuter = new G4LogicalVolume(solidGeneratorOuter, LCSt, "Generator");
    G4VSolid *solidGeneratorInner = new G4Box("GeneratorInner", 0.5 * generatorDimX - generatorWallSize, 0.5 * generatorDimY - generatorWallSize, 0.5 * generatorDimZ - generatorWallSize);
    G4LogicalVolume *logicGeneratorInner = new G4LogicalVolume(solidGeneratorInner, vacuum, "GeneratorInner");

    G4VisAttributes *genVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 0.0));
    genVisAtt->SetForceWireframe(true);
    logicGeneratorOuter->SetVisAttributes(genVisAtt);

    //WROV Arm
    G4double armRmin = 0 * cm;
    G4double armRmax = 12.5/2 * cm;
    G4double armHalfLengthZ = 68.6 / 2 * cm;
    G4VSolid *solidArm = new G4Tubs("Arm", armRmin, armRmax, armHalfLengthZ, 0., 2 * M_PI * rad);
    G4LogicalVolume *logicArm = new G4LogicalVolume(solidArm, LCSt, "Arm");

    G4VisAttributes *armAtt = new G4VisAttributes(G4Color(1.0, 1., 1.0));
    armAtt->SetForceSolid(true);
    logicArm->SetVisAttributes(armAtt);

    // Source border for presen.
    G4double pRmin = 0.*cm;
    G4double pRmax = 2.5*cm;
    G4double pSPhi = 0*deg;
    G4double pDPhi = 2*M_PI;
    G4double pSTheta = 0*deg;
    G4double pDTheta = M_PI;
    G4Sphere* solidSourceBorder = new G4Sphere("SOURCE", pRmin, pRmax, pSPhi, pDPhi, pSTheta, pDTheta);
    G4LogicalVolume* logicSourceBorder = new G4LogicalVolume(solidSourceBorder, vacuum, "SOURCE");
    
    G4VisAttributes * sourceVisAtt = new G4VisAttributes(G4Colour(1.,1.,0.));
    sourceVisAtt->SetForceWireframe(true);
    logicSourceBorder->SetVisAttributes(sourceVisAtt);

    //ROV 1 cm of polypropylene(PP) and 3 mm of steel, seawater inside
    G4double ROVDimX = 215 * cm; // Length
    G4double ROVDimY = 125.4 * cm; //Height 
    G4double ROVDimZ = 116 * cm; //Width
    G4double ROVWallThicknessPP = 1 * cm;
    G4double ROVWallThicknessSteel = 3 * mm;

    G4VSolid *solidROVPP = new G4Box("ROVPP", 0.5 * ROVDimX, 0.5 * ROVDimY, 0.5 * ROVDimZ);
    G4LogicalVolume *logicROVPP = new G4LogicalVolume(solidROVPP, polypropylene, "ROVPP");
    G4VSolid *solidROVSteel = new G4Box("ROVSteel", 0.5 * ROVDimX - ROVWallThicknessPP, 0.5 * ROVDimY - ROVWallThicknessPP, 0.5 * ROVDimZ - ROVWallThicknessPP);
    G4LogicalVolume *logicROVSteel = new G4LogicalVolume(solidROVSteel, LCSt, "ROVSteel");
    G4VSolid *solidROVInner = new G4Box("ROVSteelInner", 0.5 * ROVDimX - ROVWallThicknessPP - ROVWallThicknessSteel, 0.5 * ROVDimY - ROVWallThicknessPP - ROVWallThicknessSteel, 0.5 * ROVDimZ - ROVWallThicknessPP - ROVWallThicknessSteel);
    G4LogicalVolume *logicROVInner = new G4LogicalVolume(solidROVInner, SeaWater, "ROVSteelInner");

    G4VisAttributes *ROVVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    ROVVisAtt->SetForceWireframe(true);
    logicROVPP->SetVisAttributes(ROVVisAtt);




    // Veto construction
    // Active area of silicon detector -> 5 - 55 mm (circle radius)
    G4double vetoDimx = 8 * cm;
    G4double vetoDimY = 1 * mm; // 2cm -> max thickness for silicon detectors
    G4double vetoDimZ = 8 * cm;
    
    G4Box *solidVeto = new G4Box("DetectorSi", 0.5 * vetoDimx, 0.5 * vetoDimY, 0.5 * vetoDimZ);
    G4LogicalVolume *logicVeto = new G4LogicalVolume(solidVeto, VetoMat, "DetectorSi");
    

    G4VisAttributes *VetoAtt = new G4VisAttributes(G4Color(1.0, 1., 0.5));
    VetoAtt->SetForceSolid(true);
    logicVeto->SetVisAttributes(VetoAtt);

    


//------------Placement-----------//
    //Sand
    G4ThreeVector sandShift(0., -170 * cm, 0.);
    new G4PVPlacement(0, sandShift, logicSANDVolume, "SAND", logicWorld, false, 1, checkOverlaps);
    //target
    G4RotationMatrix *rotationTarget = new G4RotationMatrix();
    rotationTarget->rotateY(90 * deg);
    G4ThreeVector targetCoverShift = sandShift +G4ThreeVector(
        0.,
        sandYLength / 2 + targetCoverDimY / 2,
        0.);
    new G4PVPlacement(rotationTarget, targetCoverShift, logicTGVolumeCover, "TargetVolumeCover", logicWorld, false, 2, checkOverlaps);
    new G4PVPlacement(0, G4ThreeVector(), logicTGVolume, "TargetVolume", logicTGVolumeCover, false, 3, checkOverlaps);

    //construction parameters
    G4double detectorGeneratorDistance = 10 * cm;

    //Detector case
    G4RotationMatrix *rotDet = new G4RotationMatrix();
    rotDet->rotateX(90*deg);
    G4ThreeVector detectorCaseShift = targetCoverShift + G4ThreeVector(
        0.,
        targetCoverDimY / 2 + detectorCaseHalfLength,
        -detectorCaseRmax-detectorGeneratorDistance/2);
    new G4PVPlacement(rotDet, detectorCaseShift, logicDetectorCaseOuter, "DetectorCase", logicWorld, false, 4, checkOverlaps);
    new G4PVPlacement(0, G4ThreeVector(), logicDetectorCaseInner, "DetectorCaseInner", logicDetectorCaseOuter, false, 5, checkOverlaps);
    //Detector
    new G4PVPlacement(0, G4ThreeVector(), logicAbsor, "DetectorLaBr", logicDetectorCaseInner, false, 6, checkOverlaps); 

    //Generator case
    G4RotationMatrix *rotGen = new G4RotationMatrix();
    rotGen->rotateX(90*deg);
    G4ThreeVector generatorCaseShift = targetCoverShift + G4ThreeVector(
        0.,
        targetCoverDimY / 2 + generatorCaseHalfLength,
        generatorCaseRmax+detectorGeneratorDistance/2);
    new G4PVPlacement(rotGen, generatorCaseShift, logicGeneratorCaseOuter, "GeneratorCase", logicWorld, false, 7, checkOverlaps);
    new G4PVPlacement(0, G4ThreeVector(), logicGeneratorCaseInner, "GeneratorCaseInner", logicGeneratorCaseOuter, false, 8, checkOverlaps);
    //Generator
    new G4PVPlacement(0, G4ThreeVector(), logicGeneratorOuter, "Generator", logicGeneratorCaseInner, false, 9, checkOverlaps);
    new G4PVPlacement(0, G4ThreeVector(), logicGeneratorInner, "GeneratorInner", logicGeneratorOuter, false, 10, checkOverlaps);
    //Source border for presentation
    new G4PVPlacement(0, G4ThreeVector(0, 0 ,0), logicSourceBorder, "SOURCE", logicGeneratorInner, false, 20);
    //Veto  (inside generator)
    G4RotationMatrix *rotVeto = new G4RotationMatrix();
    rotVeto->rotateX(-90*deg);
    G4ThreeVector vetoShift(0., 0., 5.0 * cm );
    new G4PVPlacement(rotVeto, vetoShift, logicVeto, "DetectorSi", logicGeneratorInner, false, 21, checkOverlaps);

    std::cout << "--------------------------Source should be placed at: " << generatorCaseShift  << std::endl;
    
    //Arm - detector
    G4ThreeVector armShiftDet = detectorCaseShift + G4ThreeVector(
        0.,
        detectorCaseHalfLength + armHalfLengthZ,
        0.);
    new G4PVPlacement(rotDet, armShiftDet, logicArm, "Arm", logicWorld, false, 11, checkOverlaps);
    //Arm - generator
    G4ThreeVector armShiftGen = generatorCaseShift + G4ThreeVector(
        0.,
        generatorCaseHalfLength + armHalfLengthZ,
        0.);
    new G4PVPlacement(rotGen, armShiftGen, logicArm, "Arm", logicWorld, false, 12, checkOverlaps);

    //ROV
    G4RotationMatrix *rotROV = new G4RotationMatrix();
    rotROV->rotateZ(90*deg);
    G4ThreeVector ROVShift = (armShiftGen + armShiftDet)/2 + G4ThreeVector(
        0.,
        armHalfLengthZ + ROVDimX / 2,
        0.);
    new G4PVPlacement (rotROV, ROVShift, logicROVPP, "ROVPP", logicWorld, false, 14, checkOverlaps);
    new G4PVPlacement (0, G4ThreeVector(), logicROVSteel, "ROVSteel", logicROVPP, false, 15, checkOverlaps);
    new G4PVPlacement (0, G4ThreeVector(), logicROVInner, "ROVSteelInner", logicROVSteel, false, 16, checkOverlaps);


    
    
    
    return physWorld;
}

void DetectorConstruction::ConstructSDandField()
{
    G4SDManager *sdManager = G4SDManager::GetSDMpointer();
    sdManager->SetVerboseLevel(2);

    SensitiveSD *myDetector = new SensitiveSD("Detector");
    SetSensitiveDetector("DetectorLaBr", myDetector);
    sdManager->AddNewDetector(myDetector);

    SensitiveVetoSD *vetoDetector = new SensitiveVetoSD("Veto");
    SetSensitiveDetector("DetectorSi", vetoDetector);
    sdManager->AddNewDetector(vetoDetector);
}
