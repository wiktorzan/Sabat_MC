/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class (Mandatory)

#ifndef DETECTOR_CONSTRUCTION_HH
#define DETECTOR_CONSTRUCTION_HH

#include "DetectorMessenger.hh"
#include <G4VUserDetectorConstruction.hh>
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"

#include<string>

class G4VSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class DetectorMessenger;

enum TargetVariables {
  fWater, fMustardGas, fTNT, fClark1, fClark2
};

/// Detector construction class to define materials (with their physical properties) and detector geometry.
class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  /// constructor
  DetectorConstruction();
  /// destructor
  virtual ~DetectorConstruction();

  G4VPhysicalVolume* Construct() override;
  
  void ConstructSDandField() override;
  void SetTarget(TargetVariables target) {targetType = target;};

  void SetCADFilename(std::string name) {
    filename = name;
  };
  void SetCADFiletype(std::string type) {
    filetype = type;
  };
  
private:
  DetectorMessenger* fDetMess;
  G4ThreeVector offset;
  std::string filename;
  std::string filetype;

  TargetVariables targetType = TargetVariables::fWater;
      
  G4VSolid *cad_solid;
  G4LogicalVolume * cad_logical;
  G4VPhysicalVolume *cad_physical;
};

#endif
