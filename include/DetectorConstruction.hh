/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class (Mandatory)

#ifndef DETECTOR_CONSTRUCTION_HH
#define DETECTOR_CONSTRUCTION_HH

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
  
  void SetCADFilename(std::string name) {
        filename = name;
  };
  void SetCADFiletype(std::string type) {
        filetype = type;
  };
  
private:
  G4ThreeVector offset;
  std::string filename;
  std::string filetype;
      
  G4VSolid *cad_solid;
  G4LogicalVolume * cad_logical;
  G4VPhysicalVolume *cad_physical;
};

#endif
