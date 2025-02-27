// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
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
#include "SensitiveSD.hh"

#include <sstream>

DetectorConstruction::DetectorConstruction() :
    G4VUserDetectorConstruction()
{

}

DetectorConstruction::~DetectorConstruction()
{
    ;
}


G4VPhysicalVolume* DetectorConstruction::Construct()
{

    G4bool checkOverlaps = true;

        // Define Material

    G4NistManager* man = G4NistManager::Instance();
     
      G4Element* Li = man->FindOrBuildElement("Li");
      G4Element* N  = man->FindOrBuildElement("N");
      G4Element* La = man->FindOrBuildElement("La");
      G4Element* Br = man->FindOrBuildElement("Br");
      G4Element* Ce = man->FindOrBuildElement("Ce");
	
	  G4Material *VetoMat = man->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

      // Sand sediments ----------------------------
      G4Element* H  = man->FindOrBuildElement("H");  G4Element* B  = man->FindOrBuildElement("B");
      G4Element* C  = man->FindOrBuildElement("C");  G4Element* O  = man->FindOrBuildElement("O");
      G4Element* F  = man->FindOrBuildElement("F");  G4Element* Na = man->FindOrBuildElement("Na");
      G4Element* Mg = man->FindOrBuildElement("Mg"); G4Element* Al = man->FindOrBuildElement("Al");
      G4Element* Si = man->FindOrBuildElement("Si"); G4Element* S  = man->FindOrBuildElement("S");     // eror in code
      G4Element* Cl = man->FindOrBuildElement("Cl"); G4Element* K  = man->FindOrBuildElement("K");
      G4Element* Ca = man->FindOrBuildElement("Ca"); G4Element* V  = man->FindOrBuildElement("V");
      G4Element* Cr = man->FindOrBuildElement("Cr"); G4Element* Fe = man->FindOrBuildElement("Fe");
      G4Element* Co = man->FindOrBuildElement("Co"); G4Element* Ni = man->FindOrBuildElement("Ni");
      G4Element* Cu = man->FindOrBuildElement("Cu"); G4Element* Zn = man->FindOrBuildElement("Zn");
      G4Element* As = man->FindOrBuildElement("As"); G4Element* Sr = man->FindOrBuildElement("Sr");
      G4Element* Cd = man->FindOrBuildElement("Cd"); G4Element* Ba = man->FindOrBuildElement("Ba");
      G4Element* Hg = man->FindOrBuildElement("Hg"); G4Element* Pb = man->FindOrBuildElement("Pb");
      
      //
      // Sand Sediment - 27 elemental composition
      G4Material* SandSediment = new G4Material("SandSediment", 1.535*g/cm3, 27);
      SandSediment->AddElement(H, 3.30144*perCent);
      SandSediment->AddElement(B, 0.00003*perCent); SandSediment->AddElement(C, 0.35344*perCent);
      
      SandSediment->AddElement(O, 61.9999*perCent);
      SandSediment->AddElement(F, 0.00001*perCent); SandSediment->AddElement(Na,1.34864*perCent);
      
      SandSediment->AddElement(Mg,0.00734*perCent);
      SandSediment->AddElement(Al,2.59482*perCent); SandSediment->AddElement(Si,27.5365*perCent);

      SandSediment->AddElement(S, 0.00517*perCent);
      SandSediment->AddElement(Cl,0.11055*perCent);
      SandSediment->AddElement(K, 0.87891*perCent);
      
      SandSediment->AddElement(Ca,0.84806*perCent);
      SandSediment->AddElement(V, 0.00100*perCent);
      SandSediment->AddElement(Cr,0.00100*perCent);
      
      SandSediment->AddElement(Fe,1.00237*perCent);
      SandSediment->AddElement(Co,0.00050*perCent);
      SandSediment->AddElement(Ni,0.00100*perCent);
      
      SandSediment->AddElement(Cu,0.00100*perCent);
      SandSediment->AddElement(Zn,0.00200*perCent);
      SandSediment->AddElement(As,0.00080*perCent);
      
      SandSediment->AddElement(Br,0.00038*perCent);
      SandSediment->AddElement(Sr,0.00104*perCent);
      SandSediment->AddElement(Cd,0.00010*perCent);
      
      SandSediment->AddElement(Ba, 0.00250*perCent);
      SandSediment->AddElement(Hg, 0.000002*perCent);
      SandSediment->AddElement(Pb,0.00150*perCent);
      

      
      G4Material *vacuum = man->FindOrBuildMaterial("G4_Galactic");

    //////////////////////////// Air
    G4int ncomponents;
    G4double density,temperature,pressure,fractionmass;

                    density     = 0.001225*g/cm3;
                    pressure    = 98658.96*pascal;
                    temperature = 273*kelvin;

		    G4Material* Air = new G4Material("Air", density, ncomponents=2, kStateGas, temperature, pressure);
                    Air->AddElement(N, fractionmass=79.*perCent);
                    Air->AddElement(O, fractionmass=21.*perCent);
                    G4NistManager *nist_man = G4NistManager::Instance();
                    G4Material *AIR_mat = nist_man->FindOrBuildMaterial("Air");
    /////////////////////////////////////////////////
      G4Material* H2O =
      new G4Material("Water", 1.000*g/cm3, 2);
      H2O->AddElement(H, 2);
      H2O->AddElement(O, 1);
      H2O->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);

      G4Material* H3BO3 =
      new G4Material("Boric_Acid", 1.44*g/cm3, 3);
      H3BO3->AddElement(H, 3);
      H3BO3->AddElement(B, 1);
      H3BO3->AddElement(O, 3);
      //GetIonisation()->SetMeanExcitationEnergy(78.0*eV);

      G4Material* NaCl =
      new G4Material("Salt", 2.165*g/cm3, 2);
      NaCl->AddElement(Na, 1);
      NaCl->AddElement(Cl, 1);

      G4Material* LiF =
      new G4Material("Lithium_Fluoride", 2.64*g/cm3, 2);
      LiF->AddElement(Li, 1);
      LiF->AddElement(F, 1);

                   
      //fDefaultMaterial = Galactic;

        //SeaWater
        G4Material* SeaWater =
        new G4Material("SeaWater", 1.027*g/cm3, 2);
        SeaWater->AddMaterial(H2O,98.8*perCent);
        SeaWater->AddMaterial(NaCl,1.2*perCent);

        //  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

        // LaBr3
        G4Material* LaBr3 = new G4Material("LaBr3", 5.08*g/cm3, 2);
        LaBr3->AddElement(La, 1);
        LaBr3->AddElement(Br, 3);

        // LaBr3_Ce
        G4Material* LaBr3_Ce = new G4Material("LaBr3_Ce", 5.29*g/cm3, 2);
        LaBr3_Ce->AddMaterial(LaBr3, 95*perCent);
        LaBr3_Ce->AddElement(Ce, 5*perCent);

        //  Mustard Gas
        G4Material* C4H8Cl2S = new G4Material("Mustard_Gas", 1.27*g/cm3, 4);
        C4H8Cl2S->AddElement(C,4);
        C4H8Cl2S->AddElement(H,8);
        C4H8Cl2S->AddElement(Cl,2);
        C4H8Cl2S->AddElement(S,1);


	// SiO2
        G4Material* SiO2 = new G4Material("Silicon_Dioxide", 2.196 *g/cm3, 2);
        SiO2->AddElement(Si,1);
        SiO2->AddElement(O,2);

	// Fe
        G4Material* Iron = new G4Material("Iron", 7.850 *g/cm3, 1);
        Iron->AddElement(Fe,1);
        // Lead
        G4Material* Lead = new G4Material("Lead", 11.34 *g/cm3, 1);
        Lead->AddElement(Pb,1);
        // Paraffin
        G4Material* Paraffin = man->FindOrBuildMaterial("G4_PARAFFIN");
        //Aluminum
        G4Material* Aluminum = new G4Material("Al", 2.7 *g/cm3, 1);
        Aluminum->AddElement(Al,1);
        //Paraffin_Boric_Acid
        G4Material* Paraffin_Boric_Acid =
        new G4Material("Paraffin_Boric_Acid", 0.972*g/cm3, 2);
        Paraffin_Boric_Acid->AddMaterial(Paraffin,90*perCent);
        Paraffin_Boric_Acid->AddMaterial(H3BO3,10*perCent);

        //Cadmium
        G4Material* Cadmium = new G4Material("Cadmium", 8.65 *g/cm3, 1);
        Cadmium->AddElement(Cd,1);
        //Carbon(graphite)
        G4Material* Carbon = new G4Material("Carbon", 2.3 *g/cm3, 1);
        Carbon->AddElement(C,1);

        //Low_Carbon_Stainless_Steel
        G4Material* Low_Carbon_Steel = new G4Material("Low_Carbon_Steel", 7.8499*g/cm3, 2);
        Low_Carbon_Steel->AddMaterial(Iron,99.75*perCent);
        Low_Carbon_Steel->AddMaterial(Carbon,0.25*perCent);



	G4Material* LCSt = man->FindOrBuildMaterial("Low_Carbon_Steel");
        G4Material* MGas = man->FindOrBuildMaterial("Mustard_Gas");
        G4Material* SandB= man->FindOrBuildMaterial("Silicon_Dioxide");
        G4Material* DetM = man->FindOrBuildMaterial("LaBr3");
	
	
    // rotation of objects
	//**************************
	G4RotationMatrix* rotation = new G4RotationMatrix();
		
	rotation->rotateY(90*deg);
	G4RotationMatrix* rotation1 = new G4RotationMatrix();
	
	rotation1->rotateX(90*deg);

	G4RotationMatrix *detRot = new G4RotationMatrix();
	  detRot->rotateX(45*deg);
	 
	
	G4RotationMatrix* rotation2 = new G4RotationMatrix();
	rotation2->rotateZ(0*deg);
	rotation2->rotateX(-90*deg);
	rotation2->rotateY(180*deg);
	
	
	G4RotationMatrix* rotation_dummy = new G4RotationMatrix();
	rotation_dummy->rotateZ(0*deg);
	rotation_dummy->rotateY(180*deg);
	
	
	G4RotationMatrix* rotation3 = new G4RotationMatrix();
		//   rotation->rotateZ(3*deg);
	rotation3->rotateX(55*deg);
	rotation3->rotateY(90*deg);
	
	G4RotationMatrix* rotation4 = new G4RotationMatrix();
	rotation4->rotateX(63*deg);
	rotation4->rotateY(90*deg);
	
	G4RotationMatrix* rotation5 = new G4RotationMatrix();
	rotation5->rotateY(180.*deg);
	rotation5->rotateX(-180.*deg);
	
		//rotation for gamma guide
	G4RotationMatrix* rotation6 = new G4RotationMatrix();
		//rotation6->rotateY(130.*deg);
	rotation6->rotateX(120.*deg);
		//  rotation6->rotateZ(-40.*deg);
	
    //-------------------------

  static const double Pi = 3.14159265358979323846;

  // World
    //
      G4Box* solidWorld =
      new G4Box("World",                          //name
                 4.2*m/2,4.2*m/2,4.2*m/2);     //size

      G4LogicalVolume* logicWorld =
      new G4LogicalVolume(solidWorld,             //solid
			       SeaWater,               //material
			  // AIR_mat,
                          "World");               //name

      G4VPhysicalVolume* physWorld=
      new G4PVPlacement(0,                        //no rotation
                        G4ThreeVector(),          //at (0,0,0)
                        logicWorld,               //logical volume
                        "World",                  //name
                         0,                       //mother volume
                         false,                   //no boolean operation
                         0);                      //copy number

      G4VisAttributes* worldVisAtt= new G4VisAttributes(G4Colour(1.,1.,0.5));
      worldVisAtt->SetForceWireframe(true);
     // worldVisAtt->SetVisibility(false);
      logicWorld->SetVisAttributes(worldVisAtt);



    // SUBMARINE VOLUME ++++++++++++++++++++++++++++++++++++

      G4Box* solidSubmarineVolumeOuter = new G4Box("SubmarineVolumeOuter",1.5*m,20*cm,1.5*m);
   
      G4LogicalVolume* logicSubmarineVolume = new G4LogicalVolume(solidSubmarineVolumeOuter,
					//		  AIR_mat,
		 	                                LCSt,
							"SubmarineVolumeOuter");

      new G4PVPlacement(0,
                        G4ThreeVector(0.,0.,0),
                        logicSubmarineVolume,
                        "SubmarineVolumeOuter",
                        logicWorld,
                        false,
                        1,
                        checkOverlaps);

      //Added by sks - visualization for the submarine
      G4VisAttributes* subMarVisAttOuter= new G4VisAttributes(G4Colour(0.5,0.5,1.));
      subMarVisAttOuter->SetForceWireframe(true);
     // subMarVisAttOuter->SetForceSolid(true);
      logicSubmarineVolume->SetVisAttributes(subMarVisAttOuter);

	// Put air volue inside the submarine - need
	
	G4Box* solidSubVolume = new G4Box("SubVolume",1.5*m-3*mm,20*cm-3*mm,1.5*m-3*mm);
	G4LogicalVolume* logicSubVolume = new G4LogicalVolume(solidSubVolume, vacuum,"SubVolume");
	
    new G4PVPlacement(0,
					  G4ThreeVector(0.,0*mm,0.),
					  logicSubVolume,
					  "SubVolume",
					  logicSubmarineVolume,
					  false,
					  2,
					  checkOverlaps);
	
	G4VisAttributes* subbMarVisAtt2 = new G4VisAttributes(G4Color(1.0,0.,0.5));
	subbMarVisAtt2->SetForceWireframe(true);
	logicSubVolume->SetVisAttributes(subbMarVisAtt2);

	/*	
    // Plastic Veto detector ->
    //
	G4Box* solidVeto = new G4Box("Veto",10*cm,20*mm,10*cm);
	G4LogicalVolume* logicVeto = new G4LogicalVolume(solidVeto,VetoMat,"Veto");
	
	new G4PVPlacement(0,
					  G4ThreeVector(0.,-10*cm,0.),
					  logicVeto,
					  "Veto",
					  logicSubmarineVolume,
					  false,
					  13,
					  checkOverlaps);
	
	G4VisAttributes* VetoAtt = new G4VisAttributes(G4Color(1.0,1.,0.5));
	VetoAtt->SetForceWireframe(true);
	logicVeto->SetVisAttributes(VetoAtt);*/

    //Neutron SourceGuide++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Cone is replaced by the cuboid ( 6 faces where each face is diagonal
	G4VSolid* solidSourceGuide  = new G4Tubs("SourceGuide",9.5*cm,10.*cm,5.0*cm,0.,2*M_PI*rad);
  
//	G4Box* solidSourceGuide = new G4Box("SourceGuide",10*cm,10*cm,5*cm);

      G4LogicalVolume* logicSourceGuide = new G4LogicalVolume(solidSourceGuide,
                                              //              AIR_mat,
                                                              LCSt,
                                                              "SourceGuide");
      new G4PVPlacement(rotation1,
                        G4ThreeVector(0.,-25*cm,0),
                        logicSourceGuide,
                        "SourceGuide",
                        logicWorld,
                        false,
                        3,
                        checkOverlaps);
      //Added by sks
      G4VisAttributes* sourceGVisAtt= new G4VisAttributes(G4Colour(1.0,0.0,1.0));
      //sourceGVisAtt->SetForceSolid(true);
      sourceGVisAtt->SetForceSolid(true);
      logicSourceGuide->SetVisAttributes(sourceGVisAtt);
	

	// Target (Mustard Gas)+++++++++++++++++++++++++++++++++++++++++++++

      // Steel bos for mustard gas of 3mm thickness

      G4Box* solidTGVolumeCover = new G4Box("TargetVolumeCover",97*cm, 25*cm, 25*cm);
      G4LogicalVolume* logicTGVolumeCover = new G4LogicalVolume(solidTGVolumeCover, LCSt, "TargetVolumeCover");
      new G4PVPlacement(rotation,
			G4ThreeVector(0.,-55.*cm,0.),
			logicTGVolumeCover,
			"TargetVolumeCover",
			/*logicSANDVolume*/
			logicWorld,
			false,
			4,
			checkOverlaps);


 
	G4Box* solidTGVolume = new G4Box("TargetVolume",97*cm-3.*mm,25*cm-3.*mm,25*cm-3.*mm);
	G4LogicalVolume* logicTGVolume = new G4LogicalVolume(solidTGVolume,
						//      		     MGas,
							              SeaWater,
				   			     "TargetVolume");
	new G4PVPlacement(0,
			  G4ThreeVector(),
			  logicTGVolume,
			  "TargetVolume",
			   logicTGVolumeCover,
			   false,
			   5,
			   checkOverlaps);
	
	
	
	G4VisAttributes* mustardVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,0));
	mustardVisAtt->SetForceSolid(true);
	logicTGVolume->SetVisAttributes(mustardVisAtt);

	 //Gamma ParticleGuide //trapezoidal
	
	// Trap for gammas -
	// Made of low carbon sttel - >

	// parameter -  name, a, b, c, d, e, f, g, h, i, j, k
	
	// a - half length z - 50 mm
	// b - theta - 69 deg
	// c - phi - 90 deg

	//d,e,f are for the lower base
	// h, i, j are for the upper base
	
	
		
	G4Trap* solidParticleGuideOuterCov = new G4Trap("ParticleGuideOuterCov",5.*cm,69*deg,90*deg,10.*cm+3*mm,
						     10.*cm+3*mm,10.*cm+3*mm,0*deg,10.*cm+3*mm,10.*cm+3*mm,10.*cm+3*mm,0*deg);
	
	
	G4Trap* solidParticleGuideIns = new G4Trap("ParticleGuideIns",5.*cm+0.01*mm,    // 0.01 is used to better subtraction
											                      69*deg,
											                      90*deg,
																  10.*cm,10.*cm,10.*cm,    // Lower base
											                      0*deg,
											                      10.*cm,10.*cm,10.*cm,    // Upper base
											                      0*deg);
	
	G4VSolid* solidParticleGuideOuter = new G4SubtractionSolid("ParticleGuideOuterCov-ParticleGuideIns", solidParticleGuideOuterCov, solidParticleGuideIns);
	

	G4LogicalVolume* logicParticleGuideOuter = new G4LogicalVolume(solidParticleGuideOuter,
								  LCSt,
								  "ParticleGuideOuter");
	new G4PVPlacement(rotation2,
			  G4ThreeVector(0.*cm,-25*cm,33.5*cm),
			  logicParticleGuideOuter,
			  "ParticleGuideOuter",
			  logicWorld,
			  false,
			  6,
			  checkOverlaps);

	G4VisAttributes* particleVisAttOuter= new G4VisAttributes(G4Colour(0.0,1.0,0.5));
	particleVisAttOuter->SetForceSolid(true);
	logicParticleGuideOuter->SetVisAttributes(particleVisAttOuter);

	
// Detector:
//	G4Box* solidAbsor = new G4Box("DetectorLaBr",2.54*cm,2.54*cm,2.54*cm);
	G4VSolid* solidAbsor  = new G4Tubs("DetectorLaBr",0,2.54*cm,2.54*cm,0.,2*M_PI*rad);
	G4LogicalVolume* logicAbsor =new G4LogicalVolume(solidAbsor,           
						       	 LaBr3_Ce,             
			        			 "DetectorLaBr");
		//(G4ThreeVector((first coordinate work as z,second as y and third as z)
	G4ThreeVector position = G4ThreeVector(0.*cm,-0.19855*m+4.8*cm,+51.*cm); // 4.8 is height and 50 cm is distance from source
	
	
	new G4PVPlacement(detRot,          //no rotation
			  position,          //position
			  logicAbsor,        //logical volume
			  "DetectorLaBr",    //name
			   logicSubVolume,       //mother
			   false,             //no boulean operat
			   8,
			   checkOverlaps);   //copy number
	
	G4VisAttributes* solidVisAtt= new G4VisAttributes(G4Colour(1.0,0.0,0.0));
	solidVisAtt->SetForceSolid(true);
	logicAbsor->SetVisAttributes(solidVisAtt);
	
	
    // SAND ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	G4Box* solidSANDVolume = new G4Box("SAND",2*m,50*cm,2*m);
	G4LogicalVolume* logicSANDVolume = new G4LogicalVolume(solidSANDVolume,
							       SandSediment,
							      "SAND");
	new G4PVPlacement(0,
					  G4ThreeVector(0.,-130*cm,0.),
					  logicSANDVolume,
					  "SAND",
					  logicWorld,
					  false,
					  9,
					  checkOverlaps);
	
		
	   
	G4VisAttributes* sandVisAtt= new G4VisAttributes(G4Colour(0.5,0.5,0.5));
    sandVisAtt->SetForceSolid(true);
	//sandVisAtt->SetForceWireframe(true);
	logicSANDVolume->SetVisAttributes(sandVisAtt);
	
	/*
 //Source boarder for presen. 

      G4double pRmin=0.*cm;
      G4double pRmax=2.5*cm;
      G4double pSPhi=0*deg;
      G4double pDPhi=2*Pi;
      G4double pSTheta=0*deg;
      G4double pDTheta=Pi;

      G4Sphere* solidShape1 = new G4Sphere("SOURCE", pRmin, pRmax, pSPhi, pDPhi, pSTheta, pDTheta);
      G4LogicalVolume* logicShape1 = new G4LogicalVolume(solidShape1, AIR_mat, "SOURCE"); 
      new G4PVPlacement(0, G4ThreeVector(0,-15*cm,0), logicShape1, "SOURCE", logicSubVolume, false, 10);
      
      G4VisAttributes * sourceVisAtt = new G4VisAttributes(G4Colour(1.,1.,0.));
      sourceVisAtt->SetForceWireframe(true);
      logicShape1->SetVisAttributes(sourceVisAtt);*/
	
	 // make plastic scintillators as alpha particle tagging -
	
	 G4Box* AlphaDet = new G4Box("AlphaDetStrips", scinDim_x/2.0 ,scinDim_y/2.0 , scinDim_z/2.0 );
	 G4LogicalVolume * AlphaDetLog = new G4LogicalVolume(AlphaDet, VetoMat , "AlphaDetStrips");
	
	 G4VisAttributes* BoxVisAtt =  new G4VisAttributes(G4Colour(1.,0.4,.9));
	 BoxVisAtt->SetForceWireframe(true);
	 BoxVisAtt->SetForceSolid(true);
	 AlphaDetLog->SetVisAttributes(BoxVisAtt);
	
	
	
	for(int j=0;j<7;j++)
		{
		
		G4ThreeVector loc = G4ThreeVector(-2.1*cm+j*0.7*cm,-12.*cm,0.0);
		new G4PVPlacement(0,
						  loc,             //rotation,position
						  AlphaDetLog,            //its logical volume
						  "AlphaDetStrips",             //its name
						  logicSubVolume,      //its mother (logical) volume
						  true,                 //no boolean operation
						  j+10,                 //copy number
						  checkOverlaps);       // checking overlaps
		
		}
		
			
	
// Iron shielding

      G4Box* ironShield = new G4Box("Shield", 10.*cm/2, 5.*cm/2, 20.*cm/2);     
      G4LogicalVolume* logicIron = new G4LogicalVolume(ironShield, Iron, "Shield");           
      new G4PVPlacement(0, G4ThreeVector(0,-15.0*cm,24.*cm), logicIron, "Shield", logicSubVolume, false, 11); 

      G4VisAttributes* ironVisAtt= new G4VisAttributes(G4Colour(1.,1.,0.5));
      ironVisAtt->SetForceWireframe(true);
      ironVisAtt->SetForceSolid(true);
       logicIron->SetVisAttributes(ironVisAtt);

// Lead Cap
       G4Box* LeadCap = new G4Box("PbShield", 10.*cm/2, 5.*cm/2, 4.*cm/2);
       G4LogicalVolume* logicLead = new G4LogicalVolume(LeadCap, Lead, "PbShield");
       new G4PVPlacement(0, G4ThreeVector(0,-15.0*cm,36.*cm), logicLead, "Shield", logicSubVolume, false, 12);

       G4VisAttributes* LeadVisAtt= new G4VisAttributes(G4Colour(0.2,1.,0.5));
       LeadVisAtt->SetForceWireframe(true);
       LeadVisAtt->SetForceSolid(true);
       logicLead->SetVisAttributes(LeadVisAtt);

      
    return physWorld;
}


// Implement the following for sensitive detector
void DetectorConstruction::ConstructSDandField()
{
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    sdManager->SetVerboseLevel(2);

    SensitiveSD* myDetector = new SensitiveSD("Detector");   
	
    SetSensitiveDetector("DetectorLaBr", myDetector); 
   
    sdManager->AddNewDetector(myDetector);

}
