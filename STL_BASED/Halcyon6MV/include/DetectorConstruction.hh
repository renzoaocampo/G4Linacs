#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH
#include "G4SDManager.hh"

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"

#include "SensitiveDetector.hh"

class MyDetectorConstruction : public G4VUserDetectorConstruction
{ 
public:
    MyDetectorConstruction();
    ~MyDetectorConstruction();

    G4LogicalVolume *GetScoringVolume() const { return fScoringVolume; }

    virtual G4VPhysicalVolume *Construct();
   

private:
    G4Box *solidWorld,   *solidDetector   ;
    G4LogicalVolume *logicWorld,  *logicDetector  ,*logicmosfet,*logicVoxel ;
    G4VPhysicalVolume *physWorld, *physDetector  ;
 

    void DefineMaterials();
    virtual void ConstructSDandField();

    G4GenericMessenger *fMessenger;

    G4LogicalVolume *fScoringVolume;
 
};

#endif