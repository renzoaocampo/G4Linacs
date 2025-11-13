#ifndef PHYSICSLIST_HH
#define PHYSICSLIST_HH
#include "G4VModularPhysicsList.hh"
#include "G4EmStandardPhysics.hh" 
#include "G4ParticleTypes.hh"
#include "G4ProcessManager.hh"
#include "G4OpticalPhoton.hh"
#include "G4Decay.hh" 
#include "G4SystemOfUnits.hh"
 
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"

  
#include "G4VModularPhysicsList.hh"
#include "globals.hh"


class MyPhysicsList : public G4VModularPhysicsList {
public:
  MyPhysicsList();
  
  virtual ~MyPhysicsList();

};
 
#endif
