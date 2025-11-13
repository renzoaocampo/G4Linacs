#include "PhysicsList.hh" 

#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"

#include "G4UnitsTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

// Bosons 
#include "G4Geantino.hh"
#include "G4Gamma.hh"
#include "G4OpticalPhoton.hh"

// Leptons 
#include "G4Electron.hh"
#include "G4Positron.hh"

// Hadrons
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"

#include "G4SystemOfUnits.hh"
#include "G4EmParameters.hh"  // Faltaba incluir esto para usar emParameters

MyPhysicsList::MyPhysicsList() 
: G4VModularPhysicsList() 
{
  // Default cut value
  defaultCutValue = 0.05* mm;

  // Configurar parámetros EM
  auto emParameters = G4EmParameters::Instance();


  // Registrar física electromagnética Livermore
  RegisterPhysics(new G4EmStandardPhysics_option4());
  // emParameters->SetDefaults(); // (Opcional) para asegurarse de partir de parámetros estándar
  // emParameters->SetAuger(true); // Activar electrones Auger
  // emParameters->SetFluo(true);  // Activar fluorescencia
  // emParameters->SetDeexcitationIgnoreCut(true); // Permitir deexcitación por debajo del cut
  // emParameters->SetPixe(true); // <<<< ¡Esto activa el PIXE!
  // (Opcional) Configurar el rango de energía de corte
  // G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100*eV, 1*GeV);
}

MyPhysicsList::~MyPhysicsList() {}
