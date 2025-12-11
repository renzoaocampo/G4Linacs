// ====================================================================
// 1. MACRO CRÍTICA: Debe ser la primera línea absoluta del archivo.
// Esto activa la lectura de 'TetrahedralMesh' dentro de CADMesh.hh
// ====================================================================
#define USE_CADMESH_TETGEN 1

#include "DetectorConstruction.hh"
#include "CADMesh.hh" // Debe ir DESPUÉS del define

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4AssemblyVolume.hh"
#include "G4VisAttributes.hh"

// ====================================================================
// 2. IMPLEMENTACIÓN DE CONSTRUCTOR Y DESTRUCTOR (Faltaban en tu código)
// ====================================================================

// Constructor vacío (necesario porque lo declaraste en el .hh)
MyDetectorConstruction::MyDetectorConstruction()
 : G4VUserDetectorConstruction(), 
   fMessenger(nullptr),
   fScoringVolume(nullptr)
{
}

// Destructor vacío
MyDetectorConstruction::~MyDetectorConstruction()
{
}

// ====================================================================
// 3. MÉTODO CONSTRUCT
// ====================================================================
G4VPhysicalVolume* MyDetectorConstruction::Construct()
{
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* worldMat = nist->FindOrBuildMaterial("G4_AIR");
    G4Material* boneMat  = nist->FindOrBuildMaterial("G4_Pb");
    G4Material* waterMat = nist->FindOrBuildMaterial("G4_WATER");

    // --- World ---
    G4Box* solidWorld = new G4Box("World", 2.*m, 2.*m, 2.*m);
    G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, worldMat, "World");
    G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, true);

    // --- Carga Pelvis (Malla Tetraédrica) ---
    //
    // CAMBIO REALIZADO: Usar FromPLY en lugar de FromSTL
    auto meshPelvis = CADMesh::TetrahedralMesh::FromPLY(
        "C:\\Users\\renzo\\Downloads\\MLC_HALCYON_STL\\ply\\304.ply");
    meshPelvis->SetScale(1*cm); 
    meshPelvis->SetOffset(G4ThreeVector(0, 0, 0));
    meshPelvis->SetMaterial(boneMat);

    G4AssemblyVolume* assemblyPelvis = meshPelvis->GetAssembly();

    // Variables explícitas para evitar errores de referencia en Windows
    G4ThreeVector posPelvis = G4ThreeVector(0, 0, 0);
    G4RotationMatrix* rotPelvis = new G4RotationMatrix();
    assemblyPelvis->MakeImprint(logicWorld, posPelvis, rotPelvis);

   
     

    return physWorld;
}