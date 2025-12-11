// ====================================================================
// 1. MACRO CRÍTICA
// ====================================================================
#define USE_CADMESH_TETGEN 1

#include "DetectorConstruction.hh"
#include "CADMesh.hh"

#include "G4Region.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4AssemblyVolume.hh"
#include "G4VisAttributes.hh"
// NUEVO HEADER NECESARIO:
#include "G4LogicalVolumeStore.hh" 

#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;

// ====================================================================
// CONSTRUCTOR Y DESTRUCTOR
// ====================================================================
MyDetectorConstruction::MyDetectorConstruction()
 : G4VUserDetectorConstruction(), 
   fMessenger(nullptr),
   fScoringVolume(nullptr)
{}

MyDetectorConstruction::~MyDetectorConstruction()
{}

// ====================================================================
// MÉTODO CONSTRUCT
// ====================================================================
G4VPhysicalVolume* MyDetectorConstruction::Construct()
{
    // 1. Definir Materiales
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* worldMat = nist->FindOrBuildMaterial("G4_AIR");
    G4Material* mlcMat = nist->FindOrBuildMaterial("G4_Pb"); 
    





    // Pb/Sb para la placa
    G4double fractionmass;
    G4Element* elPb = nist->FindOrBuildElement("Pb");
    G4Element* elSb = nist->FindOrBuildElement("Sb");
    G4Material* matPbSb = new G4Material("PbSb_Alloy", 11.0*g/cm3, 2);
    matPbSb->AddElement(elPb, fractionmass=0.96);
    matPbSb->AddElement(elSb, fractionmass=0.04);


    G4Material* matW      = nist->FindOrBuildMaterial("G4_W");
    

    
    G4Material* matFe     = nist->FindOrBuildMaterial("G4_Fe");
    G4Material* matCarbon = nist->FindOrBuildMaterial("G4_C");
    G4Material* matCu     = nist->FindOrBuildMaterial("G4_Cu"); 

    // 2. Definir World
    G4Box* solidWorld = new G4Box("World", 2.*m, 2.*m, 2.*m);
    G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, worldMat, "World");
    logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible()); // World invisible
    
    G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, true);

    // 3. Estilo Visual (Sólido)
    G4VisAttributes* mlcVisAtt = new G4VisAttributes(G4Colour(0.6, 0.6, 0.7, 1.0)); 
    mlcVisAtt->SetVisibility(true);
    mlcVisAtt->SetForceSolid(true); // Relleno sólido




    G4double pb_radius = 30.15*mm;
    G4double pb_height = 1.0*mm;

    G4Tubs* solidPlatePbSb =
        new G4Tubs("PlatePbSb",0,pb_radius,pb_height/2,0,360*deg);

    G4LogicalVolume* logicPlatePbSb =
        new G4LogicalVolume(solidPlatePbSb, matPbSb,"PlatePbSb");

    new G4PVPlacement(
        0,
        G4ThreeVector(0,0,104.648*mm),
        logicPlatePbSb,"PlatePbSb",
        logicWorld,false,0,true);

    // ==========================================================
    // TARGET MCNP — superficies 13, 14, 15, 16
    // ==========================================================

    // --- Superficie 13 ---
    /*
        13 RCC 0 0 0   0 0 0.0635   0.301
    */
    {
        G4double r = 0.301*cm;
        G4double h = 0.0635*cm;

        G4Tubs* solid13 =
            new G4Tubs("TARGET_13",0,r,h/2,0,360*deg);

        G4LogicalVolume* logic13 =
            new G4LogicalVolume(solid13, matW, "TARGET_13");

        new G4PVPlacement(
            0, G4ThreeVector(0,0,h/2),
            logic13,"TARGET_13",logicWorld,false,0,true);
            
G4Region* TargetRegion = new G4Region("TargetRegion");
TargetRegion->AddRootLogicalVolume(logic13);
    }

    // --- Superficie 14 ---
    /*
        14 RCC 0 0 0.0635   0 0 1.016   0.301
    */
    {
        G4double r = 0.301*cm;
        G4double h = 1.016*cm;

        G4Tubs* solid14 =
            new G4Tubs("TARGET_14",0,r,h/2,0,360*deg);

        G4LogicalVolume* logic14 =
            new G4LogicalVolume(solid14, matCu, "TARGET_14");

        new G4PVPlacement(
            0, G4ThreeVector(0,0,0.0635*cm + h/2),
            logic14,"TARGET_14",logicWorld,false,0,true);
    }

 

    // ==========================================================
    // VISUALIZACIÓN
    // ========================================================== 
    // logicPrimColl->SetVisAttributes(new G4VisAttributes(G4Color(1,0,0)));
    logicPlatePbSb->SetVisAttributes(new G4VisAttributes(G4Color(0.3,0.3,0.3)));



























    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!! MLCs DESDE PLY !!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // 4. Bucle de Carga
    fs::path exec_cwd = fs::current_path();
    fs::path folder = exec_cwd.parent_path().parent_path() / "ply" / "MLCs";
    std::string folderPath = folder.string();

    if (!fs::exists(folder)) {
        G4cerr << "WARNING: carpeta PLY no encontrada: " << folderPath << G4endl;
    } else {
        G4cout << "Cargando PLY desde: " << folderPath << G4endl;
    }

    G4cout << "--- INICIANDO CARGA MASIVA DE PLY ---" << G4endl;

    for (const auto& entry : fs::directory_iterator(folderPath)) {
        if (entry.is_regular_file() && entry.path().extension() == ".ply") {
            
            std::string filePath = entry.path().string();
            G4cout << "Cargando: " << filePath << G4endl;

            auto mesh = CADMesh::TetrahedralMesh::FromPLY(filePath);
            
            mesh->SetScale(10.0); 
            mesh->SetOffset(G4ThreeVector(0, 0, 0));
            mesh->SetMaterial(mlcMat);

            G4AssemblyVolume* assembly = mesh->GetAssembly();

            // Colocación
            G4ThreeVector position = G4ThreeVector(0, 0, 0);
            G4RotationMatrix* rotation = new G4RotationMatrix();
            assembly->MakeImprint(logicWorld, position, rotation);
            
            delete rotation; 
        }
    }
    
    // -------------------------------------------------------------
    // 5. APLICAR VISUALIZACIÓN (SOLUCIÓN STORE)
    // -------------------------------------------------------------
    // Iteramos sobre TODOS los volúmenes lógicos creados en la simulación.
    // Si el volumen tiene el material de las MLCs, le aplicamos el estilo.
    // Esto evita el error de "GetTriplets" porque no tocamos el assembly.
    
    G4LogicalVolumeStore* store = G4LogicalVolumeStore::GetInstance();
    
    for (auto volume : *store) {
        // Verificamos si es un tetraedro de nuestras MLCs
        // (Podemos chequear por Material o por Nombre)
        if (volume->GetMaterial() == mlcMat) {
            volume->SetVisAttributes(mlcVisAtt);
        }
    }
    // -------------------------------------------------------------

    G4cout << "--- CARGA COMPLETA ---" << G4endl;

    

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    return physWorld;
}