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
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4AssemblyVolume.hh"
#include "G4VisAttributes.hh"
// NUEVO HEADER NECESARIO:
#include "G4LogicalVolumeStore.hh" 

#include <filesystem>
#include <iostream>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <cctype>
#include <memory>

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
    G4Material* mlcMat = nist->FindOrBuildMaterial("G4_W");
    




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
    G4Box* solidWorld = new G4Box("World", 5.*m, 5.*m, 5.*m);
    G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, worldMat, "World");
    logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible()); // World invisible
    
    G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, true);

    // 3. Estilo Visual (Sólido)
    G4VisAttributes* mlcVisAtt = new G4VisAttributes(G4Colour(0.6, 0.6, 0.7, 1.0)); 
    mlcVisAtt->SetVisibility(true);
    mlcVisAtt->SetForceSolid(true); // Relleno sólido

    // Atributo universal: todos los volúmenes visibles y sólidos
    G4VisAttributes* solidVisAtt = new G4VisAttributes();
    solidVisAtt->SetVisibility(true);
    solidVisAtt->SetForceSolid(true);






    

    // ==========================================================
    // COLIMADOR PRIMARIO (tu forma actual, no la MCNP)
    // ==========================================================
    G4Cons* solidPrimColl =
        new G4Cons("PrimColl",
                   0.4*cm,       // rmin1
                   3.557*cm,     // rmax1
                   1.945*cm,     // rmin2
                   3.557*cm,     // rmax2
                   5.1*cm,       // hz
                   0, 360*deg);

    G4LogicalVolume* logicPrimColl =
        new G4LogicalVolume(solidPrimColl, matW, "PrimColl");

    new G4PVPlacement(
        0, G4ThreeVector(0,0,2.7*cm),
        logicPrimColl,"PrimColl",logicWorld,false,0,true);



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

   
//*========= Revestimiento =========
G4Tubs* revestimiento = new G4Tubs("Revestimiento",
    63.26 * cm, 68.34 * cm,
    59.625 * cm / 2, 0, 360 * deg);

G4LogicalVolume* logicRevestimiento = new G4LogicalVolume(
    revestimiento, matFe, "LogicRevestimiento");

new G4PVPlacement(0, G4ThreeVector(0, 0,   20.6325 * cm),
    logicRevestimiento, "Revestimiento", logicWorld, false, 0, true);



    // ==========================================================
    // COLIMADORES PRIMARIOS (Primary Collimators)
    // Material: m91 (W puro) = G4_W
    // ==========================================================

    // Declare logical volumes at function scope
 
    
    // ==========================================================
    // VISUALIZACIÓN
    // ========================================================== 
    // logicPrimColl->SetVisAttributes(new G4VisAttributes(G4Color(1,0,0)));
    logicPlatePbSb->SetVisAttributes(new G4VisAttributes(G4Color(0.3,0.3,0.3)));
 
   

























    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!! PLYs DESDE CARPETA ply !!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // 4. Bucle de Carga (recursivo desde <exec>/../../ply)
    fs::path exec_cwd = fs::current_path();
    fs::path ply_root = exec_cwd.parent_path().parent_path() / "ply";

    if (!fs::exists(ply_root)) {
        G4cerr << "WARNING: carpeta raíz ply no encontrada: " << ply_root.string() << G4endl;
    } else {
        G4cout << "Buscando PLY recursivamente en: " << ply_root.string() << G4endl;
    }

    // Crear materiales específicos y mapa carpeta->material
    G4Element* elH = nist->FindOrBuildElement("H");
    G4Element* elC = nist->FindOrBuildElement("C");
    G4Element* elO = nist->FindOrBuildElement("O");

    // Treatment bed (por ejemplo carpeta "camilla")
    G4Material* matTreatmentBed = new G4Material("TreatmentBed", 1.0*g/cm3, 3);
    matTreatmentBed->AddElement(elH, 0.057444);
    matTreatmentBed->AddElement(elC, 0.774589);
    matTreatmentBed->AddElement(elO, 0.167968);

    // Base de mesa (por ejemplo carpeta "mesa" o "base")
    G4Material* matBaseMesa = nist->FindOrBuildMaterial("G4_Fe");

    // Mapa carpeta -> material (usar nombres en minúsculas)
    std::unordered_map<std::string, G4Material*> folderMaterialMap;
    folderMaterialMap["camilla"] = matTreatmentBed;
    folderMaterialMap["mesa"] = matBaseMesa;
    folderMaterialMap["base"] = matBaseMesa;

    G4cout << "--- INICIANDO CARGA MASIVA DE PLY (recursiva por carpetas) ---" << G4endl;

    if (fs::exists(ply_root)) {
        for (const auto& entry : fs::recursive_directory_iterator(ply_root)) {
            if (!entry.is_regular_file()) continue;
            if (entry.path().extension() != ".ply") continue;

            std::string filePath = entry.path().string();
            G4cout << "Cargando: " << filePath << G4endl;

            // Obtener carpeta padre y normalizar a minúsculas
            std::string parentFolder = entry.path().parent_path().filename().string();
            std::transform(parentFolder.begin(), parentFolder.end(), parentFolder.begin(),
                           [](unsigned char c){ return std::tolower(c); });

            // Intentar cargar la malla de forma segura (FromPLY devuelve std::shared_ptr)
            std::shared_ptr<CADMesh::TetrahedralMesh> meshPtr;
            try {
                meshPtr = CADMesh::TetrahedralMesh::FromPLY(filePath);
            } catch (const std::exception& e) {
                G4cerr << "ERROR: excepción al cargar PLY: " << e.what() << " -> " << filePath << G4endl;
                continue;
            } catch (...) {
                G4cerr << "ERROR: excepción desconocida al cargar PLY -> " << filePath << G4endl;
                continue;
            }

            if (!meshPtr) {
                G4cerr << "ERROR: mesh nula para " << filePath << G4endl;
                continue;
            }

            auto mesh = meshPtr.get();

            mesh->SetScale(10.0);
            mesh->SetOffset(G4ThreeVector(0, 0, 0));

            // Seleccionar material según carpeta
            G4Material* chosenMat = mlcMat; // por defecto
            auto it = folderMaterialMap.find(parentFolder);
            if (it != folderMaterialMap.end() && it->second) {
                chosenMat = it->second;
                G4cout << "  -> Carpeta: " << parentFolder << " -> material asignado." << G4endl;
            } else {
                G4cout << "  -> Carpeta: " << parentFolder << " -> material por defecto (MLC Pb)." << G4endl;
            }

            mesh->SetMaterial(chosenMat);

            G4AssemblyVolume* assembly = mesh->GetAssembly();
            if (!assembly) {
                G4cerr << "ERROR: assembly nulo para " << filePath << G4endl;
                continue;
            }

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
    // Aplicamos estilos visuales específicos.
    
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