// ====================================================================
// 1. MACRO CRÍTICA
// ====================================================================
#define USE_CADMESH_TETGEN 1

#include "DetectorConstruction.hh"
#include "CADMesh.hh"

#include "SensitiveDetector.hh"
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
#include <cstdlib>

namespace fs = std::filesystem;

// ====================================================================
// CONSTRUCTOR Y DESTRUCTOR
// ====================================================================
MyDetectorConstruction::MyDetectorConstruction() {} 
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

    // Valor por defecto para offset Z de la mesa (cm). Se puede sobreescribir
    // con la variable de entorno G4LINACS_MESA_Z_OFFSET_CM (valor en cm).
    double mesaOffsetCm = 20.0;
    if (const char* env = std::getenv("G4LINACS_MESA_Z_OFFSET_CM")) {
        try {
            mesaOffsetCm = std::stod(env);
            G4cout << "G4LINACS_MESA_Z_OFFSET_CM detectado: " << mesaOffsetCm << " cm" << G4endl;
        } catch (...) {
            G4cerr << "WARNING: G4LINACS_MESA_Z_OFFSET_CM no numérico, usando " << mesaOffsetCm << " cm" << G4endl;
        }
    }

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
            G4Material* chosenMat = mlcMat; // por defecto: Tungsteno (MLCs)
            auto it = folderMaterialMap.find(parentFolder);
            if (it != folderMaterialMap.end() && it->second) {
                chosenMat = it->second;
                G4cout << "  -> Carpeta: " << parentFolder << " -> material asignado." << G4endl;
            } else {
                G4cout << "  -> Carpeta: " << parentFolder << " -> material por defecto (Tungsteno - MLCs)." << G4endl;
            }

            mesh->SetMaterial(chosenMat);

            G4AssemblyVolume* assembly = mesh->GetAssembly();
            if (!assembly) {
                G4cerr << "ERROR: assembly nulo para " << filePath << G4endl;
                continue;
            }

            // Colocación
            G4ThreeVector position = G4ThreeVector(0, 0, 0);
            // Si el PLY pertenece a la carpeta camilla, desplazar en Z según variable (por defecto +20 cm)
            if (parentFolder == "camilla" || parentFolder.find("camilla") != std::string::npos) {
                position += G4ThreeVector(0, 0, mesaOffsetCm*cm);
                G4cout << "  -> Aplicando offset Z +" << mesaOffsetCm << " cm para carpeta: " << parentFolder << G4endl;
            }

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





























    
//?//?//?//?//?//?//?//?//?//?//?//?//?//?
//* FANTOMA  
//?//?//?//?//?//?//?//?//?//?//?//?//?//?

 
 
    G4double PhantomX = 50*cm;
    G4double PhantomY = 50*cm;
    G4double PhantomZ =50*cm;
    G4Material* phantomMaterial =nist->FindOrBuildMaterial("G4_WATER");
 
 
// Crear el volumen del fantoma principal de PMMA
G4Box* solidPhantom = new G4Box("solidPhantom", PhantomX/2, PhantomY/2, PhantomZ/2);
G4VSolid* solidPhantomWithHoles = solidPhantom; // Inicialmente, es el bloque completo

  
G4LogicalVolume* logicPhantom = new G4LogicalVolume(solidPhantom, phantomMaterial, "logicalPhantom");

// Rotación y desplazamiento: rotar 180 grados en X y desplazar en Z (100 mm + media altura del fantoma)
G4RotationMatrix* rotPhantom = new G4RotationMatrix();
rotPhantom->rotateX(180.0*deg);
G4double zShift = 100.0*cm + PhantomZ/2.0;

// Colocar el fantoma en el mundo con rotación y desplazamiento en Z
G4VPhysicalVolume* physPhantom = new G4PVPlacement(rotPhantom,
                                                    G4ThreeVector(0, 0, zShift),
                                                    logicPhantom,
                                                    "physPhantom",
                                                    logicWorld,
                                                    false,
                                                    0,
                                                    true);
    



//& Número de detectores en cada fila y columna d
G4int numDetectors = 41; //61; // 23;  
  
G4double separation = 0.5*cm;  
G4int numLayers = 5; //7
G4double layerSeparation =  0.5*cm; 
G4double profundidadGrupo= 0.3*cm;


  

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^v
//^^Detector   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^v
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^v
//*-INTERNOS-

 G4double SDiodeX = 0.5*cm;
 G4double SDiodeY = 0.5*cm;
 G4double SDiodeZ = 0.5*cm;   
 
// Crear el volumen sólido del detector
  solidDetector = new G4Box("solidDetector", SDiodeX/2, SDiodeY/2, SDiodeZ/2);

// Crear el volumen lógico del detector
  logicVoxel = new G4LogicalVolume(solidDetector,phantomMaterial, "logicDetector");
    G4int ndet =1;
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) {
    for(G4int i = 0; i < numDetectors; i++) {
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos = -0.5*(numDetectors-1)*separation + i*separation;
            G4double yPos = -0.5*(numDetectors-1)*separation + j*separation;
            G4double zPos = +PhantomZ/2-profundidadGrupo -layer * (layerSeparation)   ;
            
            G4VPhysicalVolume *physDetector = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logicVoxel,
                "physDetector",
                logicPhantom,
                false,
                ndet
                
            );
            ndet++;
        } 
}
}
  
 
 

 








//& Número de detectores en cada fila y columna d  
 numLayers = 6; //7
 layerSeparation = 5*cm; 
  profundidadGrupo= 5*cm;


 

 
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) {
    for(G4int i = 0; i < numDetectors; i++) {
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos = -0.5*(numDetectors-1)*separation + i*separation;
            G4double yPos = -0.5*(numDetectors-1)*separation + j*separation;
            G4double zPos = +PhantomZ/2-profundidadGrupo -layer * (layerSeparation)   ;
            
            G4VPhysicalVolume *physDetector = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logicVoxel,
                "physDetector",
                logicPhantom,
                false,
                ndet
                
            );
            ndet++;
        } 
}
}
  
 
 


    return physWorld;
}



void MyDetectorConstruction::ConstructSDandField() {
    auto sensDet = new MySensitiveDetector("SensitiveDetector", "TrackerHitsCollection");
    G4SDManager::GetSDMpointer()->AddNewDetector(sensDet);

    if (logicVoxel != nullptr) {
        logicVoxel->SetSensitiveDetector(sensDet);
    }
} 