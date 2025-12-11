// ====================================================================
// 1. MACRO CRÍTICA: Primera línea absoluta.
// ====================================================================
#define USE_CADMESH_TETGEN 1

#include "DetectorConstruction.hh"
#include "CADMesh.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4AssemblyVolume.hh"
#include "G4VisAttributes.hh"

// Librería para leer carpetas (requiere C++17)
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
    
    // Usaremos Plomo (Pb) para el MLC (Multileaf Collimator)
    G4Material* mlcMat = nist->FindOrBuildMaterial("G4_Pb"); 

    // 2. Definir World
    G4Box* solidWorld = new G4Box("World", 2.*m, 2.*m, 2.*m);
    G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, worldMat, "World");
    G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, true);

    // 3. Bucle para cargar TODOS los archivos PLY de la carpeta
    std::string folderPath = "C:\\Users\\renzo\\Downloads\\MLC_HALCYON_STL\\ply";
    
    G4cout << "--- INICIANDO CARGA MASIVA DE PLY ---" << G4endl;

    // Iterar sobre cada archivo en el directorio
    for (const auto& entry : fs::directory_iterator(folderPath)) {
        
        // Verificar que sea un archivo y que tenga extensión .ply
        if (entry.is_regular_file() && entry.path().extension() == ".ply") {
            
            // Convertir ruta a string
            std::string filePath = entry.path().string();
            G4cout << "Cargando: " << filePath << G4endl;

            // a) Cargar Malla
            auto mesh = CADMesh::TetrahedralMesh::FromPLY(filePath);
            
            // b) Configuración
            // NOTA: SetScale(1.0) mantiene el tamaño original del PLY (mm).
            // Si usas 10*cm, estás escalando x100 (muy grande). Ajusta según necesidad.
            mesh->SetScale(1.0); 
            
            // Offset: Asumimos que la posición relativa (X,Y) viene correcta desde Blender.
            // Solo desplazamos en Z si es necesario.
            mesh->SetOffset(G4ThreeVector(0, 0, 0));
            
            // c) Material
            mesh->SetMaterial(mlcMat);

            // d) Obtener Ensamblaje y Colocar
            G4AssemblyVolume* assembly = mesh->GetAssembly();

            // Variables para placement
            G4ThreeVector position = G4ThreeVector(0, 0, 0);
            G4RotationMatrix* rotation = new G4RotationMatrix();
            
            // Imprimir en el mundo
            assembly->MakeImprint(logicWorld, position, rotation);
        }
    }
    
    G4cout << "--- CARGA COMPLETA ---" << G4endl;

    return physWorld;
}