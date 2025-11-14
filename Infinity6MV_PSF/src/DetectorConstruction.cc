#include "DetectorConstruction.hh"
#include "G4SubtractionSolid.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Region.hh"
#include "G4Colour.hh"
#include "G4Polycone.hh"
#include "G4Trd.hh"
#include "G4Cons.hh"
#include "G4Trap.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4Region.hh"
#include "G4UnionSolid.hh"
#include "SensitiveDetector.hh"

#include "G4IntersectionSolid.hh"
MyDetectorConstruction::MyDetectorConstruction() {}
MyDetectorConstruction::~MyDetectorConstruction() {}
 
G4VPhysicalVolume* MyDetectorConstruction::Construct() {
    G4NistManager* nist = G4NistManager::Instance();
 


//============================================================
//                       MATERIALES
//============================================================

// ---------- Materiales NIST ----------
G4Material* columnaMat       = nist->FindOrBuildMaterial("G4_Galactic");
G4Material* al               = nist->FindOrBuildMaterial("G4_Al");
G4Material* aluminum         = nist->FindOrBuildMaterial("G4_Al");
G4Material* PMMA             = nist->FindOrBuildMaterial("G4_PLEXIGLASS"); 
G4Material* cobre            = nist->FindOrBuildMaterial("G4_Cu");
G4Material* Sii              = nist->FindOrBuildMaterial("G4_Si");
G4Material* tungsten         = nist->FindOrBuildMaterial("G4_W");
G4Material* revestimientoMaterial = nist->FindOrBuildMaterial("G4_W");
G4Material* worldMat         = nist->FindOrBuildMaterial("G4_AIR");
G4Material* stainlessSteel   = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");

// ---------- Elementos básicos ----------
G4Element* H   = nist->FindOrBuildElement("H");   // Hidrógeno
G4Element* C   = nist->FindOrBuildElement("C");   // Carbono
G4Element* O   = nist->FindOrBuildElement("O");   // Oxígeno
G4Element* Si  = nist->FindOrBuildElement("Si");  // Silicio
G4Element* Zn  = nist->FindOrBuildElement("Zn");  // Zinc
G4Element* B   = nist->FindOrBuildElement("B");   // Boro

// ---------- Elementos para acero inoxidable ----------
G4Element* Fe  = new G4Element("Iron",      "Fe", 26, 55.85 * g / mole);
G4Element* Cr  = new G4Element("Chromium",  "Cr", 24, 52.00 * g / mole);
G4Element* Ni  = new G4Element("Nickel",    "Ni", 28, 58.69 * g / mole);
G4Element* Mn  = new G4Element("Manganese", "Mn", 25, 54.94 * g / mole);

// ---------- Elementos especiales ----------
G4Element* elRe = new G4Element("Rhenium", "Re", 75., 186.207*g/mole);

// ---------- Compuestos y aleaciones personalizados ----------

// BaTiO3
G4Material* BaTiO3 = new G4Material("BaTiO3", 6.02 * g/cm3, 3);
BaTiO3->AddElement(nist->FindOrBuildElement("Ba"), 1);
BaTiO3->AddElement(nist->FindOrBuildElement("Ti"), 1);
BaTiO3->AddElement(nist->FindOrBuildElement("O"),  3);

// Epoxy
G4double densityEpoxy = 1.2 * g/cm3;  // según bibliografía
G4Material* epoxy = new G4Material("Epoxy", densityEpoxy, 3);
epoxy->AddElement(H, 0.077); // 7.7% Hidrógeno
epoxy->AddElement(C, 0.533); // 53.3% Carbono
epoxy->AddElement(O, 0.390); // 39.0% Oxígeno

// Quartz (SiO2) → MoldCompoundBlack
G4double density = 2.196 * g/cm3;
G4Element* elSi  = nist->FindOrBuildElement("Si");
G4Element* elO   = nist->FindOrBuildElement("O");

G4Material* quartz = new G4Material("quartz", density, 2);
quartz->AddElement(elSi, 1);
quartz->AddElement(elO,  2);

G4Material* MoldCompoundBlack = quartz;

// FR4 (Glass + Epoxy)
density = 1.86 * g/cm3;
G4Material* mat_FR4 = new G4Material("mat_FR4", density, 2);
mat_FR4->AddMaterial(quartz, 0.528);
mat_FR4->AddMaterial(epoxy,  0.472);

// W-Re Alloy (90% W, 10% Re)
G4Material* wre_alloy = new G4Material("WRe_Alloy", 19.25 * g/cm3, 2);
wre_alloy->AddMaterial(tungsten, 90. * perCent);
wre_alloy->AddElement(elRe,      10. * perCent);

// MLC Tungsten Alloy (95% W, 3.75% Ni, 1.25% Fe)
density = 18.0 * g/cm3;
G4Material* matTungstenMLC = new G4Material("MLCTungstenAlloy", density, 3);
matTungstenMLC->AddElement(nist->FindOrBuildElement("W"),  95.0 * perCent);
matTungstenMLC->AddElement(nist->FindOrBuildElement("Ni"), 3.75 * perCent);
matTungstenMLC->AddElement(nist->FindOrBuildElement("Fe"), 1.25 * perCent);

// ---------- Control ----------
G4bool checkOverlaps = true;



 


//?============================================================
//?                     GEOMETRÍA GEANT4 (CORREGIDO)
//?============================================================

//*========= 1. DEFINICIÓN DE PARÁMETROS GLOBALES =========
G4double SSDValue = 100 * cm; // Definido UNA sola vez aquí
G4double target_thickness = 1.1 * cm;

 
//*========= 2. MUNDO =========
G4double worldX = 13.31 * m;
G4double worldY = 13.31 * m;
G4double worldZ = 13.01 * m;

solidWorld = new G4Box("World", worldX / 2, worldY / 2, worldZ / 2);
logicWorld = new G4LogicalVolume(solidWorld, worldMat, "World");
physWorld  = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, true);

G4VisAttributes* worldVis = new G4VisAttributes(G4Colour(1.0,1.0,1.0,0.05));
worldVis->SetForceSolid(true);
logicWorld->SetVisAttributes(worldVis);

//*========= 3. BLOQUE CABEZAL (HEAD) Y OFFSET =========
// Dimensiones del contenedor de aire
G4double headSizeXY = 140.0 * cm;
G4double headSizeZ  = 52.0 * cm;
G4double headHalfZ  = headSizeZ / 2.0;

G4Box* solidHead = new G4Box("HeadBlock", headSizeXY/2, headSizeXY/2, headHalfZ);
G4LogicalVolume* logicHead = new G4LogicalVolume(solidHead, worldMat, "LogicHead"); // Aire

// Calculo del Offset para centrar el Target en Z=0 del mundo final
// La cara superior del target original está en: SSDValue + target_thickness/2
// Queremos que esa cara coincida con la cara +Z del bloque (headHalfZ)
G4double originalTargetTopZ = SSDValue + target_thickness / 2.0;
G4double zShift = originalTargetTopZ - headHalfZ;
G4ThreeVector geometryOffset(0, 0, zShift);

G4VisAttributes* headVis = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0, 0.1));
headVis->SetVisibility(true);
logicHead->SetVisAttributes(headVis);

//*========= 4. COMPONENTES (Target, Filtros, etc) =========

// --- Target ---
G4double target_radius = 1.0 * cm;
G4Tubs* solidTarget = new G4Tubs("Target", 0., target_radius, target_thickness / 2., 0. * deg, 360. * deg);
G4LogicalVolume* logicTarget = new G4LogicalVolume(solidTarget, wre_alloy, "TargetLV");
G4Region* TargetRegion = new G4Region("TargetRegion");
TargetRegion->AddRootLogicalVolume(logicTarget);

G4ThreeVector posTargetOrig(0, 0, SSDValue - 11 * mm);
new G4PVPlacement(nullptr, posTargetOrig - geometryOffset, logicTarget, "Target", logicHead, false, 0, true);

G4VisAttributes* targetVis = new G4VisAttributes(G4Colour::Gray());
targetVis->SetVisibility(true);
logicTarget->SetVisAttributes(targetVis);

// --- Disco y Filtro ---
G4double discRadius    = 4.0 * cm;
G4double discThickness = 0.5 * cm;
G4double filterHeight  = 2.77 * cm;
G4double filterZ = (SSDValue - 15.9 * cm);
G4double discZ   = filterZ - (discThickness / 2 + filterHeight / 2);

G4Tubs* solidDisc = new G4Tubs("AlDisc", 0, discRadius, discThickness / 2, 0, 360 * deg);
G4Cons* solidFilter = new G4Cons("FlatteningFilter", 0, 4.0*cm, 0, 0, filterHeight / 2, 0, 360 * deg);

G4LogicalVolume* logicDisc   = new G4LogicalVolume(solidDisc, aluminum, "AlDiscLV");
G4LogicalVolume* logicFilter = new G4LogicalVolume(solidFilter, stainlessSteel, "FlatteningFilterLV");

new G4PVPlacement(0, G4ThreeVector(0, 0, discZ) - geometryOffset,   logicDisc,   "AlDisc", logicHead, false, 0, true);
new G4PVPlacement(0, G4ThreeVector(0, 0, filterZ) - geometryOffset, logicFilter, "FlatteningFilter", logicHead, false, 0, true);

G4VisAttributes* discVis = new G4VisAttributes(G4Colour::Blue());
discVis->SetForceSolid(true); logicDisc->SetVisAttributes(discVis);
G4VisAttributes* filterVis = new G4VisAttributes(G4Colour::Cyan());
filterVis->SetForceSolid(true); logicFilter->SetVisAttributes(filterVis);

// --- Colimador Primario ---
G4double cone_height = 11.2 * cm;
G4Cons* cone_hole = new G4Cons("ConeHole", 0, 5.0*cm, 0, 1.2*cm, cone_height/2, 0, 360*deg);
G4Tubs* collimator_body = new G4Tubs("CollimatorBody", 0, 7.0*cm, cone_height/2, 0, 360*deg);
G4SubtractionSolid* collimator_hollowed = new G4SubtractionSolid("CollimatorWithHole", collimator_body, cone_hole);
G4LogicalVolume* collimator_logical = new G4LogicalVolume(collimator_hollowed, tungsten, "CollimatorLogical");

G4double z_position_collimator = cone_height / 2;
G4ThreeVector posCollOrig(0, 0, SSDValue - z_position_collimator);
new G4PVPlacement(nullptr, posCollOrig - geometryOffset, collimator_logical, "PrimaryCollimator", logicHead, false, 0, true);

G4VisAttributes* collimatorVis = new G4VisAttributes(G4Colour::Magenta());
collimatorVis->SetForceSolid(true); collimator_logical->SetVisAttributes(collimatorVis);

// --- Revestimiento y Cilindro ---
G4Tubs* revestimiento = new G4Tubs("Revestimiento", 63.26*cm, 68.34*cm, 59.625*cm/2, 0, 360*deg);
G4LogicalVolume* logicRevestimiento = new G4LogicalVolume(revestimiento, revestimientoMaterial, "LogicRevestimiento");
G4ThreeVector posRevOrig(0, 0, SSDValue - 20.6325 * cm);
new G4PVPlacement(0, posRevOrig - geometryOffset, logicRevestimiento, "Revestimiento", logicHead, false, 0, checkOverlaps);

G4Tubs* cilindroHueco = new G4Tubs("CilindroHueco", 19*cm, 68.34*cm, 1.0*cm, 0, 360*deg); // h=2cm/2
G4LogicalVolume* logicCilindroHueco = new G4LogicalVolume(cilindroHueco, revestimientoMaterial, "LogicCilindroHueco");
G4double zPosCilindro = 20.6325 * cm + 59.625 * cm / 2 + 1.0 * cm;
new G4PVPlacement(0, G4ThreeVector(0, 0, SSDValue - zPosCilindro) - geometryOffset, logicCilindroHueco, "CilindroHueco", logicHead, false, 0, checkOverlaps);

//*========= 5. MLC =========
G4double fieldSize  = 10.0 * cm;
G4double L          = fieldSize / 2.0;
G4double leafWidthY = 5.0 * mm;
G4double leafZPosition = SSDValue - 30.9 * cm;
int numLeaves = 80;
G4double interleafGap = 0.09 * mm;

// Definición sólida MLC
G4Box* blockSolid = new G4Box("Block", 17.0*cm, 4.5*cm, 4.5*cm);
G4Tubs* diskSolid = new G4Tubs("Disk", 0, 17.0*cm, 2.5*mm, 0, 360.0*deg);
G4Transform3D* transformDisk = new G4Transform3D(G4RotationMatrix(), G4ThreeVector(0, -7.5*mm, 0));
G4IntersectionSolid* intersectSolid = new G4IntersectionSolid("LeafSolid", blockSolid, diskSolid, *transformDisk);
G4LogicalVolume* leafLogic = new G4LogicalVolume(intersectSolid, matTungstenMLC, "LeafLogic");
G4VisAttributes* leafVis = new G4VisAttributes(G4Colour::Red());
leafVis->SetForceSolid(true); leafLogic->SetVisAttributes(leafVis);

auto rotX = new G4RotationMatrix();
rotX->rotateX(90.0 * deg);
rotX->rotateY(9 * mrad);

G4double mlcHalfOpening = L * (SSDValue - leafZPosition) / SSDValue;
G4double pitch = leafWidthY + interleafGap;
G4double yHalfOpening = L * (SSDValue - leafZPosition) / SSDValue;

std::vector<G4double> leafXPShift(numLeaves, 0.0);
std::vector<G4double> leafXNShift(numLeaves, 0.0);
int numActiveLeaves = 0; // VARIABLE DECLARADA AQUI
int extraLeaves = 1;

for (int i = 0; i < numLeaves; ++i) {
    G4double yLeaf = -((numLeaves - 1) * pitch / 2.0) + i * pitch;
    if (std::abs(yLeaf) <= yHalfOpening + extraLeaves * pitch) {
        numActiveLeaves++; // USO DE LA VARIABLE
        leafXPShift[i] = mlcHalfOpening;
        leafXNShift[i] = -mlcHalfOpening;
    }
}

for (int i = 0; i < numLeaves; ++i) {
    G4double y = -((numLeaves - 1) * pitch / 2.0) + i * pitch;
    // +X
    new G4PVPlacement(rotX, G4ThreeVector(17.0*cm + leafXPShift[i], y, leafZPosition) - geometryOffset,
        leafLogic, ("LeafXP_"+std::to_string(i)).c_str(), logicHead, false, i, true);
    // -X
    new G4PVPlacement(rotX, G4ThreeVector(-17.0*cm + leafXNShift[i], y, leafZPosition) - geometryOffset,
        leafLogic, ("LeafXN_"+std::to_string(i)).c_str(), logicHead, false, i + numLeaves, true);
}

//*========= 6. JAWS (CORREGIDO) =========
G4double jawLengthX   = 20.0 * cm;
G4double jawWidthY    = 20.0 * cm;
G4double jawHeightZ   = 7.80 * cm;
G4double jawZPosition = SSDValue - 43.2 * cm;
G4Material* jawMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_W");

// 1. Crear Sólidos y Lógicos PRIMERO
G4Box* solidJawYP = new G4Box("JawYP", jawLengthX/2, jawWidthY/2, jawHeightZ/2);
G4LogicalVolume* logicJawYP = new G4LogicalVolume(solidJawYP, jawMaterial, "LogicJawYP");

G4Box* solidJawYN = new G4Box("JawYN", jawLengthX/2, jawWidthY/2, jawHeightZ/2);
G4LogicalVolume* logicJawYN = new G4LogicalVolume(solidJawYN, jawMaterial, "LogicJawYN");

G4VisAttributes* jawVis = new G4VisAttributes(G4Colour(1.0, 0.5, 0.0));
jawVis->SetForceSolid(true);
logicJawYP->SetVisAttributes(jawVis);
logicJawYN->SetVisAttributes(jawVis);

// 2. Cálculos de posición
G4double a = jawHeightZ;
G4double b = jawWidthY;
G4double z0 = SSDValue - jawZPosition;
G4double D1 = 0.5 * std::sqrt(a*a + b*b);
G4double theta = (3.14159265/2) - std::atan(SSDValue/L);
G4double beta  = std::atan(SSDValue/L);
G4double alpha = std::atan(a/b);
G4double xj = z0*L/SSDValue + D1*std::cos(alpha+theta) + D1*std::sin(alpha+theta)/std::tan(beta);

G4RotationMatrix* jawRotationP = new G4RotationMatrix(); jawRotationP->rotateX(-theta*rad);
G4RotationMatrix* jawRotationN = new G4RotationMatrix(); jawRotationN->rotateX(theta*rad);

// 3. Colocación
new G4PVPlacement(jawRotationP, G4ThreeVector(0, xj, jawZPosition) - geometryOffset, logicJawYP, "PhysJawYP", logicHead, false, 0, true);
new G4PVPlacement(jawRotationN, G4ThreeVector(0, -xj, jawZPosition) - geometryOffset, logicJawYN, "PhysJawYN", logicHead, false, 1, true);


//*========= 7. COLOCACIÓN FINAL DEL CABEZAL =========
// Rotar 180 grados y posicionar en el mundo
G4RotationMatrix* headRotation = new G4RotationMatrix();
headRotation->rotateY(180.0 * deg);

// Subir el centro del bloque para que el TopTarget (local +Z) quede en World Z=0
G4ThreeVector headPosition(0, 0, headHalfZ);

new G4PVPlacement(headRotation, headPosition, logicHead, "HeadPhys", logicWorld, false, 0, true);





G4double M_PI =3.14159265358979323846;



//*========= VERIFICACION =========
G4cout << "Field Size: " << fieldSize/cm << " cm\n";
G4cout << "Jaw displacement xj: " << xj/cm << " cm\n";
G4cout << "Jaw rotation theta: " << theta*180/M_PI << " deg\n";
G4cout << "Number of active MLC leaves: " << numActiveLeaves << " out of " << numLeaves << G4endl;

for(int i=0;i<numLeaves;i++){
    if(leafXPShift[i]!=0.0 || leafXNShift[i]!=0.0){
        G4cout << "Leaf " << i << ": XP shift = " << leafXPShift[i]/cm
               << " cm, XN shift = " << leafXNShift[i]/cm << " cm" << G4endl;
    }
}

 
//? Dimensiones y material del fantoma
    G4double PhantomX = 30*cm;
    G4double PhantomY = 30*cm;
    G4double PhantomZ =30*cm;
    G4Material* phantomMaterial =nist->FindOrBuildMaterial("G4_WATER");
 
 
// Crear el volumen del fantoma principal de PMMA
G4Box* solidPhantom = new G4Box("solidPhantom", PhantomX/2, PhantomY/2, PhantomZ/2);
G4VSolid* solidPhantomWithHoles = solidPhantom; // Inicialmente, es el bloque completo

  
G4LogicalVolume* logicPhantom = new G4LogicalVolume(solidPhantom, phantomMaterial, "logicalPhantom");

// Colocar el fantoma en el mundo
G4VPhysicalVolume* physPhantom = new G4PVPlacement(0, G4ThreeVector(0, 0, 100*cm+PhantomZ/2),  
                                                    logicPhantom , 
                                                    "physPhantom", 
                                                    logicWorld, 
                                                    false, 
                                                    0, 
                                                    true);
    



 


 



 





// //TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// ////*   1.5 - 9.5 cm !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// //TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// //TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 


//& Número de detectores en cada fila y columna d
G4int numDetectors = 3;// 61; // 23;  

G4int outerNumDetectors =0  ; //31;
//& Separación
G4double separation = 0.5*cm;
G4double outerSeparation = 1.4 * cm; 
//& Capas
G4int numLayers = 20; //7
G4double layerSeparation =  0.5*cm;
//& PROFUNDIDAD PRIMERA CAPA
G4double profundidadGrupo= 1.5*cm;


 


//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*  DIODOS    //*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
G4double coberturaX = 0.450 * cm;
G4double coberturaY = 0.400 * cm;
G4double coberturaZ = 0.120 * cm;
 
 G4double SDiodeX = 0.5*cm;
 G4double SDiodeY = 0.5*cm;
 G4double SDiodeZ = 0.5*cm;
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^v
//^^cobertura    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^v
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^v
// Creación del cubo exterior (coberturadiodo)
 G4int ndeth =1;
// Crear el volumen exterior (coberturadiodo)
G4Box* solidCoberturaDiode = new G4Box("CoberturaDiode", coberturaX/2, coberturaY/2, coberturaZ/2);

// Crear el volumen del hueco (SDiode)
G4Box* solidSDiodeh = new G4Box("SDiode", SDiodeX/2, SDiodeY/2, SDiodeZ/2);

// Posición del hueco dentro de la cobertura (ajusta según sea necesario)
G4ThreeVector posicionHueco(0, 0, 0);
 

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^v
//^^SENSITIVO    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^v
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^v
//*-INTERNOS-

// Crear el volumen sólido del detector
  solidDetector = new G4Box("solidDetector", SDiodeX/2, SDiodeY/2, SDiodeZ/2);

// Crear el volumen lógico del detector
  logicmosfet = new G4LogicalVolume(solidDetector,Sii, "logicDetector");
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
                logicmosfet,
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

    if (logicmosfet != nullptr) {
        logicmosfet->SetSensitiveDetector(sensDet);
    }
}
