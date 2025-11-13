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
//?                       GEOMETRÍA GEANT4
//?============================================================

//*========= Mundo =========
G4double worldX = 13.31 * m;
G4double worldY = 13.31 * m;
G4double worldZ = 13.01 * m;

solidWorld = new G4Box("World", worldX / 2, worldY / 2, worldZ / 2);
logicWorld = new G4LogicalVolume(solidWorld, worldMat, "World");
physWorld  = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, true);

G4VisAttributes* worldVis = new G4VisAttributes(G4Colour(1.0,1.0,1.0,0.05)); // casi transparente
worldVis->SetForceSolid(true);
logicWorld->SetVisAttributes(worldVis);


//*========= Target =========
G4double target_radius    = 1.0 * cm;   // Radio del disco
G4double target_thickness = 0.09 * cm;  // Grosor en dirección del haz
G4double SSDValue         =100 * cm;

G4Tubs* solidTarget = new G4Tubs("Target",
    0., target_radius,
    target_thickness / 2.,
    0. * deg, 360. * deg);

G4LogicalVolume* logicTarget = new G4LogicalVolume(solidTarget, wre_alloy, "TargetLV");

G4Region* TargetRegion = new G4Region("TargetRegion");
TargetRegion->AddRootLogicalVolume(logicTarget);

G4ThreeVector posTarget(0, 0, SSDValue - 11 * mm);
new G4PVPlacement(nullptr, posTarget, logicTarget, "Target", logicWorld, false, 0, true);

G4VisAttributes* targetVis = new G4VisAttributes(G4Colour::Gray());
targetVis->SetVisibility(true);
logicTarget->SetVisAttributes(targetVis); 

// //*========= Disco y filtro (invertidos en -Z) =========
G4double discRadius    = 4.0 * cm;
G4double discThickness = 0.5 * cm;

G4double filterHeight = 2.77 * cm;
G4double filterRmin1  = 0.0 * cm;
G4double filterRmax1  = 4.0 * cm;   // ahora va abajo
G4double filterRmin2  = 0.0 * cm;
G4double filterRmax2  = 0.0 * cm;   // ahora va arriba

// Posiciones invertidas en Z
G4double filterZ = (SSDValue - 15.4 * cm);
G4double discZ   = filterZ - (discThickness / 2 + filterHeight / 2);

G4Tubs* solidDisc = new G4Tubs("AlDisc", 0, discRadius, discThickness / 2, 0, 360 * deg);

G4Cons* solidFilter = new G4Cons("FlatteningFilter",
    filterRmin1, filterRmax1,   // base en -Z
    filterRmin2, filterRmax2,   // punta en +Z
    filterHeight / 2,
    0, 360 * deg);

G4LogicalVolume* logicDisc   = new G4LogicalVolume(solidDisc, aluminum, "AlDiscLV");
G4LogicalVolume* logicFilter = new G4LogicalVolume(solidFilter, stainlessSteel, "FlatteningFilterLV");

new G4PVPlacement(0, G4ThreeVector(0, 0, discZ),   logicDisc,   "AlDisc",           logicWorld, false, 0, true);
new G4PVPlacement(0, G4ThreeVector(0, 0, filterZ), logicFilter, "FlatteningFilter", logicWorld, false, 0, true);


G4VisAttributes* discVis   = new G4VisAttributes(G4Colour::Blue());
discVis->SetForceSolid(true);
logicDisc->SetVisAttributes(discVis);

G4VisAttributes* filterVis = new G4VisAttributes(G4Colour::Cyan());
filterVis->SetForceSolid(true);
logicFilter->SetVisAttributes(filterVis);


//*========= Colimador primario =========
G4double cone_height      = 11.2 * cm;
G4double cone_base_radius = 5.0 * cm;
G4double cone_top_radius  = 1.2 * cm;

G4double body_height       = cone_height;
G4double body_inner_radius = 0.0 * cm;
G4double body_outer_radius = 7.0 * cm;

G4Cons* cone_hole = new G4Cons("ConeHole",
    0, cone_base_radius,
    0, cone_top_radius,
    cone_height / 2,
    0, 360 * deg);

G4Tubs* collimator_body = new G4Tubs("CollimatorBody",
    body_inner_radius, body_outer_radius,
    body_height / 2,
    0, 360 * deg);

G4SubtractionSolid* collimator_hollowed = new G4SubtractionSolid(
    "CollimatorWithHole", collimator_body, cone_hole,
    nullptr, G4ThreeVector(0, 0, 0));

G4LogicalVolume* collimator_logical = new G4LogicalVolume(
    collimator_hollowed, tungsten, "CollimatorLogical");

G4double z_position_collimator = cone_height / 2;
new G4PVPlacement(nullptr,
    G4ThreeVector(0, 0, SSDValue - z_position_collimator),
    collimator_logical, "PrimaryCollimator", logicWorld, false, 0, true);

G4VisAttributes* visAttr = new G4VisAttributes(G4Colour::Grey());
visAttr->SetForceSolid(true);
collimator_logical->SetVisAttributes(visAttr);


G4VisAttributes* collimatorVis = new G4VisAttributes(G4Colour::Magenta());
collimatorVis->SetForceSolid(true);
collimator_logical->SetVisAttributes(collimatorVis);

//*========= Revestimiento =========
G4Tubs* revestimiento = new G4Tubs("Revestimiento",
    63.26 * cm, 68.34 * cm,
    59.625 * cm / 2, 0, 360 * deg);

G4LogicalVolume* logicRevestimiento = new G4LogicalVolume(
    revestimiento, revestimientoMaterial, "LogicRevestimiento");

new G4PVPlacement(0, G4ThreeVector(0, 0, SSDValue - 20.6325 * cm),
    logicRevestimiento, "Revestimiento", logicWorld, false, 0, checkOverlaps);



G4VisAttributes* revestVis = new G4VisAttributes(G4Colour::Green());
revestVis->SetForceSolid(true);
logicRevestimiento->SetVisAttributes(revestVis);

//*========= Cilindro hueco superior =========
G4double rInterno = 19 * cm;
G4double rExterno = 68.34 * cm;
G4double hCilindro = 2.0 * cm;

G4Tubs* cilindroHueco = new G4Tubs("CilindroHueco",
    rInterno, rExterno,
    hCilindro / 2, 0, 360 * deg);

G4LogicalVolume* logicCilindroHueco = new G4LogicalVolume(
    cilindroHueco, revestimientoMaterial, "LogicCilindroHueco");

G4double zPosCilindro = 20.6325 * cm + 59.625 * cm / 2 + hCilindro / 2;

new G4PVPlacement(0, G4ThreeVector(0, 0, SSDValue - zPosCilindro),
    logicCilindroHueco, "CilindroHueco", logicWorld, false, 0, checkOverlaps);


G4VisAttributes* cilindroVis = new G4VisAttributes(G4Colour::Yellow());
cilindroVis->SetForceSolid(true);
logicCilindroHueco->SetVisAttributes(cilindroVis);

 
//*========= PARAMETROS DEL CAMPO =========
G4double fieldSize  = 5.0 * cm;      // <-- Cambiar campo (ej. 10, 20, etc.)
G4double L          = fieldSize / 2.0; // Semi-campo en superficie

//*========= MLC =========
G4double leafWidthY     = 5.0 * mm;
G4double leafHeightZ    = 90.0 * mm;
G4double leafZPosition  = SSDValue - 30.9 * cm;
int numLeaves           = 80;
G4double interleafGap   = 0.09 * mm;

G4double blockX = 34.0 * cm;
G4double blockY = 9.0 * cm;
G4double blockZ = 9.0 * cm;
auto blockSolid = new G4Box("Block", blockX / 2, blockY / 2, blockZ / 2);

G4double diskRadiusMLC = 17.0 * cm;
G4double diskHeight    = 5.0 * mm;
auto diskSolid = new G4Tubs("Disk", 0, diskRadiusMLC, diskHeight / 2, 0, 360.0 * deg);

G4ThreeVector diskOffset(0, -7.5 * mm, 0);
auto transformDisk = new G4Transform3D(G4RotationMatrix(), diskOffset);

auto intersectSolid = new G4IntersectionSolid("LeafSolid", blockSolid, diskSolid, *transformDisk);
auto leafLogic = new G4LogicalVolume(intersectSolid, matTungstenMLC, "LeafLogic");

auto rotX = new G4RotationMatrix();
rotX->rotateX(90.0 * deg);
rotX->rotateY(9 * mrad);

// === Calculos automáticos ===
G4double mlcHalfOpening = L * (SSDValue - leafZPosition) / SSDValue; // desplazamiento en X
G4double pitch = leafWidthY + interleafGap;
G4double yHalfOpening = L * (SSDValue - leafZPosition) / SSDValue;   // límite en Y para hojas activas

std::vector<G4double> leafXPShift(numLeaves, 0.0);
std::vector<G4double> leafXNShift(numLeaves, 0.0);
int numActiveLeaves = 0;

// Selección automática de hojas que se moverán (+2 hojas extra en cada extremo)
// Selección automática de hojas que se moverán (+extraLeaves en cada extremo)
int extraLeaves = 1; 
for (int i = 0; i < numLeaves; ++i) {
    G4double yLeaf = -((numLeaves - 1) * pitch / 2.0) + i * pitch;

    if (std::abs(yLeaf) <= yHalfOpening + extraLeaves * pitch) {
        numActiveLeaves++;

        // 👉 todas las hojas se abren igual (sin apertura parcial)
        G4double frac = 1.0;

        leafXPShift[i] = mlcHalfOpening * frac;
        leafXNShift[i] = -mlcHalfOpening * frac;
    }
    else {
        leafXPShift[i] = 0.0;
        leafXNShift[i] = 0.0;
    }
}


// Hojas +X
for (int i = 0; i < numLeaves; ++i) {
    G4double y = -((numLeaves - 1) * pitch / 2.0) + i * pitch;
    G4double x = blockX / 2 + leafXPShift[i];
    G4ThreeVector pos(x, y, leafZPosition);

G4cout << "positivo " << y << " cm\n";
    new G4PVPlacement(rotX, pos, leafLogic, ("LeafXPPhys_" + std::to_string(i)).c_str(),
        logicWorld, false, i, true); 

     
    G4double xn = -blockX / 2 + leafXNShift[i];  
    G4ThreeVector posn(xn, y, leafZPosition);

G4cout << "negativo " << y << " cm\n";
    new G4PVPlacement(rotX, posn, leafLogic, ("LeafXNPhys_" + std::to_string(i)).c_str(),
        logicWorld, false, i + numLeaves, true);
    
}

auto leafVis = new G4VisAttributes(G4Colour::Red());
leafVis->SetForceSolid(true);
leafLogic->SetVisAttributes(leafVis);

//*========= Jaws =========
G4double jawLengthX   = 20.0 * cm;
G4double jawWidthY    = 20.0 * cm;
G4double jawHeightZ   = 7.80 * cm;
G4double jawZPosition = SSDValue - 43.2 * cm;


G4Material* jawMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_W");

// Calculo de desplazamiento y angulo usando geometria como tu script Python
G4double a = jawHeightZ;
G4double b = jawWidthY;
G4double z0 = SSDValue-jawZPosition;
G4double D1 = 0.5*std::sqrt(a*a + b*b);
G4double SSD = SSDValue;
G4double  M_PI = 3.14159265;
G4double theta = M_PI/2 - std::atan(SSD/L);    // inclinación
G4double beta  = std::atan(SSD/L);
G4double alpha = std::atan(a/b);

G4double xj = z0*L/SSD + D1*std::cos(alpha+theta) + D1*std::sin(alpha+theta)/std::tan(beta);

// Giro de las jaws
G4RotationMatrix* jawRotationP = new G4RotationMatrix();
G4RotationMatrix* jawRotationN = new G4RotationMatrix();
jawRotationP->rotateX(-theta*rad);
jawRotationN->rotateX(theta*rad);

// Jaw +Y
{
    G4ThreeVector pos(0, xj, jawZPosition);
    G4Box* solidJawYP = new G4Box("JawYP", jawLengthX/2, jawWidthY/2, jawHeightZ/2);
    G4LogicalVolume* logicJawYP = new G4LogicalVolume(solidJawYP,jawMaterial,"LogicJawYP");
    new G4PVPlacement(jawRotationP,pos,logicJawYP,"PhysJawYP",logicWorld,false,0,true);
    G4VisAttributes* jawYPVis = new G4VisAttributes(G4Colour(1.0,0.5,0.0));
    jawYPVis->SetForceSolid(true);
    logicJawYP->SetVisAttributes(jawYPVis);
}

// Jaw -Y
{
    G4ThreeVector pos(0, -xj, jawZPosition);
    G4Box* solidJawYN = new G4Box("JawYN", jawLengthX/2, jawWidthY/2, jawHeightZ/2);
    G4LogicalVolume* logicJawYN = new G4LogicalVolume(solidJawYN,jawMaterial,"LogicJawYN");
    new G4PVPlacement(jawRotationN,pos,logicJawYN,"PhysJawYN",logicWorld,false,1,true);
    G4VisAttributes* jawYNVis = new G4VisAttributes(G4Colour(1.0,0.5,0.0));
    jawYNVis->SetForceSolid(true);
    logicJawYN->SetVisAttributes(jawYNVis);
}

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

//?//?//?//?//?//?//?//?//?//?//?//?//?//?
//?WORLD
 
//? Dimensiones y material del fantoma
    G4double PhantomX = 40.2*cm;
    G4double PhantomY = 28*cm;
    G4double PhantomZ = 0.143*m;
    G4Material* phantomMaterial =nist->FindOrBuildMaterial("G4_PLEXIGLASS");

    // Dimensiones de los huecos de aire
    G4double airHoleX = 36.2*cm;
    G4double airHoleY = .2*m; 
    G4double airHoleZ = 0.5*cm;

// Profundidades donde se necesitan los huecos
std::vector<G4double> airHoleDepths = {1.5*cm, 3.5*cm, 5.5*cm, 7.5*cm, 9.5*cm};

// Crear el volumen del fantoma principal de PMMA
G4Box* solidPhantom = new G4Box("solidPhantom", PhantomX/2, PhantomY/2, PhantomZ/2);
G4VSolid* solidPhantomWithHoles = solidPhantom; // Inicialmente, es el bloque completo

// Crear y restar cada hueco
for (G4double depth : airHoleDepths) {
    G4ThreeVector holePosition(airHoleX/2-PhantomX/2+4*cm  , 0, +PhantomZ/2 - depth+(airHoleZ/2- 0.120*cm/2- 0.007*cm-0.16*cm- 0.007*cm) );// -medioDiodo - pad -bak- 
    G4Box* solidAirHole = new G4Box("solidAirHole", airHoleX/2, airHoleY/2, airHoleZ/2);
    solidPhantomWithHoles = new G4SubtractionSolid("solidPhantomWithHoles", solidPhantomWithHoles, solidAirHole, 0, holePosition);
}

// Crear el volumen lógico del fantoma con huecos
G4LogicalVolume* logiPhantom = new G4LogicalVolume(solidPhantomWithHoles, phantomMaterial, "logicalPhantom");

// Colocar el fantoma en el mundo
G4VPhysicalVolume* physPhantom = new G4PVPlacement(0, G4ThreeVector(PhantomX/2-14*cm, 0, -PhantomZ/2), logiPhantom, "physPhantom", logicWorld, false, 0, true);
    



 


 



 





// //TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// ////*   1.5 - 9.5 cm !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// //TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// //TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 


//& Número de detectores en cada fila y columna d
G4int numDetectors =  23; // 23;  

G4int outerNumDetectors =0  ; //31;
//& Separación
G4double separation = 0.8*cm;
G4double outerSeparation = 1.4 * cm; 
//& Capas
G4int numLayers = 5; //7
G4double layerSeparation =  2.*cm;
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
 
 G4double SDiodeX = 0.265*cm;
 G4double SDiodeY = 0.265*cm;
 G4double SDiodeZ = 0.06*cm;
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

// Crear el sólido final restando el hueco de la cobertura
G4SubtractionSolid* solidCoberturaDiodeConHueco = new G4SubtractionSolid("CoberturaDiodeConHueco", solidCoberturaDiode, solidSDiodeh, 0, posicionHueco);

// Crear el volumen lógico con el material epoxi
G4LogicalVolume* logicCoberturaDiodeConHueco = new G4LogicalVolume(solidCoberturaDiodeConHueco, epoxy, "CoberturaDiodeConHueco");


for(G4int layer = 0; layer < numLayers; layer++) {
    for(G4int i = 0; i < numDetectors; i++) {
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos = -0.5*(numDetectors-1)*separation + i*separation;
            G4double yPos = -0.5*(numDetectors-1)*separation + j*separation;
            G4double zPos = +PhantomZ/2-profundidadGrupo -layer * (layerSeparation) -PhantomZ/2 ;
            
            G4VPhysicalVolume *physCoberturaDiodeConHueco = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logicCoberturaDiodeConHueco,
                "logicCoberturaDiodeConHueco",
                logicWorld,
                false,
                ndeth
                
            );
            ndeth++;
        }
}
}
  
 

//*-EXTERNOS-
  
// Crear detectores externos
for(G4int layer = 0; layer < numLayers; layer++) {
    for(G4int i = 0; i < outerNumDetectors; i++) {
        for(G4int j = 0; j < outerNumDetectors; j++) {
            // Calcular la posición del detector externo
            G4double xPos = -0.5*(outerNumDetectors-1)* outerSeparation+ i* outerSeparation;
            G4double yPos = -0.5*(outerNumDetectors-1)* outerSeparation + j* outerSeparation;

            // Asegurarse de no superponer los detectores centrales
            G4double centralHalfWidth = 0.5 * (numDetectors-1) * separation+outerSeparation  ;
            if (std::abs(xPos) >= centralHalfWidth || std::abs(yPos) >= centralHalfWidth) {
                G4double zPos = +PhantomZ/2-profundidadGrupo - layer * layerSeparation-PhantomZ/2;

            G4VPhysicalVolume *physCoberturaDiodeConHueco = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logicCoberturaDiodeConHueco,
                "logicCoberturaDiodeConHueco",
                logicWorld,
                false,
                ndeth
                
            );
            ndeth++;
            }

            
        }
    }
} 



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
            G4double zPos = +PhantomZ/2-profundidadGrupo -layer * (layerSeparation)-PhantomZ/2  ;
            
            G4VPhysicalVolume *physDetector = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logicmosfet,
                "physDetector",
                logicWorld,
                false,
                ndet
                
            );
            ndet++;
        } 
}
}
  
 

//*-EXTERNOS-
  
// Crear detectores externos
for(G4int layer = 0; layer < numLayers; layer++) {
    for(G4int i = 0; i < outerNumDetectors; i++) {
        for(G4int j = 0; j < outerNumDetectors; j++) {
            // Calcular la posición del detector externo
            G4double xPos = -0.5*(outerNumDetectors-1)* outerSeparation+ i* outerSeparation;
            G4double yPos = -0.5*(outerNumDetectors-1)* outerSeparation + j* outerSeparation;

            // Asegurarse de no superponer los detectores centrales
            G4double centralHalfWidth = 0.5 * (numDetectors-1) * separation+outerSeparation  ;
            if (std::abs(xPos) >= centralHalfWidth || std::abs(yPos) >= centralHalfWidth) {
                G4double zPos = +PhantomZ/2-profundidadGrupo - layer * layerSeparation-PhantomZ/2;

            G4VPhysicalVolume *physDetector = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logicmosfet,
                "physDetector",
                 logicWorld,
                false,
                ndet
                
            );
            ndet++;
            }

            
        }
    }
} 






//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^v
//^^^^^CONEXIONES   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^v
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^v
 G4int ndethcn =1;
//*-COBREEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEINTERNOS-
 G4double cobreConectionSDiodeX = 0.12*cm;
 G4double cobreConectionSDiodeY = 0.12*cm;
 G4double cobreConectionSDiodeZ = 0.007*cm;
 
// Crear el volumen sólido del detector
   G4Box *  cobreConectionsolidDetector = new G4Box("cobreConectionsolidDetector", cobreConectionSDiodeX/2, cobreConectionSDiodeY/2, cobreConectionSDiodeZ/2);

// Crear el volumen lógico del detector
   G4LogicalVolume * logiccobreConectionDetector = new G4LogicalVolume(cobreConectionsolidDetector,cobre, "logiccobreConectionDetector");
 
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) {
    for(G4int i = 0; i < numDetectors; i++) {
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos = -0.305*cm-0.5*(numDetectors-1)*separation + i*separation;
            G4double yPos = -0.5*(numDetectors-1)*separation + j*separation;
            G4double zPos = +PhantomZ/2-profundidadGrupo -layer * (layerSeparation) -SDiodeZ -cobreConectionSDiodeZ/2-PhantomZ/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
            
            G4VPhysicalVolume *physcobreConectionDetector = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobreConectionDetector, "physcobreConectionDetector",
                logicWorld, false, ndethcn  );
            ndethcn++;
        }
}
}

G4RotationMatrix* rotation = new G4RotationMatrix();
rotation->rotateZ(-45.0 * deg);                                                                                                                                                                                                                           



//& debajo de los pads que conecta a solo algunos pads

                                                                                                                                                                                                                                                                            //! abajito xd PISTAS DEL DIODO IZQ 1
                                                                                                                                                                                                                                                                                        cobreConectionSDiodeX = 0.12*mm;
                                                                                                                                                                                                                                                                                            cobreConectionSDiodeY =1.733*mm;
                                                                                                                                                                                                                                                                                            cobreConectionSDiodeZ = 0.0035*cm;
                                                                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                            // Crear el volumen sólido del detector
                                                                                                                                                                                                                                                                                            cobreConectionsolidDetector = new G4Box("cobreConectionsolidDetector", cobreConectionSDiodeX/2, cobreConectionSDiodeY/2, cobreConectionSDiodeZ/2);

                                                                                                                                                                                                                                                                                            // Crear el volumen lógico del detector
                                                                                                                                                                                                                                                                                            logiccobreConectionDetector = new G4LogicalVolume(cobreConectionsolidDetector,cobre, "logiccobreConectionDetector");
                                                                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                            // Colocar los detectores en la cuadrícula y capas
                                                                                                                                                                                                                                                                                            for(G4int layer = 0; layer < numLayers; layer++) { 
                                                                                                                                                                                                                                                                                                    for(G4int j = 0; j < numDetectors; j++) {
                                                                                                                                                                                                                                                                                                        G4double xPos = -0.11*mm-0.305*cm -0.5*(22)*0.8*cm + 0*0.8*cm;
                                                                                                                                                                                                                                                                                                        G4double yPos =-0.09*mm -1.2*mm/2-cobreConectionSDiodeY/2 -0.5*(numDetectors-1)*separation + j*separation;
                                                                                                                                                                                                                                                                                                        G4double zPos = +PhantomZ/2-profundidadGrupo -layer * (layerSeparation) -SDiodeZ -cobreConectionSDiodeZ/2-PhantomZ/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                        G4VPhysicalVolume *physcobreConectionDetector = new G4PVPlacement(
                                                                                                                                                                                                                                                                                                            0,
                                                                                                                                                                                                                                                                                                            G4ThreeVector(xPos, yPos, zPos),
                                                                                                                                                                                                                                                                                                            logiccobreConectionDetector, "physcobreConectionDetector",
                                                                                                                                                                                                                                                                                                            logicWorld, false, ndethcn  );
                                                                                                                                                                                                                                                                                                        ndethcn++;
                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                }
                                                                                                                                                                                                                                                                                }





                                                                                                                                                                                                                                                                            //regresaractual
                                                                                                                                                                                                                                                                            //! abajito xd PISTAS DEL DIODO IZQ 2
                                                                                                                                                                                                                                                                                        cobreConectionSDiodeX = 0.12*mm;
                                                                                                                                                                                                                                                                                            cobreConectionSDiodeY =1.219*mm;
                                                                                                                                                                                                                                                                                            cobreConectionSDiodeZ = 0.0035*cm;
                                                                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                            // Crear el volumen sólido del detector
                                                                                                                                                                                                                                                                                            cobreConectionsolidDetector = new G4Box("cobreConectionsolidDetector", cobreConectionSDiodeX/2, cobreConectionSDiodeY/2, cobreConectionSDiodeZ/2);

                                                                                                                                                                                                                                                                                            // Crear el volumen lógico del detector
                                                                                                                                                                                                                                                                                            logiccobreConectionDetector = new G4LogicalVolume(cobreConectionsolidDetector,cobre, "logiccobreConectionDetector");
                                                                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                            // Colocar los detectores en la cuadrícula y capas
                                                                                                                                                                                                                                                                                            for(G4int layer = 0; layer < numLayers; layer++) {
                                                                                                                                                                                                                                                                                                for(G4int j = 0; j < numDetectors; j++) { 
                                                                                                                                                                                                                                                                                                        G4double xPos =0.862*mm/2 -0.305*cm  -0.5*(22)*0.8*cm + 0*0.8*cm;
                                                                                                                                                                                                                                                                                                        G4double yPos =-0.862*mm/2-1.733*mm -cobreConectionSDiodeY/2 -0.5*(numDetectors-1)*separation + j*separation;
                                                                                                                                                                                                                                                                                                        G4double zPos = +PhantomZ/2-profundidadGrupo -layer * (layerSeparation) -SDiodeZ -cobreConectionSDiodeZ/2-PhantomZ/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                        G4VPhysicalVolume *physcobreConectionDetector = new G4PVPlacement(
                                                                                                                                                                                                                                                                                                            rotation,
                                                                                                                                                                                                                                                                                                            G4ThreeVector(xPos, yPos, zPos),
                                                                                                                                                                                                                                                                                                            logiccobreConectionDetector, "physcobreConectionDetector",
                                                                                                                                                                                                                                                                                                            logicWorld, false, ndethcn  );
                                                                                                                                                                                                                                                                                                        ndethcn++;
                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                }
                                                                                                                                                                                                                                                                                }



//! abajito xd PISTAS DEL DIODO IZQ 1
                                                                                                                                                                                                                                                                                        cobreConectionSDiodeX = 0.12*mm;
                                                                                                                                                                                                                                                                                            cobreConectionSDiodeY =1.622*mm;
                                                                                                                                                                                                                                                                                            cobreConectionSDiodeZ = 0.0035*cm;
                                                                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                            // Crear el volumen sólido del detector
                                                                                                                                                                                                                                                                                            cobreConectionsolidDetector = new G4Box("cobreConectionsolidDetector", cobreConectionSDiodeX/2, cobreConectionSDiodeY/2, cobreConectionSDiodeZ/2);

                                                                                                                                                                                                                                                                                            // Crear el volumen lógico del detector
                                                                                                                                                                                                                                                                                            logiccobreConectionDetector = new G4LogicalVolume(cobreConectionsolidDetector,cobre, "logiccobreConectionDetector");
                                                                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                            // Colocar los detectores en la cuadrícula y capas
                                                                                                                                                                                                                                                                                            for(G4int layer = 0; layer < numLayers; layer++) { 
                                                                                                                                                                                                                                                                                                    for(G4int j = 0; j < numDetectors; j++) {
                                                                                                                                                                                                                                                                                                        G4double xPos = -0.305*cm -0.5*(22)*0.8*cm + 2*0.8*cm;
                                                                                                                                                                                                                                                                                                        G4double yPos =-1.2*mm/2-cobreConectionSDiodeY/2 -0.5*(numDetectors-1)*separation + j*separation;
                                                                                                                                                                                                                                                                                                        G4double zPos = +PhantomZ/2-profundidadGrupo -layer * (layerSeparation) -SDiodeZ -cobreConectionSDiodeZ/2-PhantomZ/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                        G4VPhysicalVolume *physcobreConectionDetector = new G4PVPlacement(
                                                                                                                                                                                                                                                                                                            0,
                                                                                                                                                                                                                                                                                                            G4ThreeVector(xPos, yPos, zPos),
                                                                                                                                                                                                                                                                                                            logiccobreConectionDetector, "physcobreConectionDetector",
                                                                                                                                                                                                                                                                                                            logicWorld, false, ndethcn  );
                                                                                                                                                                                                                                                                                                        ndethcn++;
                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                }
                                                                                                                                                                                                                                                                                }





                                                                                                                                                                                                                                                                            //regresaractual
                                                                                                                                                                                                                                                                            //! abajito xd PISTAS DEL DIODO IZQ 2
                                                                                                                                                                                                                                                                                        cobreConectionSDiodeX = 0.12*mm;
                                                                                                                                                                                                                                                                                            cobreConectionSDiodeY =1.078*mm;
                                                                                                                                                                                                                                                                                            cobreConectionSDiodeZ = 0.0035*cm;
                                                                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                            // Crear el volumen sólido del detector
                                                                                                                                                                                                                                                                                            cobreConectionsolidDetector = new G4Box("cobreConectionsolidDetector", cobreConectionSDiodeX/2, cobreConectionSDiodeY/2, cobreConectionSDiodeZ/2);

                                                                                                                                                                                                                                                                                            // Crear el volumen lógico del detector
                                                                                                                                                                                                                                                                                            logiccobreConectionDetector = new G4LogicalVolume(cobreConectionsolidDetector,cobre, "logiccobreConectionDetector");
                                                                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                            // Colocar los detectores en la cuadrícula y capas
                                                                                                                                                                                                                                                                                            for(G4int layer = 0; layer < numLayers; layer++) {
                                                                                                                                                                                                                                                                                                for(G4int j = 0; j < numDetectors; j++) { 
                                                                                                                                                                                                                                                                                                        G4double xPos =0.762*mm/2 -0.305*cm  -0.5*(22)*0.8*cm + 2*0.8*cm;
                                                                                                                                                                                                                                                                                                        G4double yPos =-0.762*mm/2-1.733*mm -cobreConectionSDiodeY/2 -0.5*(numDetectors-1)*separation + j*separation;
                                                                                                                                                                                                                                                                                                        G4double zPos = +PhantomZ/2-profundidadGrupo -layer * (layerSeparation) -SDiodeZ -cobreConectionSDiodeZ/2-PhantomZ/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                        G4VPhysicalVolume *physcobreConectionDetector = new G4PVPlacement(
                                                                                                                                                                                                                                                                                                            rotation,
                                                                                                                                                                                                                                                                                                            G4ThreeVector(xPos, yPos, zPos),
                                                                                                                                                                                                                                                                                                            logiccobreConectionDetector, "physcobreConectionDetector",
                                                                                                                                                                                                                                                                                                            logicWorld, false, ndethcn  );
                                                                                                                                                                                                                                                                                                        ndethcn++;
                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                }
                                                                                                                                                                                                                                                                                }





//! abajito xd PISTAS DEL DIODO IZQ 1
                                                                                                                                                                                                                                                                                        cobreConectionSDiodeX = 0.12*mm;
                                                                                                                                                                                                                                                                                            cobreConectionSDiodeY =1.346*mm;
                                                                                                                                                                                                                                                                                            cobreConectionSDiodeZ = 0.0035*cm;
                                                                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                            // Crear el volumen sólido del detector
                                                                                                                                                                                                                                                                                            cobreConectionsolidDetector = new G4Box("cobreConectionsolidDetector", cobreConectionSDiodeX/2, cobreConectionSDiodeY/2, cobreConectionSDiodeZ/2);

                                                                                                                                                                                                                                                                                            // Crear el volumen lógico del detector
                                                                                                                                                                                                                                                                                            logiccobreConectionDetector = new G4LogicalVolume(cobreConectionsolidDetector,cobre, "logiccobreConectionDetector");
                                                                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                            // Colocar los detectores en la cuadrícula y capas
                                                                                                                                                                                                                                                                                            for(G4int layer = 0; layer < numLayers; layer++) { 
                                                                                                                                                                                                                                                                                                    for(G4int j = 0; j < numDetectors; j++) {
                                                                                                                                                                                                                                                                                                        G4double xPos = -0.305*cm -0.5*(22)*0.8*cm + 4*0.8*cm;
                                                                                                                                                                                                                                                                                                        G4double yPos =-1.2*mm/2-cobreConectionSDiodeY/2 -0.5*(numDetectors-1)*separation + j*separation;
                                                                                                                                                                                                                                                                                                        G4double zPos = +PhantomZ/2-profundidadGrupo -layer * (layerSeparation) -SDiodeZ -cobreConectionSDiodeZ/2-PhantomZ/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                        G4VPhysicalVolume *physcobreConectionDetector = new G4PVPlacement(
                                                                                                                                                                                                                                                                                                            0,
                                                                                                                                                                                                                                                                                                            G4ThreeVector(xPos, yPos, zPos),
                                                                                                                                                                                                                                                                                                            logiccobreConectionDetector, "physcobreConectionDetector",
                                                                                                                                                                                                                                                                                                            logicWorld, false, ndethcn  );
                                                                                                                                                                                                                                                                                                        ndethcn++;
                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                }
                                                                                                                                                                                                                                                                                }





                                                                                                                                                                                                                                                                            //regresaractual
                                                                                                                                                                                                                                                                            //! abajito xd PISTAS DEL DIODO IZQ 2
                                                                                                                                                                                                                                                                                        cobreConectionSDiodeX = 0.12*mm;
                                                                                                                                                                                                                                                                                            cobreConectionSDiodeY =1.078*mm;
                                                                                                                                                                                                                                                                                            cobreConectionSDiodeZ = 0.0035*cm;
                                                                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                            // Crear el volumen sólido del detector
                                                                                                                                                                                                                                                                                            cobreConectionsolidDetector = new G4Box("cobreConectionsolidDetector", cobreConectionSDiodeX/2, cobreConectionSDiodeY/2, cobreConectionSDiodeZ/2);

                                                                                                                                                                                                                                                                                            // Crear el volumen lógico del detector
                                                                                                                                                                                                                                                                                            logiccobreConectionDetector = new G4LogicalVolume(cobreConectionsolidDetector,cobre, "logiccobreConectionDetector");
                                                                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                            // Colocar los detectores en la cuadrícula y capas
                                                                                                                                                                                                                                                                                            for(G4int layer = 0; layer < numLayers; layer++) {
                                                                                                                                                                                                                                                                                                for(G4int j = 0; j < numDetectors; j++) { 
                                                                                                                                                                                                                                                                                                        G4double xPos =0.762*mm/2 -0.305*cm  -0.5*(22)*0.8*cm + 4*0.8*cm;
                                                                                                                                                                                                                                                                                                        G4double yPos =-0.11*mm-0.762*mm/2-1.346*mm -cobreConectionSDiodeY/2 -0.5*(numDetectors-1)*separation + j*separation;
                                                                                                                                                                                                                                                                                                        G4double zPos = +PhantomZ/2-profundidadGrupo -layer * (layerSeparation) -SDiodeZ -cobreConectionSDiodeZ/2-PhantomZ/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                        G4VPhysicalVolume *physcobreConectionDetector = new G4PVPlacement(
                                                                                                                                                                                                                                                                                                            rotation,
                                                                                                                                                                                                                                                                                                            G4ThreeVector(xPos, yPos, zPos),
                                                                                                                                                                                                                                                                                                            logiccobreConectionDetector, "physcobreConectionDetector",
                                                                                                                                                                                                                                                                                                            logicWorld, false, ndethcn  );
                                                                                                                                                                                                                                                                                                        ndethcn++;
                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                }
                                                                                                                                                                                                                                                                                }











//! abajito xd PISTAS DEL DIODO IZQ 1
                                                                                                                                                                                                                                                                                        cobreConectionSDiodeX = 0.12*mm;
                                                                                                                                                                                                                                                                                            cobreConectionSDiodeY =0.964*mm;
                                                                                                                                                                                                                                                                                            cobreConectionSDiodeZ = 0.0035*cm;
                                                                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                            // Crear el volumen sólido del detector
                                                                                                                                                                                                                                                                                            cobreConectionsolidDetector = new G4Box("cobreConectionsolidDetector", cobreConectionSDiodeX/2, cobreConectionSDiodeY/2, cobreConectionSDiodeZ/2);

                                                                                                                                                                                                                                                                                            // Crear el volumen lógico del detector
                                                                                                                                                                                                                                                                                            logiccobreConectionDetector = new G4LogicalVolume(cobreConectionsolidDetector,cobre, "logiccobreConectionDetector");
                                                                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                            // Colocar los detectores en la cuadrícula y capas
                                                                                                                                                                                                                                                                                            for(G4int layer = 0; layer < numLayers; layer++) { 
                                                                                                                                                                                                                                                                                                    for(G4int j = 0; j < numDetectors; j++) {
                                                                                                                                                                                                                                                                                                        G4double xPos = -0.305*cm -0.5*(22)*0.8*cm + 6*0.8*cm;
                                                                                                                                                                                                                                                                                                        G4double yPos =-1.2*mm/2-cobreConectionSDiodeY/2 -0.5*(numDetectors-1)*separation + j*separation;
                                                                                                                                                                                                                                                                                                        G4double zPos = +PhantomZ/2-profundidadGrupo -layer * (layerSeparation) -SDiodeZ -cobreConectionSDiodeZ/2-PhantomZ/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                        G4VPhysicalVolume *physcobreConectionDetector = new G4PVPlacement(
                                                                                                                                                                                                                                                                                                            0,
                                                                                                                                                                                                                                                                                                            G4ThreeVector(xPos, yPos, zPos),
                                                                                                                                                                                                                                                                                                            logiccobreConectionDetector, "physcobreConectionDetector",
                                                                                                                                                                                                                                                                                                            logicWorld, false, ndethcn  );
                                                                                                                                                                                                                                                                                                        ndethcn++;
                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                }
                                                                                                                                                                                                                                                                                }





                                                                                                                                                                                                                                                                            //regresaractual
                                                                                                                                                                                                                                                                            //! abajito xd PISTAS DEL DIODO IZQ 2
                                                                                                                                                                                                                                                                                        cobreConectionSDiodeX = 0.12*mm;
                                                                                                                                                                                                                                                                                            cobreConectionSDiodeY =1.257*mm;
                                                                                                                                                                                                                                                                                            cobreConectionSDiodeZ = 0.0035*cm;
                                                                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                            // Crear el volumen sólido del detector
                                                                                                                                                                                                                                                                                            cobreConectionsolidDetector = new G4Box("cobreConectionsolidDetector", cobreConectionSDiodeX/2, cobreConectionSDiodeY/2, cobreConectionSDiodeZ/2);

                                                                                                                                                                                                                                                                                            // Crear el volumen lógico del detector
                                                                                                                                                                                                                                                                                            logiccobreConectionDetector = new G4LogicalVolume(cobreConectionsolidDetector,cobre, "logiccobreConectionDetector");
                                                                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                            // Colocar los detectores en la cuadrícula y capas
                                                                                                                                                                                                                                                                                            for(G4int layer = 0; layer < numLayers; layer++) {
                                                                                                                                                                                                                                                                                                for(G4int j = 0; j < numDetectors; j++) { 
                                                                                                                                                                                                                                                                                                        G4double xPos =0.889*mm/2 -0.305*cm  -0.5*(22)*0.8*cm + 6*0.8*cm;
                                                                                                                                                                                                                                                                                                        G4double yPos =-0.07*mm-0.889*mm/2-0.964*mm -cobreConectionSDiodeY/2 -0.5*(numDetectors-1)*separation + j*separation;
                                                                                                                                                                                                                                                                                                        G4double zPos = +PhantomZ/2-profundidadGrupo -layer * (layerSeparation) -SDiodeZ -cobreConectionSDiodeZ/2-PhantomZ/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                        G4VPhysicalVolume *physcobreConectionDetector = new G4PVPlacement(
                                                                                                                                                                                                                                                                                                            rotation,
                                                                                                                                                                                                                                                                                                            G4ThreeVector(xPos, yPos, zPos),
                                                                                                                                                                                                                                                                                                            logiccobreConectionDetector, "physcobreConectionDetector",
                                                                                                                                                                                                                                                                                                            logicWorld, false, ndethcn  );
                                                                                                                                                                                                                                                                                                        ndethcn++;
                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                }
                                                                                                                                                                                                                                                                                }








//! abajito xd PISTAS DEL DIODO IZQ 1
                                                                                                                                                                                                                                                                                        cobreConectionSDiodeX = 0.12*mm;
                                                                                                                                                                                                                                                                                            cobreConectionSDiodeY =0.835*mm;
                                                                                                                                                                                                                                                                                            cobreConectionSDiodeZ = 0.0035*cm;
                                                                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                            // Crear el volumen sólido del detector
                                                                                                                                                                                                                                                                                            cobreConectionsolidDetector = new G4Box("cobreConectionsolidDetector", cobreConectionSDiodeX/2, cobreConectionSDiodeY/2, cobreConectionSDiodeZ/2);

                                                                                                                                                                                                                                                                                            // Crear el volumen lógico del detector
                                                                                                                                                                                                                                                                                            logiccobreConectionDetector = new G4LogicalVolume(cobreConectionsolidDetector,cobre, "logiccobreConectionDetector");
                                                                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                            // Colocar los detectores en la cuadrícula y capas
                                                                                                                                                                                                                                                                                            for(G4int layer = 0; layer < numLayers; layer++) { 
                                                                                                                                                                                                                                                                                                    for(G4int j = 0; j < numDetectors; j++) {
                                                                                                                                                                                                                                                                                                        G4double xPos = -0.305*cm -0.5*(22)*0.8*cm + 8*0.8*cm;
                                                                                                                                                                                                                                                                                                        G4double yPos =-1.2*mm/2-cobreConectionSDiodeY/2 -0.5*(numDetectors-1)*separation + j*separation;
                                                                                                                                                                                                                                                                                                        G4double zPos = +PhantomZ/2-profundidadGrupo -layer * (layerSeparation) -SDiodeZ -cobreConectionSDiodeZ/2-PhantomZ/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                        G4VPhysicalVolume *physcobreConectionDetector = new G4PVPlacement(
                                                                                                                                                                                                                                                                                                            0,
                                                                                                                                                                                                                                                                                                            G4ThreeVector(xPos, yPos, zPos),
                                                                                                                                                                                                                                                                                                            logiccobreConectionDetector, "physcobreConectionDetector",
                                                                                                                                                                                                                                                                                                            logicWorld, false, ndethcn  );
                                                                                                                                                                                                                                                                                                        ndethcn++;
                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                }
                                                                                                                                                                                                                                                                                }





                                                                                                                                                                                                                                                                            //regresaractual
                                                                                                                                                                                                                                                                            //! abajito xd PISTAS DEL DIODO IZQ 2
                                                                                                                                                                                                                                                                                        cobreConectionSDiodeX = 0.12*mm;
                                                                                                                                                                                                                                                                                            cobreConectionSDiodeY =1.068*mm;
                                                                                                                                                                                                                                                                                            cobreConectionSDiodeZ = 0.0035*cm;
                                                                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                            // Crear el volumen sólido del detector
                                                                                                                                                                                                                                                                                            cobreConectionsolidDetector = new G4Box("cobreConectionsolidDetector", cobreConectionSDiodeX/2, cobreConectionSDiodeY/2, cobreConectionSDiodeZ/2);

                                                                                                                                                                                                                                                                                            // Crear el volumen lógico del detector
                                                                                                                                                                                                                                                                                            logiccobreConectionDetector = new G4LogicalVolume(cobreConectionsolidDetector,cobre, "logiccobreConectionDetector");
                                                                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                            // Colocar los detectores en la cuadrícula y capas
                                                                                                                                                                                                                                                                                            for(G4int layer = 0; layer < numLayers; layer++) {
                                                                                                                                                                                                                                                                                                for(G4int j = 0; j < numDetectors; j++) { 
                                                                                                                                                                                                                                                                                                        G4double xPos =0.762*mm/2 -0.305*cm  -0.5*(22)*0.8*cm + 8*0.8*cm;
                                                                                                                                                                                                                                                                                                        G4double yPos =-0.11*mm-0.762*mm/2-0.835*mm -cobreConectionSDiodeY/2 -0.5*(numDetectors-1)*separation + j*separation;
                                                                                                                                                                                                                                                                                                        G4double zPos = +PhantomZ/2-profundidadGrupo -layer * (layerSeparation) -SDiodeZ -cobreConectionSDiodeZ/2-PhantomZ/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                        G4VPhysicalVolume *physcobreConectionDetector = new G4PVPlacement(
                                                                                                                                                                                                                                                                                                            rotation,
                                                                                                                                                                                                                                                                                                            G4ThreeVector(xPos, yPos, zPos),
                                                                                                                                                                                                                                                                                                            logiccobreConectionDetector, "physcobreConectionDetector",
                                                                                                                                                                                                                                                                                                            logicWorld, false, ndethcn  );
                                                                                                                                                                                                                                                                                                        ndethcn++;
                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                }
                                                                                                                                                                                                                                                                                }












//! abajito xd PISTAS DEL DIODO IZQ 1
                                                                                                                                                                                                                                                                                        cobreConectionSDiodeX = 0.12*mm;
                                                                                                                                                                                                                                                                                            cobreConectionSDiodeY =0.722*mm;
                                                                                                                                                                                                                                                                                            cobreConectionSDiodeZ = 0.0035*cm;
                                                                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                            // Crear el volumen sólido del detector
                                                                                                                                                                                                                                                                                            cobreConectionsolidDetector = new G4Box("cobreConectionsolidDetector", cobreConectionSDiodeX/2, cobreConectionSDiodeY/2, cobreConectionSDiodeZ/2);

                                                                                                                                                                                                                                                                                            // Crear el volumen lógico del detector
                                                                                                                                                                                                                                                                                            logiccobreConectionDetector = new G4LogicalVolume(cobreConectionsolidDetector,cobre, "logiccobreConectionDetector");
                                                                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                            // Colocar los detectores en la cuadrícula y capas
                                                                                                                                                                                                                                                                                            for(G4int layer = 0; layer < numLayers; layer++) { 
                                                                                                                                                                                                                                                                                                    for(G4int j = 0; j < numDetectors; j++) {
                                                                                                                                                                                                                                                                                                        G4double xPos = -0.305*cm -0.5*(22)*0.8*cm + 10*0.8*cm;
                                                                                                                                                                                                                                                                                                        G4double yPos =-1.2*mm/2-cobreConectionSDiodeY/2 -0.5*(numDetectors-1)*separation + j*separation;
                                                                                                                                                                                                                                                                                                        G4double zPos = +PhantomZ/2-profundidadGrupo -layer * (layerSeparation) -SDiodeZ -cobreConectionSDiodeZ/2-PhantomZ/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                        G4VPhysicalVolume *physcobreConectionDetector = new G4PVPlacement(
                                                                                                                                                                                                                                                                                                            0,
                                                                                                                                                                                                                                                                                                            G4ThreeVector(xPos, yPos, zPos),
                                                                                                                                                                                                                                                                                                            logiccobreConectionDetector, "physcobreConectionDetector",
                                                                                                                                                                                                                                                                                                            logicWorld, false, ndethcn  );
                                                                                                                                                                                                                                                                                                        ndethcn++;
                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                }
                                                                                                                                                                                                                                                                                }





                                                                                                                                                                                                                                                                            //regresaractual
                                                                                                                                                                                                                                                                            //! abajito xd PISTAS DEL DIODO IZQ 2
                                                                                                                                                                                                                                                                                        cobreConectionSDiodeX = 0.12*mm;
                                                                                                                                                                                                                                                                                            cobreConectionSDiodeY =0.898*mm;
                                                                                                                                                                                                                                                                                            cobreConectionSDiodeZ = 0.0035*cm;
                                                                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                            // Crear el volumen sólido del detector
                                                                                                                                                                                                                                                                                            cobreConectionsolidDetector = new G4Box("cobreConectionsolidDetector", cobreConectionSDiodeX/2, cobreConectionSDiodeY/2, cobreConectionSDiodeZ/2);

                                                                                                                                                                                                                                                                                            // Crear el volumen lógico del detector
                                                                                                                                                                                                                                                                                            logiccobreConectionDetector = new G4LogicalVolume(cobreConectionsolidDetector,cobre, "logiccobreConectionDetector");
                                                                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                            // Colocar los detectores en la cuadrícula y capas
                                                                                                                                                                                                                                                                                            for(G4int layer = 0; layer < numLayers; layer++) {
                                                                                                                                                                                                                                                                                                for(G4int j = 0; j < numDetectors; j++) { 
                                                                                                                                                                                                                                                                                                        G4double xPos =0.635*mm/2 -0.305*cm  -0.5*(22)*0.8*cm + 10*0.8*cm;
                                                                                                                                                                                                                                                                                                        G4double yPos =-0.16*mm-0.722*mm/2-0.722*mm -cobreConectionSDiodeY/2 -0.5*(numDetectors-1)*separation + j*separation;
                                                                                                                                                                                                                                                                                                        G4double zPos = +PhantomZ/2-profundidadGrupo -layer * (layerSeparation) -SDiodeZ -cobreConectionSDiodeZ/2-PhantomZ/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                        G4VPhysicalVolume *physcobreConectionDetector = new G4PVPlacement(
                                                                                                                                                                                                                                                                                                            rotation,
                                                                                                                                                                                                                                                                                                            G4ThreeVector(xPos, yPos, zPos),
                                                                                                                                                                                                                                                                                                            logiccobreConectionDetector, "physcobreConectionDetector",
                                                                                                                                                                                                                                                                                                            logicWorld, false, ndethcn  );
                                                                                                                                                                                                                                                                                                        ndethcn++;
                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                }
                                                                                                                                                                                                                                                                                }




//& debajo de los pads que conecta a solo algunos pads






















    //!PISTAS DEL DIODO IZQ 1
            cobreConectionSDiodeX = 0.15*mm;
                cobreConectionSDiodeY = 2.555*mm;
                cobreConectionSDiodeZ = 0.0035*cm;
                
                // Crear el volumen sólido del detector
                cobreConectionsolidDetector = new G4Box("cobreConectionsolidDetector", cobreConectionSDiodeX/2, cobreConectionSDiodeY/2, cobreConectionSDiodeZ/2);

                // Crear el volumen lógico del detector
                logiccobreConectionDetector = new G4LogicalVolume(cobreConectionsolidDetector,cobre, "logiccobreConectionDetector");
                
                // Colocar los detectores en la cuadrícula y capas
                for(G4int layer = 0; layer < numLayers; layer++) {
                    for(G4int i = 0; i < numDetectors; i++) {
                        for(G4int j = 0; j < numDetectors; j++) {
                            G4double xPos = -0.305*cm-0.5*(numDetectors-1)*separation + i*separation;
                            G4double yPos =0.12*cm/2+ cobreConectionSDiodeY/2-0.5*(numDetectors-1)*separation + j*separation;
                            G4double zPos = +PhantomZ/2-profundidadGrupo -layer * (layerSeparation) -SDiodeZ -cobreConectionSDiodeZ/2-PhantomZ/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                            
                            G4VPhysicalVolume *physcobreConectionDetector = new G4PVPlacement(
                                0,
                                G4ThreeVector(xPos, yPos, zPos),
                                logiccobreConectionDetector, "physcobreConectionDetector",
                                logicWorld, false, ndethcn  );
                            ndethcn++;
                        }
    }
    }






G4RotationMatrix* rotation2 = new G4RotationMatrix();
rotation2->rotateZ( 45.0 * deg);
G4RotationMatrix* rotation3 = new G4RotationMatrix();
rotation3->rotateZ( 90.0 * deg);
                        //!PISTAS DEL DIODO IZQ 2
                                    cobreConectionSDiodeX = 0.15*mm;
                                        cobreConectionSDiodeY = 1.195*mm;
                                        cobreConectionSDiodeZ = 0.0035*cm;
                                        
                                        // Crear el volumen sólido del detector
                                        cobreConectionsolidDetector = new G4Box("cobreConectionsolidDetector", cobreConectionSDiodeX/2, cobreConectionSDiodeY/2, cobreConectionSDiodeZ/2);

                                        // Crear el volumen lógico del detector
                                        logiccobreConectionDetector = new G4LogicalVolume(cobreConectionsolidDetector,cobre, "logiccobreConectionDetector");
                                        
                                        // Colocar los detectores en la cuadrícula y capas
                                        for(G4int layer = 0; layer < numLayers; layer++) {
                                            for(G4int i = 0; i < numDetectors; i++) {
                                                for(G4int j = 0; j < numDetectors; j++) {
                                                    G4double xPos = + 0.845*mm/2 -0.305*cm-0.5*(numDetectors-1)*separation + i*separation;
                                                    G4double yPos =0.12*cm/2+ 2.555*mm + 0.845*mm/2-0.5*(numDetectors-1)*separation + j*separation;
                                                    G4double zPos = +PhantomZ/2-profundidadGrupo -layer * (layerSeparation) -SDiodeZ -cobreConectionSDiodeZ/2-PhantomZ/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                                    
                                                    G4VPhysicalVolume *physcobreConectionDetector = new G4PVPlacement(
                                                        rotation2,
                                                        G4ThreeVector(xPos, yPos, zPos),
                                                        logiccobreConectionDetector, "physcobreConectionDetector",
                                                        logicWorld, false, ndethcn  );
                                                    ndethcn++;
                                                }
                            }
                            }

                                                        //!PISTAS DEL DIODO IZQ 3
                                                                    cobreConectionSDiodeX = 0.15*mm;
                                                                        cobreConectionSDiodeY = 1.105*mm;
                                                                        cobreConectionSDiodeZ = 0.0035*cm;
                                                                        
                                                                        // Crear el volumen sólido del detector
                                                                        cobreConectionsolidDetector = new G4Box("cobreConectionsolidDetector", cobreConectionSDiodeX/2, cobreConectionSDiodeY/2, cobreConectionSDiodeZ/2);

                                                                        // Crear el volumen lógico del detector
                                                                        logiccobreConectionDetector = new G4LogicalVolume(cobreConectionsolidDetector,cobre, "logiccobreConectionDetector");
                                                                        
                                                                        // Colocar los detectores en la cuadrícula y capas
                                                                        for(G4int layer = 0; layer < numLayers; layer++) {
                                                                            for(G4int i = 0; i < numDetectors; i++) {
                                                                                for(G4int j = 0; j < numDetectors; j++) {//al girar el Y se convierte en x por eso se suma el Ytt
                                                                                    G4double xPos =cobreConectionSDiodeY/2 + 0.845*mm  -0.305*cm-0.5*(numDetectors-1)*separation + i*separation;
                                                                                    G4double yPos =+0.12*cm/2+ 2.555*mm + 0.845*mm -0.5*(numDetectors-1)*separation + j*separation;
                                                                                    G4double zPos = +PhantomZ/2-profundidadGrupo -layer * (layerSeparation) -SDiodeZ -cobreConectionSDiodeZ/2-PhantomZ/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                                                                    
                                                                                    G4VPhysicalVolume *physcobreConectionDetector = new G4PVPlacement(
                                                                                        rotation3,
                                                                                        G4ThreeVector(xPos, yPos, zPos),
                                                                                        logiccobreConectionDetector, "physcobreConectionDetector",
                                                                                        logicWorld, false, ndethcn  );
                                                                                    ndethcn++;
                                                                                }
                                                            }
                                                            }






//! +0.305*cm+0.12*cm/2 + 5*cm

      //!PISTAS DEL DIODO DER  

cobreConectionSDiodeX = 0.12*cm;
cobreConectionSDiodeY = 0.18*cm;
cobreConectionSDiodeZ = 0.007*cm;

// Crear el volumen sólido del detector
    cobreConectionsolidDetector = new G4Box("cobreConectionsolidDetector", cobreConectionSDiodeX/2, cobreConectionSDiodeY/2, cobreConectionSDiodeZ/2);

// Crear el volumen lógico del detector
   logiccobreConectionDetector = new G4LogicalVolume(cobreConectionsolidDetector,cobre, "logiccobreConectionDetector");
 
for(G4int layer = 0; layer < numLayers; layer++) {
    for(G4int i = 0; i < numDetectors; i++) {
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos = +0.305*cm-0.5*(numDetectors-1)*separation + i*separation;
            G4double yPos = -0.5*(numDetectors-1)*separation + j*separation;
            G4double zPos = +PhantomZ/2-profundidadGrupo -layer * (layerSeparation) -SDiodeZ -cobreConectionSDiodeZ/2-PhantomZ/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
            
            G4VPhysicalVolume *physcobreConectionDetector = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobreConectionDetector, "physcobreConectionDetector",
                logicWorld, false, ndethcn   );
            ndethcn++;
        }
}
}
   


            //!PISTAS DEL DIODO DER 1
                    cobreConectionSDiodeX = 0.2*mm;
                        cobreConectionSDiodeY = 2.245*mm;
                        cobreConectionSDiodeZ = 0.0035*cm;
                        
                        // Crear el volumen sólido del detector
                        cobreConectionsolidDetector = new G4Box("cobreConectionsolidDetector", cobreConectionSDiodeX/2, cobreConectionSDiodeY/2, cobreConectionSDiodeZ/2);

                        // Crear el volumen lógico del detector
                        logiccobreConectionDetector = new G4LogicalVolume(cobreConectionsolidDetector,cobre, "logiccobreConectionDetector");
                        
                        // Colocar los detectores en la cuadrícula y capas
                        for(G4int layer = 0; layer < numLayers; layer++) {
                            for(G4int i = 0; i < numDetectors; i++) {
                                for(G4int j = 0; j < numDetectors; j++) {
                                    G4double xPos = +0.305*cm-0.5*(numDetectors-1)*separation + i*separation;
                                    G4double yPos =0.18*cm/2+ cobreConectionSDiodeY/2-0.5*(numDetectors-1)*separation + j*separation;
                                    G4double zPos = +PhantomZ/2-profundidadGrupo -layer * (layerSeparation) -SDiodeZ -cobreConectionSDiodeZ/2-PhantomZ/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                    
                                    G4VPhysicalVolume *physcobreConectionDetector = new G4PVPlacement(
                                        0,
                                        G4ThreeVector(xPos, yPos, zPos),
                                        logiccobreConectionDetector, "physcobreConectionDetector",
                                        logicWorld, false, ndethcn  );
                                    ndethcn++;
                                }
            }
            }



                        //!PISTAS DEL DIODO DER 1
                                cobreConectionSDiodeX = 0.2*mm;
                                    cobreConectionSDiodeY = 1.209*mm;
                                    cobreConectionSDiodeZ = 0.0035*cm;
                                    
                                    // Crear el volumen sólido del detector
                                    cobreConectionsolidDetector = new G4Box("cobreConectionsolidDetector", cobreConectionSDiodeX/2, cobreConectionSDiodeY/2, cobreConectionSDiodeZ/2);

                                    // Crear el volumen lógico del detector
                                    logiccobreConectionDetector = new G4LogicalVolume(cobreConectionsolidDetector,cobre, "logiccobreConectionDetector");
                                    
                                    // Colocar los detectores en la cuadrícula y capas
                                    for(G4int layer = 0; layer < numLayers; layer++) {
                                        for(G4int i = 0; i < numDetectors; i++) {
                                            for(G4int j = 0; j < numDetectors; j++) {
                                                G4double xPos = - 0.845*mm/2 +0.305*cm-0.5*(numDetectors-1)*separation + i*separation;
                                                G4double yPos = 2.245*mm+ 0.845*mm/2+ 0.18*cm/2 -0.5*(numDetectors-1)*separation + j*separation;
                                                G4double zPos = +PhantomZ/2-profundidadGrupo -layer * (layerSeparation) -SDiodeZ -cobreConectionSDiodeZ/2-PhantomZ/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                                
                                                G4VPhysicalVolume *physcobreConectionDetector = new G4PVPlacement(
                                                    rotation,
                                                    G4ThreeVector(xPos, yPos, zPos),
                                                    logiccobreConectionDetector, "physcobreConectionDetector",
                                                    logicWorld, false, ndethcn  );
                                                ndethcn++;
                                            }
                        }
                        }



                                                //!PISTAS DEL DIODO DER 1
                                                        cobreConectionSDiodeX = 0.2*mm;
                                                            cobreConectionSDiodeY = 1.095*mm;
                                                            cobreConectionSDiodeZ = 0.0035*cm;
                                                            
                                                            // Crear el volumen sólido del detector
                                                            cobreConectionsolidDetector = new G4Box("cobreConectionsolidDetector", cobreConectionSDiodeX/2, cobreConectionSDiodeY/2, cobreConectionSDiodeZ/2);

                                                            // Crear el volumen lógico del detector
                                                            logiccobreConectionDetector = new G4LogicalVolume(cobreConectionsolidDetector,cobre, "logiccobreConectionDetector");
                                                            
                                                            // Colocar los detectores en la cuadrícula y capas
                                                            for(G4int layer = 0; layer < numLayers; layer++) {//al girar el Y se convierte en x por eso se suma el Ytt
                                                                for(G4int i = 0; i < numDetectors; i++) {
                                                                    for(G4int j = 0; j < numDetectors; j++) {
                                                                        G4double xPos =-1*cobreConectionSDiodeY/2 - 0.845*mm  +0.305*cm-0.5*(numDetectors-1)*separation + i*separation;
                                                                        G4double yPos = 2.245*mm+ 0.845*mm + 0.18*cm/2 -0.5*(numDetectors-1)*separation + j*separation;
                                                                        G4double zPos = +PhantomZ/2-profundidadGrupo -layer * (layerSeparation) -SDiodeZ -cobreConectionSDiodeZ/2-PhantomZ/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                                                        
                                                                        G4VPhysicalVolume *physcobreConectionDetector = new G4PVPlacement(
                                                                            rotation3,
                                                                            G4ThreeVector(xPos, yPos, zPos),
                                                                            logiccobreConectionDetector, "physcobreConectionDetector",
                                                                            logicWorld, false, ndethcn  );
                                                                        ndethcn++;
                                                                    }
                                                }
                                                }





//? Esta va solo en la parte derecha



                                                //!PISTAS DEL DIODO DER 1
                                                        cobreConectionSDiodeX = 52.234*mm;
                                                            cobreConectionSDiodeY = 0.2*mm;
                                                            cobreConectionSDiodeZ = 0.0035*cm;
                                                            
                                                            // Crear el volumen sólido del detector
                                                            cobreConectionsolidDetector = new G4Box("cobreConectionsolidDetector", cobreConectionSDiodeX/2, cobreConectionSDiodeY/2, cobreConectionSDiodeZ/2);

                                                            // Crear el volumen lógico del detector
                                                            logiccobreConectionDetector = new G4LogicalVolume(cobreConectionsolidDetector,cobre, "logiccobreConectionDetector");
                                                            
                                                            // Colocar los detectores en la cuadrícula y capas
                                                            for(G4int layer = 0; layer < numLayers; layer++) {//al girar el Y se convierte en x por eso se suma el Ytt
                                                                for(G4int j = 0; j < numDetectors; j++) {
                                                                     
                                                                        G4double xPos =2.343*mm+8.8*cm+cobreConectionSDiodeX/2;
                                                                        G4double yPos = 4*mm+ j*separation -0.5*(numDetectors-1)*separation;
                                                                        G4double zPos = +PhantomZ/2-profundidadGrupo -layer * (layerSeparation) -SDiodeZ -cobreConectionSDiodeZ/2-PhantomZ/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                                                        
                                                                        G4VPhysicalVolume *physcobreConectionDetector = new G4PVPlacement(
                                                                            0,
                                                                            G4ThreeVector(xPos, yPos, zPos),
                                                                            logiccobreConectionDetector, "physcobreConectionDetector",
                                                                            logicWorld, false, ndethcn  );
                                                                        ndethcn++;
                                                                    
                                                }
                                                }



//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^v
//^^^^PATITAS    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^v
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^v


//*-PATA IZQUIERDA   EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEINTERNOS-
 G4double cobrePatitaBaseSDiodeX = 0.55*mm;
 G4double cobrePatitaBaseSDiodeY =0.9*mm;
 G4double cobrePatitaBaseSDiodeZ = 0.1*mm;
  G4int ndetpat =1;
// Crear el volumen sólido del detector
   G4Box *  cobrePatitaBasesolidDetector = new G4Box("solidPatitaBase", cobrePatitaBaseSDiodeX/2, cobrePatitaBaseSDiodeY/2, cobrePatitaBaseSDiodeZ/2);

// Crear el volumen lógico del detector
   G4LogicalVolume * logiccobrePatitaBaseDetector = new G4LogicalVolume(cobrePatitaBasesolidDetector,cobre, "logicPatitaBaseDetector");
 
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) {
    for(G4int i = 0; i < numDetectors; i++) {
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos = -0.55*mm -0.450*cm/2-cobrePatitaBaseSDiodeX /2 -0.5*(numDetectors-1)*separation + i*separation;
            G4double yPos = -0.5*(numDetectors-1)*separation + j*separation;
            G4double zPos = +PhantomZ/2-profundidadGrupo -layer * (layerSeparation) -SDiodeZ +cobrePatitaBaseSDiodeZ/2-PhantomZ/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
            
            G4VPhysicalVolume *physcobrePatitaBaseDetector = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobrePatitaBaseDetector, "physcobrePatitaBaseDetector",
                logicWorld, false,   ndetpat  );
             ndetpat++;
        }
}
}

  cobrePatitaBaseSDiodeX = 0.55*mm;
 cobrePatitaBaseSDiodeY = 1.7*mm;
 cobrePatitaBaseSDiodeZ = 0.1*mm;

   cobrePatitaBasesolidDetector = new G4Box("solidPatitaBase", cobrePatitaBaseSDiodeX/2, cobrePatitaBaseSDiodeY/2, cobrePatitaBaseSDiodeZ/2);

// Crear el volumen lógico del detector
  logiccobrePatitaBaseDetector = new G4LogicalVolume(cobrePatitaBasesolidDetector,cobre, "logicPatitaBaseDetector");
 
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) {
    for(G4int i = 0; i < numDetectors; i++) {
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =  -0.450*cm/2-cobrePatitaBaseSDiodeX /2 -0.5*(numDetectors-1)*separation + i*separation;
            G4double yPos = -0.5*(numDetectors-1)*separation + j*separation;
            G4double zPos = +PhantomZ/2-profundidadGrupo -layer * (layerSeparation)  +cobrePatitaBaseSDiodeZ/2-PhantomZ/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
            
            G4VPhysicalVolume *physcobrePatitaBaseDetector = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobrePatitaBaseDetector, "physcobrePatitaBaseDetector",
                logicWorld, false,  ndetpat  );
             ndetpat++;
        }
}
}  


 G4double cobrePatitaBaseSDiodeX2 = 0.1*mm;
 G4double cobrePatitaBaseSDiodeY2 = 1.7*mm;
 G4double cobrePatitaBaseSDiodeZ2 = 0.5*mm;
 
// Crear el volumen sólido del detector
   G4Box *  cobrePatitaBasesolidDetector2 = new G4Box("solidPatitaBase2", cobrePatitaBaseSDiodeX2/2, cobrePatitaBaseSDiodeY2/2, cobrePatitaBaseSDiodeZ2/2);

// Crear el volumen lógico del detector
   G4LogicalVolume * logiccobrePatitaBaseDetector2 = new G4LogicalVolume(cobrePatitaBasesolidDetector2,cobre, "logicPatitaBaseDetector2");
 
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) {
    for(G4int i = 0; i < numDetectors; i++) {
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =  -cobrePatitaBaseSDiodeX  -0.450*cm/2-cobrePatitaBaseSDiodeX2 /2 -0.5*(numDetectors-1)*separation + i*separation;
            G4double yPos = -0.5*(numDetectors-1)*separation + j*separation;
            G4double zPos =  +PhantomZ/2-profundidadGrupo -layer * (layerSeparation) -SDiodeZ  +cobrePatitaBaseSDiodeZ2/2+cobrePatitaBaseSDiodeZ-PhantomZ/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
            
            G4VPhysicalVolume *physcobrePatitaBaseDetector2 = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobrePatitaBaseDetector2, "physcobrePatitaBaseDetector2",
                logicWorld, false,  ndetpat );
             ndetpat++;
        }
}
}  



//*-PATA derecha   EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEINTERNOS-
 
  
 
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) {
    for(G4int i = 0; i < numDetectors; i++) {
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos = +0.55*mm +0.450*cm/2+cobrePatitaBaseSDiodeX /2 -0.5*(numDetectors-1)*separation + i*separation;
            G4double yPos = -0.5*(numDetectors-1)*separation + j*separation;
            G4double zPos =+PhantomZ/2-profundidadGrupo -layer * (layerSeparation) -SDiodeZ +cobrePatitaBaseSDiodeZ/2 -PhantomZ/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
            
            G4VPhysicalVolume *physcobrePatitaBaseDetector = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobrePatitaBaseDetector, "physcobrePatitaBaseDetector",
                logicWorld, false,  ndetpat  );
             ndetpat++;
        }
}
}


// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) {
    for(G4int i = 0; i < numDetectors; i++) {
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =  +0.450*cm/2+cobrePatitaBaseSDiodeX /2 -0.5*(numDetectors-1)*separation + i*separation;
            G4double yPos = -0.5*(numDetectors-1)*separation + j*separation;
            G4double zPos =+PhantomZ/2-profundidadGrupo -layer * (layerSeparation)  +cobrePatitaBaseSDiodeZ/2 -PhantomZ/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
            
            G4VPhysicalVolume *physcobrePatitaBaseDetector = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobrePatitaBaseDetector, "physcobrePatitaBaseDetector",
                logicWorld, false,   ndetpat  );
             ndetpat++;
        }
}
}  

 
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) {
    for(G4int i = 0; i < numDetectors; i++) {
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =  +cobrePatitaBaseSDiodeX  +0.450*cm/2+cobrePatitaBaseSDiodeX2 /2 -0.5*(numDetectors-1)*separation + i*separation;
            G4double yPos = -0.5*(numDetectors-1)*separation + j*separation;
            G4double zPos =  +PhantomZ/2-profundidadGrupo -layer * (layerSeparation) -SDiodeZ  +cobrePatitaBaseSDiodeZ2/2+cobrePatitaBaseSDiodeZ-PhantomZ/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
            
            G4VPhysicalVolume *physcobrePatitaBaseDetector2 = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobrePatitaBaseDetector2, "physcobrePatitaBaseDetector2",
                logicWorld, false,  ndetpat  );
            ndetpat++;
        }
}
}  





 
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//* 
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*


//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//* CAPACITORES       //*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*


//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^v
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^v
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^v


//*-INTERNOS-
// Dimensiones del capacitor
G4double CapacitorX = 0.16 *cm;
G4double CapacitorY = 0.08*cm;
G4double CapacitorZ = 0.08*cm;
// Separación entre capacitores 
G4double shiftCapa=0*cm; //-0.13
// Crear el volumen sólido del capacitor
  G4Box * solidcapacitor = new G4Box("solidcapacitor", CapacitorX/2, CapacitorY/2, CapacitorZ/2);

// Crear el volumen lógico del capacitor
  G4LogicalVolume * logiccapacitor = new G4LogicalVolume(solidcapacitor, BaTiO3, "logiccapacitor");
    G4int ncapacitor=0;
// Colocar los capacitores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) {
    for(G4int i = 0; i < numDetectors; i++) {
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPosC = -0.5*(numDetectors-1)*separation + i*separation;
            G4double yPosC =SDiodeY/2  +(separation -SDiodeY)/2  -0.5*(numDetectors-1)*separation + j*separation+shiftCapa;
            G4double zPosC =-PhantomZ/2-0.02*cm +PhantomZ/2-profundidadGrupo -layer * ( layerSeparation)  ;
            
            G4VPhysicalVolume *physcapacitor = new G4PVPlacement(
                0,
                G4ThreeVector(xPosC, yPosC, zPosC),
                logiccapacitor,
                "physcapacitor",
                logicWorld,
                false,
                ncapacitor
                
            );
            ncapacitor ++;
        }
}
}
 



 
 

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^v
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^v
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^v



//*-COBREEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEINTERNOS-
// Dimensiones del capacitor
G4double cobreConectionCapacitorX = 0.08 *cm;
G4double cobreConectionCapacitorY = 0.09*cm;
G4double cobreConectionCapacitorZ = 0.007*cm;
// Separación entre capacitores 

// Crear el volumen sólido del capacitor
  G4Box * solidcobreConectioncapacitor = new G4Box("solidcobreConectioncapacitor", cobreConectionCapacitorX/2, cobreConectionCapacitorY/2, cobreConectionCapacitorZ/2);

// Crear el volumen lógico del capacitor
  G4LogicalVolume * logiccobreConectioncapacitor = new G4LogicalVolume(solidcobreConectioncapacitor, cobre, "logiccobreConectioncapacitor");
   
// Colocar los capacitores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) {
    for(G4int i = 0; i < numDetectors; i++) {
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPosC =-0.065*cm -0.5*(numDetectors-1)*separation + i*separation;
            G4double yPosC =SDiodeY/2  +(separation -SDiodeY)/2  -0.5*(numDetectors-1)*separation + j*separation+shiftCapa;
            G4double zPosC =-PhantomZ/2-0.02*cm +PhantomZ/2-profundidadGrupo -layer * ( layerSeparation) -cobreConectionCapacitorZ /2-CapacitorZ/2 ;
            
            G4VPhysicalVolume *physcobreConectioncapacitor = new G4PVPlacement(
                0,
                G4ThreeVector(xPosC, yPosC, zPosC),
                logiccobreConectioncapacitor,
                "physcobreConectioncapacitor",
                logicWorld,
                false,
                ncapacitor
                
            );
            ncapacitor ++;
        }
}
}

for(G4int layer = 0; layer < numLayers; layer++) {
    for(G4int i = 0; i < numDetectors; i++) {
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPosC =+0.065*cm -0.5*(numDetectors-1)*separation + i*separation;
            G4double yPosC =SDiodeY/2  +(separation -SDiodeY)/2  -0.5*(numDetectors-1)*separation + j*separation+shiftCapa;
            G4double zPosC =-PhantomZ/2-0.02*cm +PhantomZ/2-profundidadGrupo -layer * ( layerSeparation) -cobreConectionCapacitorZ /2-CapacitorZ/2 ;
            
            G4VPhysicalVolume *physcobreConectioncapacitor = new G4PVPlacement(
                0,
                G4ThreeVector(xPosC, yPosC, zPosC),
                logiccobreConectioncapacitor,
                "physcobreConectioncapacitor",
                logicWorld,
                false,
               ncapacitor
                
            );
           ncapacitor ++;
        }
}
}


 

for(G4int layer = 0; layer < numLayers; layer++) {
    for(G4int i = 0; i < outerNumDetectors; i++) {
        for(G4int j = 0; j < outerNumDetectors; j++) {
            // Calcular la posición del detector externo
            G4double xPos =-0.15*cm  -0.5*(outerNumDetectors-1)* outerSeparation+ i* outerSeparation;
            G4double yPos = 0.265*cm/2  +(separation -0.265*cm)/2  -0.5*(outerNumDetectors-1)* outerSeparation + j* outerSeparation+shiftCapa;

            // Asegurarse de no superponer los detectores centrales
            G4double centralHalfWidthx = +0.5 * (numDetectors-2) * separation+outerSeparation ;
            G4double centralHalfWidthy =  0.265*cm/2  +(separation -0.265*cm)/2  +0.5*(numDetectors-1)*separation + separation-0.5*cm ;
             
            if (std::abs(xPos) >= centralHalfWidthx || std::abs(yPos) >=  centralHalfWidthy) {
               
            G4double zPos = -PhantomZ/2-0.02*cm+PhantomZ/2-profundidadGrupo - layer * layerSeparation -cobreConectionCapacitorZ/2 -CapacitorZ/2;
 
            G4VPhysicalVolume *physcapacitor = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobreConectioncapacitor,
               "physcobreConectioncapacitor",
                logicWorld,
                false,
               ncapacitor
                
            );
            ncapacitor++;
            }

        }
    }
}


for(G4int layer = 0; layer < numLayers; layer++) {
    for(G4int i = 0; i < outerNumDetectors; i++) {
        for(G4int j = 0; j < outerNumDetectors; j++) {
            // Calcular la posición del detector externo
            G4double xPos =+0.15*cm  -0.5*(outerNumDetectors-1)* outerSeparation+ i* outerSeparation;
            G4double yPos = 0.265*cm/2  +(separation -0.265*cm)/2  -0.5*(outerNumDetectors-1)* outerSeparation + j* outerSeparation+shiftCapa;

            // Asegurarse de no superponer los detectores centrales
            G4double centralHalfWidthx = +0.5 * (numDetectors-2) * separation+outerSeparation ;
            G4double centralHalfWidthy =  0.265*cm/2  +(separation -0.265*cm)/2  +0.5*(numDetectors-1)*separation + separation-0.5*cm ;
             
            if (std::abs(xPos) >= centralHalfWidthx || std::abs(yPos) >=  centralHalfWidthy) {
               
            G4double zPos =-PhantomZ/2 -0.02*cm+PhantomZ/2-profundidadGrupo - layer * layerSeparation -cobreConectionCapacitorZ/2 -CapacitorZ/2;
 
            G4VPhysicalVolume *physcapacitor = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobreConectioncapacitor,
               "physcobreConectioncapacitor",
                logicWorld,
                false,
                ncapacitor
                
            );
           ncapacitor++;
            }

        }
    }
}


//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//* 
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*





//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--Baquelita       --//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//* 
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*







//!//!//!//!
//! COBRE//!
//!//!//!//!
G4double CobreConectionDiodeZ=0.007*cm;
//!//!//!//!
//!//!//!//!


G4double BakeliteWidth =36.2*cm ;
G4double BakeliteHeight = 20 *cm ;
G4double BakeliteThickness = 0.16*cm;

// Separación entre capacitores 

// Crear el volumen sólido del capacitor
  G4Box * solidBakelite = new G4Box("solidcapacitor", BakeliteWidth/2, BakeliteHeight/2, BakeliteThickness/2);

// Crear el volumen lógico del capacitor
  G4LogicalVolume * logicBakelite = new G4LogicalVolume(solidBakelite,mat_FR4, "logiccapacitor");
for (int i = 0; i < numLayers; i++) {
    // Calcula la posición Z de cada capa en función de `i` y `layerSeparation`
    G4double posZ = PhantomZ / 2 - profundidadGrupo - i * (layerSeparation) - BakeliteThickness / 2 - SDiodeZ - CobreConectionDiodeZ -PhantomZ/2;

    // Crear la baquelita en la posición calculada
    G4VPhysicalVolume *physBakelite = new G4PVPlacement(0,
                                                        G4ThreeVector(BakeliteWidth/2-10*cm , 0., posZ),
                                                        logicBakelite,
                                                        "physBakelite",
                                                        logicWorld,
                                                        false,
                                                        i + 1,  // ID único para cada capa
                                                        true);
}
 

//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//* 
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*




// //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^v
// //^^^^^PISTAS DEBAJO DIODO CARA SUPERIOR   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^v
// //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^v
 
//*-PISTACOBRE 1-
 G4double Pista1X = 22.0726*cm;
 G4double Pista1Y = 0.012*cm;
 G4double Pista1Z = 1*0.0035*cm;
 G4double shiftX =.9897*cm ;
 G4double shiftY =-.3264*cm ;
// Crear el volumen sólido del detector
   G4Box *  cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

// Crear el volumen lógico del detector
   G4LogicalVolume * logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
 
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) { 
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =+Pista1X/2-BakeliteHeight/2+shiftX  ;
            G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ +Pista1Z/2-CobreConectionDiodeZ  ;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
            
            G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobrePista1, "physcobrePista1",
                logicWorld, false, ndethcn  );
            ndethcn++;
        
}
} 


                            //*-PISTACOBRE 1- 1-
                             
                             
                             Pista1X =-0.10*mm+6.240*mm; //! urgentencamniar al empezar a simular
                             Pista1Y = 0.012*cm;
                             Pista1Z = 1*0.0035*cm;
                             shiftX =22.0726*cm    +.9897*cm -Pista1X/2+4.412*mm/2  ;
                             shiftY =-.3264*cm +4.412*mm/2;
                            // Crear el volumen sólido del detector
                             cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

                            // Crear el volumen lógico del detector
                              logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
                            
                            // Colocar los detectores en la cuadrícula y capas
                            for(G4int layer = 0; layer < numLayers; layer++) { 
                                    for(G4int j = 0; j < numDetectors; j++) {
                                        G4double xPos = +Pista1X/2-BakeliteHeight/2+shiftX   ;
                                        G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
                                        G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ +Pista1Z/2-CobreConectionDiodeZ  ;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                        
                                        G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                                             rotation,
                                            G4ThreeVector(xPos, yPos, zPos),
                                            logiccobrePista1, "physcobrePista1",
                                            logicWorld, false, ndethcn  );
                                        ndethcn++;
                                    
                            }
                            } 
                    //*-PISTACOBRE 1- 2-
                                
                                
                                Pista1X =8.964*mm; //! urgentencamniar al empezar a simular
                                Pista1Y = 0.012*cm;
                                Pista1Z = 1*0.0035*cm;
                                shiftX =22.0726*cm    +.9897*cm  +4.412*mm   ;
                                shiftY =-.3264*cm +4.412*mm ;
                                // Crear el volumen sólido del detector
                                cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

                                // Crear el volumen lógico del detector
                                logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
                                
                                // Colocar los detectores en la cuadrícula y capas
                                for(G4int layer = 0; layer < numLayers; layer++) { 
                                        for(G4int j = 0; j < numDetectors; j++) {
                                            G4double xPos = +Pista1X/2-BakeliteHeight/2+shiftX   ;
                                            G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
                                            G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ +Pista1Z/2-CobreConectionDiodeZ  ;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                            
                                            G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                                                0,
                                                G4ThreeVector(xPos, yPos, zPos),
                                                logiccobrePista1, "physcobrePista1",
                                                logicWorld, false, ndethcn  );
                                            ndethcn++;
                                        
                                }
                                } 

//*-PISTACOBRE 2-
  Pista1X = 204.724*mm;
  Pista1Y = 0.012*cm;
 Pista1Z = 1*0.0035*cm;
  shiftX =25.794*mm ;
  shiftY =-3.009*mm ;
// Crear el volumen sólido del detector
    cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

// Crear el volumen lógico del detector
   logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
 
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) { 
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =+Pista1X/2-BakeliteHeight/2+shiftX  ;
            G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ +Pista1Z/2-CobreConectionDiodeZ  ;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
            
            G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobrePista1, "physcobrePista1",
                logicWorld, false, ndethcn  );
            ndethcn++;
        
}
} 
 

                                        //*-PISTACOBRE 2- 1-
                                    
                                    
                                    Pista1X =-0.10*mm+6.286*mm; //! urgentencamniar al empezar a simular
                                    Pista1Y = 0.012*cm;
                                    Pista1Z = 1*0.0035*cm;
                                    shiftX =204.724*mm + 25.794*mm -Pista1X/2+4.445*mm/2  ;
                                    shiftY =-3.009*mm  +4.445*mm/2;
                                    // Crear el volumen sólido del detector
                                    cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

                                    // Crear el volumen lógico del detector
                                    logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
                                    
                                    // Colocar los detectores en la cuadrícula y capas
                                    for(G4int layer = 0; layer < numLayers; layer++) { 
                                            for(G4int j = 0; j < numDetectors; j++) {
                                                G4double xPos = +Pista1X/2-BakeliteHeight/2+shiftX   ;
                                                G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
                                                G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ +Pista1Z/2-CobreConectionDiodeZ  ;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                                
                                                G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                                                    rotation,
                                                    G4ThreeVector(xPos, yPos, zPos),
                                                    logiccobrePista1, "physcobrePista1",
                                                    logicWorld, false, ndethcn  );
                                                ndethcn++;
                                            
                                    }
                                    } 


                                                    //*-PISTACOBRE 2- 1-
                                                    
                                                    
                                                    Pista1X =6.858*mm; //! urgentencamniar al empezar a simular
                                                    Pista1Y = 0.012*cm;
                                                    Pista1Z = 1*0.0035*cm;
                                                    shiftX =204.724*mm + 25.794*mm  +4.445*mm   ;
                                                    shiftY =-3.009*mm  +4.445*mm ;
                                                    // Crear el volumen sólido del detector
                                                    cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

                                                    // Crear el volumen lógico del detector
                                                    logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
                                                    
                                                    // Colocar los detectores en la cuadrícula y capas
                                                    for(G4int layer = 0; layer < numLayers; layer++) { 
                                                            for(G4int j = 0; j < numDetectors; j++) {
                                                                G4double xPos = +Pista1X/2-BakeliteHeight/2+shiftX   ;
                                                                G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
                                                                G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ +Pista1Z/2-CobreConectionDiodeZ  ;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                                                
                                                                G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                                                                    0,
                                                                    G4ThreeVector(xPos, yPos, zPos),
                                                                    logiccobrePista1, "physcobrePista1",
                                                                    logicWorld, false, ndethcn  );
                                                                ndethcn++;
                                                            
                                                    }
                                                    } 

//*-PISTACOBRE 3-
  Pista1X = 188.634*mm;
  Pista1Y = 0.012*cm;
 Pista1Z = 1*0.0035*cm;
  shiftX =41.774*mm ;
  shiftY =-2.755*mm ;
// Crear el volumen sólido del detector
    cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

// Crear el volumen lógico del detector
   logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
 
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) { 
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =+Pista1X/2-BakeliteHeight/2+shiftX  ;
            G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation ;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ +Pista1Z/2-CobreConectionDiodeZ  ;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
            
            G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobrePista1, "physcobrePista1",
                logicWorld, false, ndethcn  );
            ndethcn++;
        
}
} 




                            //*-PISTACOBRE 3- 1-
                             
                             
                             Pista1X =-0.10*mm+6.286*mm; //! urgentencamniar al empezar a simular
                             Pista1Y = 0.012*cm;
                             Pista1Z = 1*0.0035*cm;
                             shiftX =188.634*mm + 41.774*mm -Pista1X/2+4.445*mm/2  ;
                             shiftY =-2.755*mm  +4.445*mm/2;
                            // Crear el volumen sólido del detector
                             cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

                            // Crear el volumen lógico del detector
                              logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
                            
                            // Colocar los detectores en la cuadrícula y capas
                            for(G4int layer = 0; layer < numLayers; layer++) { 
                                    for(G4int j = 0; j < numDetectors; j++) {
                                        G4double xPos = +Pista1X/2-BakeliteHeight/2+shiftX   ;
                                        G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
                                        G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ +Pista1Z/2-CobreConectionDiodeZ  ;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                        
                                        G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                                             rotation,
                                            G4ThreeVector(xPos, yPos, zPos),
                                            logiccobrePista1, "physcobrePista1",
                                            logicWorld, false, ndethcn  );
                                        ndethcn++;
                                    
                            }
                            } 



                        //*-PISTACOBRE 3- 1-
                                                    
                                                    
                                                    Pista1X =6.819*mm; //! urgentencamniar al empezar a simular
                                                    Pista1Y = 0.012*cm;
                                                    Pista1Z = 1*0.0035*cm;
                                                    shiftX =188.634*mm + 41.774*mm +4.445*mm   ;
                                                    shiftY =-2.755*mm  +4.445*mm ;
                                                    // Crear el volumen sólido del detector
                                                    cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

                                                    // Crear el volumen lógico del detector
                                                    logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
                                                    
                                                    // Colocar los detectores en la cuadrícula y capas
                                                    for(G4int layer = 0; layer < numLayers; layer++) { 
                                                            for(G4int j = 0; j < numDetectors; j++) {
                                                                G4double xPos = +Pista1X/2-BakeliteHeight/2+shiftX   ;
                                                                G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
                                                                G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ +Pista1Z/2-CobreConectionDiodeZ  ;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                                                
                                                                G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                                                                    0,
                                                                    G4ThreeVector(xPos, yPos, zPos),
                                                                    logiccobrePista1, "physcobrePista1",
                                                                    logicWorld, false, ndethcn  );
                                                                ndethcn++;
                                                            
                                                    }
                                                    } 








//*-PISTACOBRE 4-
  Pista1X = 172.417*mm;
  Pista1Y = 0.012*cm;
 Pista1Z = 1*0.0035*cm;
  shiftX =57.903*mm ;
  shiftY =-2.501*mm ;
// Crear el volumen sólido del detector
    cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

// Crear el volumen lógico del detector
   logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
 
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) { 
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =+Pista1X/2-BakeliteHeight/2+shiftX  ;
            G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation ;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ +Pista1Z/2-CobreConectionDiodeZ  ;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
            
            G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobrePista1, "physcobrePista1",
                logicWorld, false, ndethcn  );
            ndethcn++;
        
}
} 


//*-PISTACOBRE 4- 1-
                             
                             
                             Pista1X =-0.10*mm+6.286*mm; //! urgentencamniar al empezar a simular
                             Pista1Y = 0.012*cm;
                             Pista1Z = 1*0.0035*cm;
                             shiftX =172.417*mm + 57.903*mm -Pista1X/2+4.445*mm/2  ;
                             shiftY =-2.501*mm  +4.445*mm/2;
                            // Crear el volumen sólido del detector
                             cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

                            // Crear el volumen lógico del detector
                              logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
                            
                            // Colocar los detectores en la cuadrícula y capas
                            for(G4int layer = 0; layer < numLayers; layer++) { 
                                    for(G4int j = 0; j < numDetectors; j++) {
                                        G4double xPos = +Pista1X/2-BakeliteHeight/2+shiftX   ;
                                        G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
                                        G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ +Pista1Z/2-CobreConectionDiodeZ  ;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                        
                                        G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                                             rotation,
                                            G4ThreeVector(xPos, yPos, zPos),
                                            logiccobrePista1, "physcobrePista1",
                                            logicWorld, false, ndethcn  );
                                        ndethcn++;
                                    
                            }
                            } 





        //*-PISTACOBRE 4- 1-
                                    
                                    
                                    Pista1X =6.819*mm; //! urgentencamniar al empezar a simular
                                    Pista1Y = 0.012*cm;
                                    Pista1Z = 1*0.0035*cm;
                                    shiftX =172.417*mm + 57.903*mm  +4.445*mm   ;
                                    shiftY =-2.501*mm  +4.445*mm ;
                                    // Crear el volumen sólido del detector
                                    cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

                                    // Crear el volumen lógico del detector
                                    logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
                                    
                                    // Colocar los detectores en la cuadrícula y capas
                                    for(G4int layer = 0; layer < numLayers; layer++) { 
                                            for(G4int j = 0; j < numDetectors; j++) {
                                                G4double xPos = +Pista1X/2-BakeliteHeight/2+shiftX   ;
                                                G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
                                                G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ +Pista1Z/2-CobreConectionDiodeZ  ;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                                
                                                G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                                                    0,
                                                    G4ThreeVector(xPos, yPos, zPos),
                                                    logiccobrePista1, "physcobrePista1",
                                                    logicWorld, false, ndethcn  );
                                                ndethcn++;
                                            
                                    }
                                    } 











//*-PISTACOBRE 5-
  Pista1X = 156.454*mm;
  Pista1Y = 0.012*cm;
 Pista1Z = 1*0.0035*cm;
  shiftX =73.778*mm ;
  shiftY =-2.247*mm ;
// Crear el volumen sólido del detector
    cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

// Crear el volumen lógico del detector
   logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
 
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) { 
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =+Pista1X/2-BakeliteHeight/2+shiftX  ;
            G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation ;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ +Pista1Z/2-CobreConectionDiodeZ  ;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
            
            G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobrePista1, "physcobrePista1",
                logicWorld, false, ndethcn  );
            ndethcn++;
        
}
} 



//*-PISTACOBRE 5- 1-
                             
                             
                             Pista1X =-0.10*mm+6.286*mm; //! urgentencamniar al empezar a simular
                             Pista1Y = 0.012*cm;
                             Pista1Z = 1*0.0035*cm;
                             shiftX =156.454*mm +73.778*mm -Pista1X/2+4.445*mm/2  ;
                             shiftY =-2.247*mm  +4.445*mm/2;
                            // Crear el volumen sólido del detector
                             cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

                            // Crear el volumen lógico del detector
                              logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
                            
                            // Colocar los detectores en la cuadrícula y capas
                            for(G4int layer = 0; layer < numLayers; layer++) { 
                                    for(G4int j = 0; j < numDetectors; j++) {
                                        G4double xPos = +Pista1X/2-BakeliteHeight/2+shiftX   ;
                                        G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
                                        G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ +Pista1Z/2-CobreConectionDiodeZ  ;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                        
                                        G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                                             rotation,
                                            G4ThreeVector(xPos, yPos, zPos),
                                            logiccobrePista1, "physcobrePista1",
                                            logicWorld, false, ndethcn  );
                                        ndethcn++;
                                    
                            }
                            } 




            //*-PISTACOBRE 5- 2-
                                        
                                        
                                        Pista1X =6.819*mm; //! urgentencamniar al empezar a simular
                                        Pista1Y = 0.012*cm;
                                        Pista1Z = 1*0.0035*cm;
                                        shiftX =156.454*mm +73.778*mm +4.445*mm   ;
                                        shiftY =-2.247*mm  +4.445*mm  ;
                                        // Crear el volumen sólido del detector
                                        cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

                                        // Crear el volumen lógico del detector
                                        logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
                                        
                                        // Colocar los detectores en la cuadrícula y capas
                                        for(G4int layer = 0; layer < numLayers; layer++) { 
                                                for(G4int j = 0; j < numDetectors; j++) {
                                                    G4double xPos = +Pista1X/2-BakeliteHeight/2+shiftX   ;
                                                    G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
                                                    G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ +Pista1Z/2-CobreConectionDiodeZ  ;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                                    
                                                    G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                                                        0,
                                                        G4ThreeVector(xPos, yPos, zPos),
                                                        logiccobrePista1, "physcobrePista1",
                                                        logicWorld, false, ndethcn  );
                                                    ndethcn++;
                                                
                                        }
                                        } 




//*-PISTACOBRE 6-
  Pista1X = 140.491*mm;
  Pista1Y = 0.012*cm;
 Pista1Z = 1*0.0035*cm;
  shiftX =89.653*mm ;
  shiftY =-1.993*mm ;
// Crear el volumen sólido del detector
    cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

// Crear el volumen lógico del detector
   logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
 
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) { 
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =+Pista1X/2-BakeliteHeight/2+shiftX  ;
            G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation ;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ +Pista1Z/2-CobreConectionDiodeZ  ;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
            
            G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobrePista1, "physcobrePista1",
                logicWorld, false, ndethcn  );
            ndethcn++;
        
}
} 





//*-PISTACOBRE 6- 1-
                             
                             
                             Pista1X =-0.10*mm+6.286*mm; //! urgentencamniar al empezar a simular
                             Pista1Y = 0.012*cm;
                             Pista1Z = 1*0.0035*cm;
                             shiftX =140.491*mm +89.653*mm -Pista1X/2+4.445*mm/2  ;
                             shiftY =-1.993*mm  +4.445*mm/2;
                            // Crear el volumen sólido del detector
                             cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

                            // Crear el volumen lógico del detector
                              logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
                            
                            // Colocar los detectores en la cuadrícula y capas
                            for(G4int layer = 0; layer < numLayers; layer++) { 
                                    for(G4int j = 0; j < numDetectors; j++) {
                                        G4double xPos = +Pista1X/2-BakeliteHeight/2+shiftX   ;
                                        G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
                                        G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ +Pista1Z/2-CobreConectionDiodeZ  ;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                        
                                        G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                                             rotation,
                                            G4ThreeVector(xPos, yPos, zPos),
                                            logiccobrePista1, "physcobrePista1",
                                            logicWorld, false, ndethcn  );
                                        ndethcn++;
                                    
                            }
                            } 


            //*-PISTACOBRE 6- 1-
                                        
                                        
                                        Pista1X =6.848*mm; //! urgentencamniar al empezar a simular
                                        Pista1Y = 0.012*cm;
                                        Pista1Z = 1*0.0035*cm;
                                        shiftX =140.491*mm +89.653*mm +4.445*mm   ;
                                        shiftY =-1.993*mm  +4.445*mm ;
                                        // Crear el volumen sólido del detector
                                        cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

                                        // Crear el volumen lógico del detector
                                        logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
                                        
                                        // Colocar los detectores en la cuadrícula y capas
                                        for(G4int layer = 0; layer < numLayers; layer++) { 
                                                for(G4int j = 0; j < numDetectors; j++) {
                                                    G4double xPos = +Pista1X/2-BakeliteHeight/2+shiftX   ;
                                                    G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
                                                    G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ +Pista1Z/2-CobreConectionDiodeZ  ;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                                    
                                                    G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                                                        0,
                                                        G4ThreeVector(xPos, yPos, zPos),
                                                        logiccobrePista1, "physcobrePista1",
                                                        logicWorld, false, ndethcn  );
                                                    ndethcn++;
                                                
                                        }
                                        } 












 

 


// //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^v
// //^^^^^PISTAS CARA INFERIOR   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^v
// //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^v
 
//*-PISTACOBRE 0-
 Pista1X = 176*mm;
  Pista1Y = .2*mm;
  Pista1Z = 1*0.0035*cm;
 shiftX =14.342*mm ;
 shiftY =-4*mm ;
// Crear el volumen sólido del detector
   cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

// Crear el volumen lógico del detector
   logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
 
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) { 
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =+Pista1X/2-BakeliteHeight/2+shiftX  ;
            G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ -CobreConectionDiodeZ - BakeliteThickness -Pista1Z/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
            
            G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobrePista1, "physcobrePista1",
                logicWorld, false, ndethcn  );
            ndethcn++;
        
}
} 


//*-PISTACOBRE 1-
 Pista1X =207.91*mm;
  Pista1Y = .12*mm;
  Pista1Z = 1*0.0035*cm;
 shiftX =22.713*mm ;
 shiftY =-3.264*mm ;
// Crear el volumen sólido del detector
   cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

// Crear el volumen lógico del detector
   logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
 
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) { 
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =+Pista1X/2-BakeliteHeight/2+shiftX  ;
            G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ -CobreConectionDiodeZ - BakeliteThickness -Pista1Z/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
            
            G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobrePista1, "physcobrePista1",
                logicWorld, false, ndethcn  );
            ndethcn++;
        
}
} 


                                //*-Palitodoblado 1-
                                Pista1X =8.151*mm;
                                Pista1Y = .12*mm;
                                Pista1Z = 1*0.0035*cm;
                                shiftX =22.713*mm -207.91*mm/2+5.764*mm/2;
                                shiftY =-3.264*mm+5.764*mm/2 +0.12*mm;
                                // Crear el volumen sólido del detector
                                cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

                                // Crear el volumen lógico del detector
                                logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
                                
                                // Colocar los detectores en la cuadrícula y capas
                                for(G4int layer = 0; layer < numLayers; layer++) { 
                                        for(G4int j = 0; j < numDetectors; j++) {
                                            G4double xPos =+207.91*mm/2-BakeliteHeight/2+shiftX  ;
                                            G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
                                            G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ -CobreConectionDiodeZ - BakeliteThickness -Pista1Z/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                            
                                            G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                                                rotation,
                                                G4ThreeVector(xPos, yPos, zPos),
                                                logiccobrePista1, "physcobrePista1",
                                                logicWorld, false, ndethcn  );
                                            ndethcn++;
                                        
                                }
                                } 

                                



//*-PISTACOBRE 2-
 Pista1X =(188.341+3.569)*mm;
  Pista1Y = .12*mm;
  Pista1Z = 1*0.0035*cm;
 shiftX =38.459*mm ;
 shiftY =-3.010*mm ;
// Crear el volumen sólido del detector
   cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

// Crear el volumen lógico del detector
   logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
 
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) { 
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =+Pista1X/2-BakeliteHeight/2+shiftX  ;
            G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ -CobreConectionDiodeZ - BakeliteThickness -Pista1Z/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
            
            G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobrePista1, "physcobrePista1",
                logicWorld, false, ndethcn  );
            ndethcn++;
        
}
} 







                                //*-Palitodoblado 2-
                                Pista1X =7.791*mm;
                                Pista1Y = .12*mm;
                                Pista1Z = 1*0.0035*cm;
                                shiftX =38.459*mm -(188.341+3.569)*mm/2+5.764*mm/2;
                                shiftY =-3.010*mm+5.764*mm/2 +0.12*mm;
                                // Crear el volumen sólido del detector
                                cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

                                // Crear el volumen lógico del detector
                                logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
                                
                                // Colocar los detectores en la cuadrícula y capas
                                for(G4int layer = 0; layer < numLayers; layer++) { 
                                        for(G4int j = 0; j < numDetectors; j++) {
                                            G4double xPos =+(188.341+3.569)*mm/2-BakeliteHeight/2+shiftX  ;
                                            G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
                                            G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ -CobreConectionDiodeZ - BakeliteThickness -Pista1Z/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                            
                                            G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                                                rotation,
                                                G4ThreeVector(xPos, yPos, zPos),
                                                logiccobrePista1, "physcobrePista1",
                                                logicWorld, false, ndethcn  );
                                            ndethcn++;
                                        
                                }
                                } 





//*-PISTACOBRE 3-
 Pista1X =(176.242)*mm;
  Pista1Y = .12*mm;
  Pista1Z = 1*0.0035*cm;
 shiftX =54.219*mm ;
 shiftY =-2.756*mm ;
// Crear el volumen sólido del detector
   cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

// Crear el volumen lógico del detector
   logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
 
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) { 
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =+Pista1X/2-BakeliteHeight/2+shiftX  ;
            G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ -CobreConectionDiodeZ - BakeliteThickness -Pista1Z/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
            
            G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobrePista1, "physcobrePista1",
                logicWorld, false, ndethcn  );
            ndethcn++;
        
}
} 









                                //*-Palitodoblado 3-
                                Pista1X =7.433*mm;
                                Pista1Y = .12*mm;
                                Pista1Z = 1*0.0035*cm;
                                shiftX =54.219*mm -(176.242)*mm/2+5.256*mm/2;
                                shiftY =-2.756*mm+5.256*mm/2 +0.12*mm;
                                // Crear el volumen sólido del detector
                                cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

                                // Crear el volumen lógico del detector
                                logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
                                
                                // Colocar los detectores en la cuadrícula y capas
                                for(G4int layer = 0; layer < numLayers; layer++) { 
                                        for(G4int j = 0; j < numDetectors; j++) {
                                            G4double xPos =+(176.242)*mm/2-BakeliteHeight/2+shiftX  ;
                                            G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
                                            G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ -CobreConectionDiodeZ - BakeliteThickness -Pista1Z/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                            
                                            G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                                                rotation,
                                                G4ThreeVector(xPos, yPos, zPos),
                                                logiccobrePista1, "physcobrePista1",
                                                logicWorld, false, ndethcn  );
                                            ndethcn++;
                                        
                                }
                                } 











//*-PISTACOBRE 4-
 Pista1X =(160.408)*mm;
  Pista1Y = .12*mm;
  Pista1Z = 1*0.0035*cm;
 shiftX =69.951*mm ;
 shiftY =-2.502*mm ;
// Crear el volumen sólido del detector
   cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

// Crear el volumen lógico del detector
   logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
 
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) { 
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =+Pista1X/2-BakeliteHeight/2+shiftX  ;
            G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ -CobreConectionDiodeZ - BakeliteThickness -Pista1Z/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
            
            G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobrePista1, "physcobrePista1",
                logicWorld, false, ndethcn  );
            ndethcn++;
        
}
} 




                                //*-Palitodoblado 4-
                                Pista1X =7.073*mm;
                                Pista1Y = .12*mm;
                                Pista1Z = 1*0.0035*cm;
                                shiftX =69.951*mm -(160.408)*mm/2+5.002*mm/2;
                                shiftY =-2.502*mm+5.002*mm/2 +0.12*mm;
                                // Crear el volumen sólido del detector
                                cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

                                // Crear el volumen lógico del detector
                                logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
                                
                                // Colocar los detectores en la cuadrícula y capas
                                for(G4int layer = 0; layer < numLayers; layer++) { 
                                        for(G4int j = 0; j < numDetectors; j++) {
                                            G4double xPos =+(160.408)*mm/2-BakeliteHeight/2+shiftX  ;
                                            G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
                                            G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ -CobreConectionDiodeZ - BakeliteThickness -Pista1Z/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                            
                                            G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                                                rotation,
                                                G4ThreeVector(xPos, yPos, zPos),
                                                logiccobrePista1, "physcobrePista1",
                                                logicWorld, false, ndethcn  );
                                            ndethcn++;
                                        
                                }
                                } 







//*-PISTACOBRE 5-
 Pista1X =(144.575)*mm;
  Pista1Y = .12*mm;
  Pista1Z = 1*0.0035*cm;
 shiftX =85.697*mm ;
 shiftY =-2.248*mm ;
// Crear el volumen sólido del detector
   cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

// Crear el volumen lógico del detector
   logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
 
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) { 
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =+Pista1X/2-BakeliteHeight/2+shiftX  ;
            G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ -CobreConectionDiodeZ - BakeliteThickness -Pista1Z/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
            
            G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobrePista1, "physcobrePista1",
                logicWorld, false, ndethcn  );
            ndethcn++;
        
}
} 



                                //*-Palitodoblado 5-
                                Pista1X =6.714*mm;
                                Pista1Y = .12*mm;
                                Pista1Z = 1*0.0035*cm;
                                shiftX =85.697*mm -(144.575)*mm/2+4.747*mm/2;
                                shiftY =-2.248*mm+4.747*mm/2 +0.12*mm;
                                // Crear el volumen sólido del detector
                                cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

                                // Crear el volumen lógico del detector
                                logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
                                
                                // Colocar los detectores en la cuadrícula y capas
                                for(G4int layer = 0; layer < numLayers; layer++) { 
                                        for(G4int j = 0; j < numDetectors; j++) {
                                            G4double xPos =+(144.575)*mm/2-BakeliteHeight/2+shiftX  ;
                                            G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
                                            G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ -CobreConectionDiodeZ - BakeliteThickness -Pista1Z/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                            
                                            G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                                                rotation,
                                                G4ThreeVector(xPos, yPos, zPos),
                                                logiccobrePista1, "physcobrePista1",
                                                logicWorld, false, ndethcn  );
                                            ndethcn++;
                                        
                                }
                                } 



//*-PISTACOBRE 6-
 Pista1X =(124.587+4.085)*mm;
  Pista1Y = .12*mm;
  Pista1Z = 1*0.0035*cm;
 shiftX =101.443*mm ;
 shiftY =-1.994*mm ;
// Crear el volumen sólido del detector
   cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

// Crear el volumen lógico del detector
   logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
 
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) { 
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =+Pista1X/2-BakeliteHeight/2+shiftX  ;
            G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ -CobreConectionDiodeZ - BakeliteThickness -Pista1Z/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
            
            G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobrePista1, "physcobrePista1",
                logicWorld, false, ndethcn  );
            ndethcn++;
        
}
} 

 //*-Palitodoblado 6-
                                Pista1X =6.355*mm;
                                Pista1Y = .12*mm;
                                Pista1Z = 1*0.0035*cm;
                                shiftX =101.443*mm -(124.587+4.085)*mm/2+4.494*mm/2;
                                shiftY =-1.994*mm+4.494*mm/2 +0.12*mm;
                                // Crear el volumen sólido del detector
                                cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

                                // Crear el volumen lógico del detector
                                logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
                                
                                // Colocar los detectores en la cuadrícula y capas
                                for(G4int layer = 0; layer < numLayers; layer++) { 
                                        for(G4int j = 0; j < numDetectors; j++) {
                                            G4double xPos =+(124.587+4.085)*mm/2-BakeliteHeight/2+shiftX  ;
                                            G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
                                            G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ -CobreConectionDiodeZ - BakeliteThickness -Pista1Z/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                            
                                            G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                                                rotation,
                                                G4ThreeVector(xPos, yPos, zPos),
                                                logiccobrePista1, "physcobrePista1",
                                                logicWorld, false, ndethcn  );
                                            ndethcn++;
                                        
                                }
                                } 

//*-PISTACOBRE 7-
 Pista1X =(116.752+4.086)*mm;
  Pista1Y = .12*mm;
  Pista1Z = 1*0.0035*cm;
 shiftX =109.189*mm ;
 shiftY =-1.740*mm ;
// Crear el volumen sólido del detector
   cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

// Crear el volumen lógico del detector
   logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
 
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) { 
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =+Pista1X/2-BakeliteHeight/2+shiftX  ;
            G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ -CobreConectionDiodeZ - BakeliteThickness -Pista1Z/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
            
            G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobrePista1, "physcobrePista1",
                logicWorld, false, ndethcn  );
            ndethcn++;
        
}
} 

 //*-Palitodoblado 6-
                                Pista1X =5.996*mm;
                                Pista1Y = .12*mm;
                                Pista1Z = 1*0.0035*cm;
                                shiftX =109.189*mm -(116.752+4.086)*mm/2+4.240*mm/2;
                                shiftY =-1.740*mm+4.240*mm/2 +0.12*mm;
                                // Crear el volumen sólido del detector
                                cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

                                // Crear el volumen lógico del detector
                                logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
                                
                                // Colocar los detectores en la cuadrícula y capas
                                for(G4int layer = 0; layer < numLayers; layer++) { 
                                        for(G4int j = 0; j < numDetectors; j++) {
                                            G4double xPos =+(116.752+4.086)*mm/2-BakeliteHeight/2+shiftX  ;
                                            G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
                                            G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ -CobreConectionDiodeZ - BakeliteThickness -Pista1Z/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
                                            
                                            G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                                                rotation,
                                                G4ThreeVector(xPos, yPos, zPos),
                                                logiccobrePista1, "physcobrePista1",
                                                logicWorld, false, ndethcn  );
                                            ndethcn++;
                                        
                                }
                                } 

//*-PISTACOBRE 8-
 Pista1X =(94.511+18.493)*mm;
  Pista1Y = .12*mm;
  Pista1Z = 1*0.0035*cm;
 shiftX =116.935*mm ;
 shiftY =-1.486*mm ;
// Crear el volumen sólido del detector
   cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

// Crear el volumen lógico del detector
   logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
 
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) { 
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =+Pista1X/2-BakeliteHeight/2+shiftX  ;
            G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ -CobreConectionDiodeZ - BakeliteThickness -Pista1Z/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
            
            G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobrePista1, "physcobrePista1",
                logicWorld, false, ndethcn  );
            ndethcn++;
        
}
} 






//*-PISTACOBRE 9-
 Pista1X =(12.725+62.738+29.708)*mm;
  Pista1Y = .12*mm;
  Pista1Z = 1*0.0035*cm;
 shiftX =124.680*mm ;
 shiftY =-1.232*mm ;
// Crear el volumen sólido del detector
   cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

// Crear el volumen lógico del detector
   logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
 
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) { 
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =+Pista1X/2-BakeliteHeight/2+shiftX  ;
            G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ -CobreConectionDiodeZ - BakeliteThickness -Pista1Z/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
            
            G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobrePista1, "physcobrePista1",
                logicWorld, false, ndethcn  );
            ndethcn++;
        
}
} 





//*-PISTACOBRE 10-
 Pista1X =(12.725+ 59.309+25.273)*mm;
  Pista1Y = .12*mm;
  Pista1Z = 1*0.0035*cm;
 shiftX =132.427*mm ;
 shiftY =-0.978*mm ;
// Crear el volumen sólido del detector
   cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

// Crear el volumen lógico del detector
   logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
 
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) { 
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =+Pista1X/2-BakeliteHeight/2+shiftX  ;
            G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ -CobreConectionDiodeZ - BakeliteThickness -Pista1Z/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
            
            G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobrePista1, "physcobrePista1",
                logicWorld, false, ndethcn  );
            ndethcn++;
        
}
} 







//*-PISTACOBRE 11-
 Pista1X =(16.027+49.657+23.788)*mm;
  Pista1Y = .12*mm;
  Pista1Z = 1*0.0035*cm;
 shiftX =140.174*mm ;
 shiftY =-0.724*mm ;
// Crear el volumen sólido del detector
   cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

// Crear el volumen lógico del detector
   logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
 
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) { 
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =+Pista1X/2-BakeliteHeight/2+shiftX  ;
            G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ -CobreConectionDiodeZ - BakeliteThickness -Pista1Z/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
            
            G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobrePista1, "physcobrePista1",
                logicWorld, false, ndethcn  );
            ndethcn++;
        
}
} 



//*-PISTACOBRE 12-
 Pista1X =(8.155+53.213+20.271)*mm;
  Pista1Y = .12*mm;
  Pista1Z = 1*0.0035*cm;
 shiftX =147.919*mm ;
 shiftY =-0.470*mm ;
// Crear el volumen sólido del detector
   cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

// Crear el volumen lógico del detector
   logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
 
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) { 
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =+Pista1X/2-BakeliteHeight/2+shiftX  ;
            G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ -CobreConectionDiodeZ - BakeliteThickness -Pista1Z/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
            
            G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobrePista1, "physcobrePista1",
                logicWorld, false, ndethcn  );
            ndethcn++;
        
}
} 



//*-PISTACOBRE 13-
 Pista1X =(9.045+37.303+27.549)*mm;
  Pista1Y = .12*mm;
  Pista1Z = 1*0.0035*cm;
 shiftX =155.665*mm ;
 shiftY =-0.216*mm ;
// Crear el volumen sólido del detector
   cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

// Crear el volumen lógico del detector
   logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
 
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) { 
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =+Pista1X/2-BakeliteHeight/2+shiftX  ;
            G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ -CobreConectionDiodeZ - BakeliteThickness -Pista1Z/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
            
            G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobrePista1, "physcobrePista1",
                logicWorld, false, ndethcn  );
            ndethcn++;
        
}
} 


//*-PISTACOBRE 14-
 Pista1X =(15.905+25.654+24.384)*mm;
  Pista1Y = .12*mm;
  Pista1Z = 1*0.0035*cm;
 shiftX =163.410*mm ;
 shiftY =0.038*mm ;
// Crear el volumen sólido del detector
   cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

// Crear el volumen lógico del detector
   logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
 
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) { 
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =+Pista1X/2-BakeliteHeight/2+shiftX  ;
            G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ -CobreConectionDiodeZ - BakeliteThickness -Pista1Z/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
            
            G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobrePista1, "physcobrePista1",
                logicWorld, false, ndethcn  );
            ndethcn++;
        
}
} 




//*-PISTACOBRE 15-
 Pista1X =( 8.158+23.241+26.709)*mm;
  Pista1Y = .12*mm;
  Pista1Z = 1*0.0035*cm;
 shiftX =171.420*mm ;
 shiftY =0.292*mm ;
// Crear el volumen sólido del detector
   cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

// Crear el volumen lógico del detector
   logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
 
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) { 
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =+Pista1X/2-BakeliteHeight/2+shiftX  ;
            G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ -CobreConectionDiodeZ - BakeliteThickness -Pista1Z/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
            
            G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobrePista1, "physcobrePista1",
                logicWorld, false, ndethcn  );
            ndethcn++;
        
}
} 




//*-PISTACOBRE 15-
 Pista1X =(13.873+8.636+27.764)*mm;
  Pista1Y = .12*mm;
  Pista1Z = 1*0.0035*cm;
 shiftX =179.040*mm ;
 shiftY =0.546*mm ;
// Crear el volumen sólido del detector
   cobrePista1 = new G4Box("cobreConectionsolidDetector", Pista1X/2, Pista1Y/2, Pista1Z/2);

// Crear el volumen lógico del detector
   logiccobrePista1 = new G4LogicalVolume(cobrePista1,cobre, "logiccobrePista1");
 
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) { 
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =+Pista1X/2-BakeliteHeight/2+shiftX  ;
            G4double yPos =shiftY -0.5*(numDetectors-1)*separation + j*separation  ;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) -SDiodeZ -CobreConectionDiodeZ - BakeliteThickness -Pista1Z/2;// SDiodeZ SE PODRA CAMBIAR A MITAD DEE ALTURA DE COBERTURA DIODO
            
            G4VPhysicalVolume *physcobrePista1 = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logiccobrePista1, "physcobrePista1",
                logicWorld, false, ndethcn  );
            ndethcn++;
        
}
}  






//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--MOSFET   -----//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//* 
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*
//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*//*--------------------//*



//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^v
//^^^^ PISTAS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^v
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^v




// //!!!!
// //????????????
// //* PAD MOSFET
//  //!!!!
// //????????????
//&      ------------ centro diodo exterior + posicion pad + grosorPAD/2+ 5cm+ mitadmosfet: 
G4double shiftmosfeX=8.8*cm +55.874*mm;//             +   0.305*cm        +0.12*cm/2 + 5*cm+0.28*cm/2;
// //*-INTERNOS-
 G4double COBREMOSFETX =0.08  *cm;
 G4double COBREMOSFETY = 0.1*cm;
 G4double COBREMOSFETZ = 0.007*cm;
 G4int nmosfet=0;
// Crear el volumen sólido del detector
  G4Box *  COBREMOSFET = new G4Box("COBREMOSFET", COBREMOSFETX/2, COBREMOSFETY/2, COBREMOSFETZ/2);

// Crear el volumen lógico del detector
 G4LogicalVolume *  logicCOBREMOSFET = new G4LogicalVolume(COBREMOSFET,cobre, "logicCOBREMOSFET");
  
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) {
    for(G4int i = 0; i < numDetectors; i++) {
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =shiftmosfeX  + i*separation/2;
            G4double yPos =0.115*cm -0.5*(numDetectors-1)*separation + j*separation;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) - SDiodeZ-COBREMOSFETZ/2;//se puede cambiar zdiode tambien xd
            
            G4VPhysicalVolume *physCOBREMOSFET = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logicCOBREMOSFET,
                "physCOBREMOSFET",
                logicWorld,
                false,
                nmosfet
                
            );
           nmosfet++;
        }
}
}

// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) {
    for(G4int i = 0; i < numDetectors; i++) {
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos = -0.095*cm+shiftmosfeX   + i*separation/2;
            G4double yPos =-0.23*cm+0.115*cm -0.5*(numDetectors-1)*separation + j*separation;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) - SDiodeZ-COBREMOSFETZ/2;//se puede cambiar zdiode tambien xd
            
            G4VPhysicalVolume *physCOBREMOSFET = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logicCOBREMOSFET,
                "physCOBREMOSFET",
                logicWorld,
                false,
                nmosfet
                
            );
           nmosfet++;
        }
}
}





// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) {
    for(G4int i = 0; i < numDetectors; i++) {
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos = +0.095*cm+shiftmosfeX   + i*separation/2;
            G4double yPos =-0.23*cm+0.115*cm -0.5*(numDetectors-1)*separation + j*separation;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) - SDiodeZ-COBREMOSFETZ/2;//se puede cambiar zdiode tambien xd
            
            G4VPhysicalVolume *physCOBREMOSFET = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logicCOBREMOSFET,
                "physCOBREMOSFET",
                logicWorld,
                false,
                nmosfet
                
            );
           nmosfet++;
        }
}
}





// //!!!!
// //????????????
// //* PATAS MOSFET
//  //!!!!
// //????????????



// //*-INTERNOS-
 G4double PCsBX =0.043  *cm;
 G4double PCsBY = 0.043*cm;
 G4double PCsBZ = 0.014*cm;
  
 
 G4int count_PCsB=0;
// Crear el volumen sólido del detector
  G4Box * PCsB = new G4Box("PCsB", PCsBX/2, PCsBY/2, PCsBZ/2);

// Crear el volumen lógico del detector
 G4LogicalVolume *  logicPCsB = new G4LogicalVolume(PCsB,cobre, "logicPCsB");
  
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) {
    for(G4int i = 0; i < numDetectors; i++) {
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos = shiftmosfeX   + i*separation/2;
            G4double yPos =0.019*cm +0.115*cm -0.5*(numDetectors-1)*separation + j*separation;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) - SDiodeZ+PCsBZ/2;//se puede cambiar zdiode tambien xd
            
            G4VPhysicalVolume *physPCsB = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logicPCsB,
                "physPCsB",
                logicWorld,
                false,
                nmosfet
                
            );
           nmosfet++;
        }
}
}

// //*-INTERNOS-
 G4double PCsMX =0.043  *cm;
 G4double PCsMZ = 0.0511*cm; 
 G4double PCsMY = 0.014*cm;
  
 
 G4int count_PCsM=0;
// Crear el volumen sólido del detector
  G4Box * PCsM = new G4Box("PCsM", PCsMX/2, PCsMY/2, PCsMZ/2);

// Crear el volumen lógico del detector
 G4LogicalVolume *  logicPCsM = new G4LogicalVolume(PCsM,cobre, "logicPCsM");
  
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) {
    for(G4int i = 0; i < numDetectors; i++) {
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos = shiftmosfeX  + i*separation/2;
            G4double yPos =0.019*cm +0.115*cm -0.5*(numDetectors-1)*separation + j*separation-PCsBY/2+PCsMY/2;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) - SDiodeZ+PCsMZ/2+PCsBZ;//se puede cambiar zdiode tambien xd
            
            G4VPhysicalVolume *physPCsM = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logicPCsM,
                "physPCsM",
                logicWorld,
                false,
                nmosfet
                
            );
           nmosfet++;
        }
}
}


 G4double PCsTX =0.043  *cm;
 G4double PCsTY = 0.055*cm;
 G4double PCsTZ = 0.014*cm;
  
 
 G4int count_PCsT=0;
// Crear el volumen sólido del detector
  G4Box * PCsT = new G4Box("PCsT", PCsTX/2, PCsTY/2, PCsTZ/2);

// Crear el volumen lógico del detector
 G4LogicalVolume *  logicPCsT = new G4LogicalVolume(PCsT,cobre, "logicPCsT");
  
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) {
    for(G4int i = 0; i < numDetectors; i++) {
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos = shiftmosfeX   + i*separation/2;
            G4double yPos =0.019*cm +0.115*cm -0.5*(numDetectors-1)*separation + j*separation-PCsTX;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) - SDiodeZ+PCsBZ/2+PCsMZ;//se puede cambiar zdiode tambien xd
            
            G4VPhysicalVolume *physPCsB = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logicPCsB,
                "physPCsB",
                logicWorld,
                false,
                nmosfet
                
            );
           nmosfet++;
        }
}
}















 
// //*-INTERNOS-
 G4double PCiiBX =0.043  *cm;
 G4double PCiiBY = 0.043*cm;
 G4double PCiiBZ = 0.014*cm;
  
 
 G4int count_PCiiB=0;
// Crear el volumen sólido del detector
  G4Box * PCiiB = new G4Box("PCiiB", PCsBX/2, PCsBY/2, PCsBZ/2);

// Crear el volumen lógico del detector
 G4LogicalVolume *  logicPCiiB = new G4LogicalVolume(PCiiB,cobre, "logicPCiiB");
  
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) {
    for(G4int i = 0; i < numDetectors; i++) {
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos = -0.095*cm+shiftmosfeX  + i*separation/2;
            G4double yPos =0.019*cm -0.23*cm+0.115*cm -0.5*(numDetectors-1)*separation + j*separation-PCiiBY+PCiiBZ;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) - SDiodeZ+PCiiBZ/2;//se puede cambiar zdiode tambien xd
            
            G4VPhysicalVolume *physPCiiB = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logicPCiiB,
                "physPCiiB",
                logicWorld,
                false,
                nmosfet
                
            );
           nmosfet++;
        }
}
}
for(G4int layer = 0; layer < numLayers; layer++) {
    for(G4int i = 0; i < numDetectors; i++) {
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos = +0.095*cm+shiftmosfeX  + i*separation/2;
            G4double yPos =0.019*cm -0.23*cm+0.115*cm -0.5*(numDetectors-1)*separation + j*separation-PCiiBY+PCiiBZ;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) - SDiodeZ+PCiiBZ/2;//se puede cambiar zdiode tambien xd
            
            G4VPhysicalVolume *physPCiiB = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logicPCiiB,
                "physPCiiB",
                logicWorld,
                false,
                nmosfet
                
            );
           nmosfet++;
        }
}
}
// //*-INTERNOS-
 G4double PCiiMX =0.043  *cm;
 G4double PCiiMZ = 0.0511*cm; 
 G4double PCiiMY = 0.014*cm;
  
 
 G4int count_PCiiM=0;
// Crear el volumen sólido del detector
  G4Box * PCiiM = new G4Box("PCiiM", PCiiMX/2, PCiiMY/2, PCiiMZ/2);

// Crear el volumen lógico del detector
 G4LogicalVolume *  logicPCiiM = new G4LogicalVolume(PCiiM,cobre, "logicPCiiM");
  
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) {
    for(G4int i = 0; i < numDetectors; i++) {
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =-0.095*cm+ shiftmosfeX   + i*separation/2;
            G4double yPos =0.019*cm -0.23*cm+0.115*cm -0.5*(numDetectors-1)*separation + j*separation-PCiiBY/2+PCiiMY/2;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) - SDiodeZ+PCiiMZ/2+PCiiBZ;//se puede cambiar zdiode tambien xd
            
            G4VPhysicalVolume *physPCiiM = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logicPCiiM,
                "physPCiiM",
                logicWorld,
                false,
                nmosfet
                
            );
           nmosfet++;
        }
}
}
for(G4int layer = 0; layer < numLayers; layer++) {
    for(G4int i = 0; i < numDetectors; i++) {
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =+0.095*cm+ shiftmosfeX   + i*separation/2;
            G4double yPos =0.019*cm -0.23*cm+0.115*cm -0.5*(numDetectors-1)*separation + j*separation-PCiiBY/2+PCiiMY/2;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) - SDiodeZ+PCiiMZ/2+PCiiBZ;//se puede cambiar zdiode tambien xd
            
            G4VPhysicalVolume *physPCiiM = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logicPCiiM,
                "physPCiiM",
                logicWorld,
                false,
                nmosfet
                
            );
           nmosfet++;
        }
}
}

 G4double PCiiTX =0.043  *cm;
 G4double PCiiTY = 0.044*cm;
 G4double PCiiTZ = 0.014*cm;
  
 
 G4int count_PCiiT=0;
// Crear el volumen sólido del detector
  G4Box * PCiiT = new G4Box("PCsT", PCiiTX/2, PCiiTY/2, PCiiTZ/2);

// Crear el volumen lógico del detector
 G4LogicalVolume *  logicPCiiT = new G4LogicalVolume(PCiiT,cobre, "logicPCsT");
  
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) {
    for(G4int i = 0; i < numDetectors; i++) {
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =-0.095*cm+ shiftmosfeX  + i*separation/2;
            G4double yPos = -0.5*(numDetectors-1)*separation + j*separation-PCiiTY/2 -0.06*cm;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) - SDiodeZ+PCiiBZ/2+PCiiMZ;//se puede cambiar zdiode tambien xd
            
            G4VPhysicalVolume *physPCiiB = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logicPCiiB,
                "physPCiiB",
                logicWorld,
                false,
                nmosfet
                
            );
           nmosfet++;
        }
}
}
for(G4int layer = 0; layer < numLayers; layer++) {
    for(G4int i = 0; i < numDetectors; i++) {
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos =+0.095*cm+ shiftmosfeX  + i*separation/2;
            G4double yPos = -0.5*(numDetectors-1)*separation + j*separation-PCiiTY/2 -0.06*cm;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) - SDiodeZ+PCiiBZ/2+PCiiMZ;//se puede cambiar zdiode tambien xd
            
            G4VPhysicalVolume *physPCiiB = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logicPCiiB,
                "physPCiiB",
                logicWorld,
                false,
                nmosfet
                
            );
           nmosfet++;
        }
}
}















// //!!!!
// //????????????
// //*  MOSFET
//  //!!!!
// //????????????
//*-mosfet-

 G4double mosfetX = 0.28*cm;
 G4double mosfetY = 0.12*cm;
 G4double mosfetZ = 0.105*cm;


//*-INTERNOS-
 G4double PATASUPCOBREMOSFETX = 0.50*mm;
 G4double PATASUPCOBREMOSFETY = 0.35*mm;
 G4double PATASUPCOBREMOSFETZ = 0.20*mm;
 

 
// Crear el volumen sólido del detector
  G4Box * solidmosfet = new G4Box("solidmosfet", mosfetX/2, mosfetY/2, mosfetZ/2);

// Crear el volumen lógico del detector
  logicDetector = new G4LogicalVolume(solidmosfet,MoldCompoundBlack, "logicmosfet");
    
// Colocar los detectores en la cuadrícula y capas
for(G4int layer = 0; layer < numLayers; layer++) {
    for(G4int i = 0; i < numDetectors; i++) {
        for(G4int j = 0; j < numDetectors; j++) {
            G4double xPos = shiftmosfeX  + i*separation/2;
            G4double yPos =  -0.5*(numDetectors-1)*separation + j*separation;
            G4double zPos = -profundidadGrupo -layer * (layerSeparation) - SDiodeZ+ mosfetZ/2 ;
            
            G4VPhysicalVolume *physmosfet = new G4PVPlacement(
                0,
                G4ThreeVector(xPos, yPos, zPos),
                logicDetector,
                "physmosfet",
                logicWorld,
                false,
                 nmosfet
                
            );
             nmosfet++;
        }
}
}
  


//TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//* &&&&&&&&&& CAJA PMMA &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


// lateral
G4double zLateralBoxPMMA = 143 * mm+10*mm;
G4double xLateralBoxPMMA = 650 * mm;
G4double yLateralBoxPMMA = 10*mm; 
G4double xLateralBoxPMMAPos = xLateralBoxPMMA/2-15*cm ;
G4double yLateralBoxPMMAPos = PhantomY/2+yLateralBoxPMMA/2;
G4double zLateralBoxPMMAPos = -zLateralBoxPMMA/2;
 
 
// lateralbox1
G4Box* lateralBoxPMMA= new G4Box("lateralBoxPMMA", xLateralBoxPMMA/2, yLateralBoxPMMA/2,  zLateralBoxPMMA/2);
G4LogicalVolume* logicLateralBoxPMMA  = new G4LogicalVolume(lateralBoxPMMA, PMMA, "logicLateralBoxPMMA");
G4VPhysicalVolume *physLateralBoxPMMA = new G4PVPlacement(0, G4ThreeVector(xLateralBoxPMMAPos, yLateralBoxPMMAPos, zLateralBoxPMMAPos), logicLateralBoxPMMA, "logicLateralBoxPMMA", logicWorld,false, 0, true);




// Lateralbox 2
G4double yLateralBoxPMMAPos2 = -PhantomY/2-yLateralBoxPMMA/2;
G4VPhysicalVolume *physLateralBoxPMMA2 = new G4PVPlacement(0, G4ThreeVector(xLateralBoxPMMAPos, yLateralBoxPMMAPos2, zLateralBoxPMMAPos), logicLateralBoxPMMA, "logicLateralBoxPMMA2", logicWorld,false, 0, true);




//Frontal

G4double zFrontalBoxPMMA = 143 * mm+10*mm;
G4double xFrontalBoxPMMA = 10* mm;
G4double yFrontalBoxPMMA = 280 *mm; 
G4double xFrontalBoxPMMAPos =  -14*cm-1*xFrontalBoxPMMA/2 ;
G4double yFrontalBoxPMMAPos =0 ;
G4double zFrontalBoxPMMAPos = -zFrontalBoxPMMA/2;

G4Box* frontalBoxPMMA= new G4Box("frontalBoxPMMA", xFrontalBoxPMMA/2, yFrontalBoxPMMA/2,  zFrontalBoxPMMA/2);
G4LogicalVolume* logicFrontalBoxPMMA  = new G4LogicalVolume(frontalBoxPMMA, PMMA, "logicFrontalBoxPMMA");
G4VPhysicalVolume *physFrontalBoxPMMA = new G4PVPlacement(0, G4ThreeVector(xFrontalBoxPMMAPos, yFrontalBoxPMMAPos, zFrontalBoxPMMAPos), logicFrontalBoxPMMA, "logicFrontalBoxPMMA", logicWorld,false, 0, true);


G4double xFrontalBoxPMMAPos2 =  -14*cm-3*xFrontalBoxPMMA/2+xLateralBoxPMMA ;
G4VPhysicalVolume *physFrontalBoxPMMA2 = new G4PVPlacement(0, G4ThreeVector(xFrontalBoxPMMAPos2, yFrontalBoxPMMAPos, zFrontalBoxPMMAPos), logicFrontalBoxPMMA, "logicFrontalBoxPMMA2", logicWorld,false, 0, true);




//Frontal

G4double zBaseBoxPMMA = 10 * mm;
G4double xBaseBoxPMMA = 630* mm;
G4double yBaseBoxPMMA = 280* mm; 
G4double xBaseBoxPMMAPos =  -14*cm+ xBaseBoxPMMA/2;
G4double yBaseBoxPMMAPos =0 ;
G4double zBaseBoxPMMAPos = - PhantomZ-zBaseBoxPMMA/2;

G4Box* baseBoxPMMA= new G4Box("baseBoxPMMA", xBaseBoxPMMA/2, yBaseBoxPMMA/2,  zBaseBoxPMMA/2);
G4LogicalVolume* logicBaseBoxPMMA  = new G4LogicalVolume(baseBoxPMMA, PMMA, "logicBaseBoxPMMA");
G4VPhysicalVolume *physBaseBoxPMMA = new G4PVPlacement(0, G4ThreeVector(xBaseBoxPMMAPos, yBaseBoxPMMAPos, zBaseBoxPMMAPos), logicBaseBoxPMMA, "logicBaseBoxPMMA", logicWorld,false, 0, true);




G4Material* CuCu = nist->FindOrBuildMaterial("G4_Cu");


//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//* &&&&&&&&&& COBRE CAJA FARADAY &&&&&&&&&&&&&&&&&&&&&&&&
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

// Parámetros configurables
G4double thickness = 2.0 * mm;
G4double width = 10.0 * cm;
G4double height = BakeliteHeight; 
G4double posX = 9.0 * cm+width/2 ;
G4double posY = 0.0 * cm;
G4double posZ = PhantomZ/2-thickness/2-1*cm;

 
// Lámina de aluminio
G4Box* solidPlate = new G4Box("AlPlate", width/2, height/2, thickness/2);
G4LogicalVolume* logicPlate = new G4LogicalVolume(solidPlate, CuCu, "AlPlate");
G4VPhysicalVolume *physAlPlate = new G4PVPlacement(0, G4ThreeVector(posX, posY, posZ), logicPlate, "AlPlate", logiPhantom, false, 0, true);









// Parámetros configurables
G4double Cuthickness = 2 * mm;
G4double Cuheight = 19.00 * cm;
G4double Cuwidth = Cuheight / 2;
G4double CuposX = 8.8 * cm + 5 * cm + Cuwidth / 2;
G4double CuposY = 0.0 * cm;
G4double initialCuposZ = -1.49 * cm + mosfetZ / 2;
G4double CUseparation = 2.0 * cm;



 
  initialCuposZ = -1.39 * cm + mosfetZ / 2;
  Cuthickness = 2 * mm;

// Definir la lámina de cobre
G4Box* CusolidPlate2 = new G4Box("CuPlate2", Cuwidth / 2, Cuheight / 2, Cuthickness / 2);
G4LogicalVolume* CulogicPlate2 = new G4LogicalVolume(CusolidPlate2, CuCu, "CuPlate2");

// Colocar 5 placas con separación de 2 cm en Z
for (G4int i = 0; i < 5; i++) {
    G4double CuposZ = initialCuposZ - i * CUseparation;
    G4VPhysicalVolume* cuPlate2PhysicalVolume = new G4PVPlacement(0, G4ThreeVector(CuposX, CuposY, CuposZ), CulogicPlate2, "CuPlate2", logicWorld, false, i, true);
      
}




// Parámetros para CuPlate3
G4double CuSizeZ_CuPlate3 = 3* mm; // Tamaño en Z (grosor) para CuPlate3
G4double CuSizeY_CuPlate3 =  19* cm; // Tamaño en Y (altura) para CuPlate3 (reusando valor de CuPlate2)
G4double CuSizeX_CuPlate3 = 5*0.15 *mm; // Tamaño en X (ancho) para CuPlate3 (calculado)
G4double CuposX_CuPlate3 = 8.8 * cm + 5 * cm -CuSizeX_CuPlate3/2; // Posición X para CuPlate3
G4double CuposY_CuPlate3 = 0.0 * cm; // Posición Y para CuPlate3
G4double initialCuposZ_CuPlate3 = -1.5 * cm + mosfetZ ; // Posición inicial Z para CuPlate3
G4double CUseparation_CuPlate3 = 2.0 * cm; // Separación en Z para CuPlate3 (reusando valor)


// Definir la lámina de cobre para CuPlate3
G4Box* CusolidPlate3 = new G4Box("CuPlate3", CuSizeX_CuPlate3 / 2, CuSizeY_CuPlate3 / 2, CuSizeZ_CuPlate3 / 2);
G4LogicalVolume* CulogicPlate3 = new G4LogicalVolume(CusolidPlate3, CuCu, "CuPlate3"); // Usando material de cobre

// Colocar 5 placas de CuPlate3 con separación en Z
for (G4int i = 0; i < 5; i++) {
    G4double CuposZ_CuPlate3_Loop = initialCuposZ_CuPlate3 - i * CUseparation_CuPlate3;
    G4VPhysicalVolume* cuPlate3PhysicalVolume = new G4PVPlacement(0, G4ThreeVector(CuposX_CuPlate3, CuposY_CuPlate3, CuposZ_CuPlate3_Loop), CulogicPlate3, "CuPlate3_Physical", logicWorld, false, i, true);

}
























// //*
// //*
// //!COLORES
// //*
// //*
 
// G4VisAttributes* visAttributesSDIODE = new G4VisAttributes(G4Colour(0.1, 0.1, 0.1, 1.0));  // RGB + Opacidad

// // Hacer que la geometría sea sólida (opaca)
// visAttributesSDIODE->SetVisibility(true);   // Hacer visible
// visAttributesSDIODE->SetForceSolid(true);   // Hacer sólida (no transparente)


//     logicDetector->SetVisAttributes(visAttributesSDIODE);


 
 
//  // Definir el color cobrizo (usando valores RGB)
// G4VisAttributes* copperColor = new G4VisAttributes(G4Colour(0.72, 0.45, 0.20)); // Color cobre (R=0.72, G=0.45, B=0.20)
// copperColor->SetVisibility(true); // Asegurarse de que el volumen sea visible
// copperColor->SetForceSolid(true); // Mostrar el objeto sólido, no solo como un contorno

// // // Asignar el color al volumen lógico de las conexiones de cobre
//      logiccobreConectionDetector ->SetVisAttributes(copperColor); 
//     logiccobreConectioncapacitor->SetVisAttributes(copperColor);
// logiccobrePista1->SetVisAttributes(copperColor);

// G4VisAttributes* visAttributesBakelite = new G4VisAttributes(G4Colour(0.1, 0.5, 0.1, 1.0));  // Verde oscuro, sin transparencia

// // Hacer visible y sólido el volumen de Bakelita
// visAttributesBakelite->SetVisibility(true);    // Hacer visible
// visAttributesBakelite->SetForceSolid(true);    // Hacer sólido (no transparente)

// // Asignar los atributos visuales al volumen lógico de la Bakelita
//     logicBakelite->SetVisAttributes(visAttributesBakelite); 
//     // logicCOBREMOSFET->SetVisAttributes(copperColor);

 

// G4VisAttributes* visAttributesmosfet= new G4VisAttributes(G4Colour(0.3, 0.2, 0.5, 1.0));  // RGB + Opacidad
// // Hacer que la geometría sea sólida (opaca)
// visAttributesmosfet->SetVisibility(true);   // Hacer visible
// visAttributesmosfet->SetForceSolid(true);   // Hacer sólida (no transparente)
 
//     logicmosfet->SetVisAttributes(visAttributesmosfet);
 

// G4VisAttributes* visAttributesCapacitor = new G4VisAttributes(G4Colour(1.0, 0.8, 0.2, 1.0));  // Naranja amarillento, sin transparencia

// // Hacer visible y sólido el capacitor
// visAttributesCapacitor->SetVisibility(true);   // Hacer visible
// visAttributesCapacitor->SetForceSolid(true);   // Hacer sólido (no transparente)

// // Asignar los atributos visuales al volumen lógico del capacitor
//     logiccapacitor->SetVisAttributes(visAttributesCapacitor);




// // Crear los atributos de visualización para el color plata claro
// G4VisAttributes* silverLightVisAttr = new G4VisAttributes(G4Colour(0.75, 0.75, 0.75, 1.0)); // Plata claro, opaco
// silverLightVisAttr->SetVisibility(true);   // Hacerlo visible
// silverLightVisAttr->SetForceSolid(true);   // Forzar la visualización en sólido

// // Asignar los atributos de visualización al volumen lógico logiccobrePatitaBaseDetector2
// logiccobrePatitaBaseDetector2->SetVisAttributes(silverLightVisAttr);

// // Asignar los atributos de visualización al volumen lógico logiccobrePatitaBaseDetector
// logiccobrePatitaBaseDetector->SetVisAttributes(silverLightVisAttr);


// // Crear los atributos de visualización para la cobertura del diodo
// G4VisAttributes* coberturaDiodeVisAttr = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0, 0.2)); // Blanco con 50% de transparencia
// coberturaDiodeVisAttr->SetVisibility(true);    // Hacerlo visible
// coberturaDiodeVisAttr->SetForceSolid(true);    // Forzar la visualización en sólido

// // Asignar los atributos de visualización al volumen lógico de la cobertura
// logicCoberturaDiodeConHueco->SetVisAttributes(coberturaDiodeVisAttr);

// // PATASUPlogicCOBREMOSFET->SetVisAttributes(silverLightVisAttr);
// logicPCsB ->SetVisAttributes(silverLightVisAttr);
// logicPCsM ->SetVisAttributes(silverLightVisAttr);
// logicPCsT ->SetVisAttributes(silverLightVisAttr);
// logicPCiiT ->SetVisAttributes(silverLightVisAttr);
// logicPCiiB ->SetVisAttributes(silverLightVisAttr);
// logicPCiiM ->SetVisAttributes(silverLightVisAttr);

// // Crear los atributos de visualización para los volúmenes

// // Atributos para logicWorld - color azul claro transparente
// G4VisAttributes* logicWorldVisAttr = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.1 )); // Azul claro, 50% transparencia
// logicWorldVisAttr->SetVisibility(true);    // Hacerlo visible
// logicWorldVisAttr->SetForceSolid(true);    // Forzar la visualización en sólido

// // Atributos para logiPhantom - color verde claro transparente
// G4VisAttributes* logiPhantomVisAttr = new G4VisAttributes(G4Colour(1.0, 1.0, 0.9, 0.1)); // Verde claro, 70% transparencia
// logiPhantomVisAttr->SetVisibility(true);   // Hacerlo visible
// logiPhantomVisAttr->SetForceSolid(true);   // Forzar la visualización en sólido

// // Atributos para medioVacioLogic - color rojo claro transparente
// // G4VisAttributes* medioVacioLogicVisAttr = new G4VisAttributes(G4Colour(0.1, 0.1, 0.3, 0.)); // Rojo claro, 80% transparencia
// // medioVacioLogicVisAttr->SetVisibility(true); // Hacerlo visible
// // medioVacioLogicVisAttr->SetForceSolid(true); // Forzar la visualización en sólido

// // Asignar los atributos a los volúmenes lógicos
// logicWorld->SetVisAttributes(logicWorldVisAttr);
// logiPhantom->SetVisAttributes(logiPhantomVisAttr);
// // medioVacioLogic->SetVisAttributes(medioVacioLogicVisAttr);
//  logicCOBREMOSFET->SetVisAttributes(copperColor); 


 





 





 
    return physWorld;
}

void MyDetectorConstruction::ConstructSDandField() {
    auto sensDet = new MySensitiveDetector("SensitiveDetector", "TrackerHitsCollection");
    G4SDManager::GetSDMpointer()->AddNewDetector(sensDet);

    if (logicmosfet != nullptr) {
        logicmosfet->SetSensitiveDetector(sensDet);
    }
}
