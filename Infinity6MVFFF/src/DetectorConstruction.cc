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

    G4Material* matCu     = nist->FindOrBuildMaterial("G4_Cu");
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

   G4double fractionmass;
    G4Element* elPb = nist->FindOrBuildElement("Pb");
    G4Element* elSb = nist->FindOrBuildElement("Sb");
    G4Material* matPbSb = new G4Material("PbSb_Alloy", 11.0*g/cm3, 2);
    matPbSb->AddElement(elPb, fractionmass=0.96);
    matPbSb->AddElement(elSb, fractionmass=0.04);

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
G4double target_radius    = 10.0 * mm;   // Radio del disco
G4double target_thickness =0.09*cm; //1.1 * cm (no sirvio, la pdd bajo demasiado) //0.09cm va bien, podria mejoraren FFF todo se malogra ;  // Grosor en dirección del haz
G4double SSDValue         =100 * cm;

G4Tubs* solidTarget = new G4Tubs("Target",
    0., target_radius,
    target_thickness / 2.,
    0. * deg, 360. * deg);

G4LogicalVolume* logicTarget = new G4LogicalVolume(solidTarget, wre_alloy, "TargetLV");

G4Region* TargetRegion = new G4Region("TargetRegion");
TargetRegion->AddRootLogicalVolume(logicTarget);

G4ThreeVector posTarget(0, 0, SSDValue - 1.1*cm);
new G4PVPlacement(nullptr, posTarget, logicTarget, "Target", logicWorld, false, 0, true);

G4VisAttributes* targetVis = new G4VisAttributes(G4Colour::Gray());
targetVis->SetVisibility(true);
logicTarget->SetVisAttributes(targetVis); 
 

 {
        G4double r = target_radius;
        G4double h = 1.016*cm;

        G4Tubs* solid14 =
            new G4Tubs("TARGET_14",0,r,h/2,0,360*deg);

        G4LogicalVolume* logic14 =
            new G4LogicalVolume(solid14, matCu, "TARGET_14");

        new G4PVPlacement(
            0, G4ThreeVector(0,0, SSDValue - 1.1*cm - target_thickness/2  - h/2),
            logic14,"TARGET_14",logicWorld,false,0,true);

    
    

    G4double pb_radius = 30.15*mm;
    G4double pb_height = 1.0*mm;

    G4Tubs* solidPlatePbSb =
        new G4Tubs("PlatePbSb",0,pb_radius,pb_height/2,0,360*deg);

    G4LogicalVolume* logicPlatePbSb =
        new G4LogicalVolume(solidPlatePbSb, matPbSb,"PlatePbSb");

    new G4PVPlacement(
        0,
        G4ThreeVector(0,0,SSDValue-104.648*mm),
        logicPlatePbSb,"PlatePbSb",
        logicWorld,false,0,true);


    }

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
G4double fieldSize  = 15.0 * cm;      // <-- Cambiar campo (ej. 10, 20, etc.)
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
    G4double PhantomX = 30*cm;
    G4double PhantomY = 30*cm;
    G4double PhantomZ =30*cm;
    G4Material* phantomMaterial =nist->FindOrBuildMaterial("G4_WATER");
 
 
// Crear el volumen del fantoma principal de PMMA
G4Box* solidPhantom = new G4Box("solidPhantom", PhantomX/2, PhantomY/2, PhantomZ/2);
G4VSolid* solidPhantomWithHoles = solidPhantom; // Inicialmente, es el bloque completo

  
G4LogicalVolume* logicPhantom = new G4LogicalVolume(solidPhantom, phantomMaterial, "logicalPhantom");

// Colocar el fantoma en el mundo
G4VPhysicalVolume* physPhantom = new G4PVPlacement(0, G4ThreeVector(0, 0, -PhantomZ/2),  
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
G4int numDetectors =  1; // 23;  

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
  logicmosfet = new G4LogicalVolume(solidDetector,phantomMaterial, "logicDetector");
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
