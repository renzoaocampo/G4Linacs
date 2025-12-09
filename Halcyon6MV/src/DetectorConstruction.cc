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
#include "G4UnionSolid.hh"
#include "SensitiveDetector.hh"

MyDetectorConstruction::MyDetectorConstruction() {}
MyDetectorConstruction::~MyDetectorConstruction() {}
 
G4VPhysicalVolume* MyDetectorConstruction::Construct() {
    G4NistManager* nist = G4NistManager::Instance();
G4Material* tungsten = nist->FindOrBuildMaterial("G4_W");
G4Material* copper = nist->FindOrBuildMaterial("G4_Cu");
G4Material* iron = nist->FindOrBuildMaterial("G4_Fe");
G4Material* air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");


//*=========Material ========= 
    // Materiales
    G4Material* worldMat = nist->FindOrBuildMaterial("G4_AIR");
    G4Material* detectorMat = nist->FindOrBuildMaterial("G4_WATER"); 
G4Material* stainlessSteel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
G4Material* aluminum = nist->FindOrBuildMaterial("G4_Al"); 

// Alternativa: aleación W-Re (90% W, 10% Re)
G4Element* elRe = new G4Element("Rhenium", "Re", 75., 186.207*g/mole);
G4Material* wre_alloy = new G4Material("WRe_Alloy", 19.25*g/cm3, 2);
wre_alloy->AddMaterial(tungsten, 90.*perCent);
wre_alloy->AddElement(elRe, 10.*perCent);
G4double density = 18.0 * g/cm3;
G4Material* matTungstenMLC = new G4Material("MLCTungstenAlloy", density, 3);
matTungstenMLC->AddElement(G4NistManager::Instance()->FindOrBuildElement("W"), 95.0 * perCent);
matTungstenMLC->AddElement(G4NistManager::Instance()->FindOrBuildElement("Ni"), 3.75 * perCent);
matTungstenMLC->AddElement(G4NistManager::Instance()->FindOrBuildElement("Fe"), 1.25 * perCent);
G4bool checkOverlaps = true;


G4Material* revestimientoMaterial = nist->FindOrBuildMaterial("G4_W");
//*=========Mundo=========
    G4double worldSize = 5.0 * m;
    solidWorld = new G4Box("World", worldSize / 2, worldSize / 2, worldSize / 2);
    logicWorld = new G4LogicalVolume(solidWorld, worldMat, "World");
    physWorld = new G4PVPlacement(
        0, G4ThreeVector(), 
        logicWorld, "World", 0, false, 0, true);

G4LogicalVolume* motherVolume = logicWorld;

G4LogicalVolume* motherLogical = logicWorld;













    // ==========================================================
    // MATERIALES
    // ==========================================================
    G4Material* matAir    = nist->FindOrBuildMaterial("G4_AIR");
    G4Material* matWater  = nist->FindOrBuildMaterial("G4_WATER");
    G4Material* matW      = nist->FindOrBuildMaterial("G4_W");
    G4Material* matFe     = nist->FindOrBuildMaterial("G4_Fe");
    G4Material* matCarbon = nist->FindOrBuildMaterial("G4_C");
    G4Material* matCu     = nist->FindOrBuildMaterial("G4_Cu");

    // Pb/Sb para la placa
    G4double fractionmass;
    G4Element* elPb = nist->FindOrBuildElement("Pb");
    G4Element* elSb = nist->FindOrBuildElement("Sb");
    G4Material* matPbSb = new G4Material("PbSb_Alloy", 11.0*g/cm3, 2);
    matPbSb->AddElement(elPb, fractionmass=0.96);
    matPbSb->AddElement(elSb, fractionmass=0.04);

   
    // ==========================================================
    // MESA
    // ==========================================================
    G4Box* solidCama = new G4Box("Cama", 27*cm, 90*cm, 2.5*cm);
    G4LogicalVolume* logicCama =
        new G4LogicalVolume(solidCama, matCarbon, "Cama");

    new G4PVPlacement(
        0, G4ThreeVector(0,20*cm,132.5*cm),
        logicCama,"Cama",logicWorld,false,0,true);

    G4Box* solidBase =
        new G4Box("BaseMesa", 25*cm, 65*cm, 42.75*cm);
    G4LogicalVolume* logicBase =
        new G4LogicalVolume(solidBase, matFe, "BaseMesa");

    new G4PVPlacement(
        0, G4ThreeVector(0,135*cm,178.25*cm),
        logicBase,"BaseMesa",logicWorld,false,0,true);


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


    // ==========================================================
    // PLACA Pb/Sb — MCNP RCC superficie 92
    // ==========================================================
    /*
        c Placa Pb/Sb
        903 53 -11.34 92
        92 RCC 0 0 104.648   0 0 1.0   30.15
          → cilindro: altura = 1 mm, radio = 30.15 mm
    */

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


    // --- Superficie 15 ---
    /*
        15 RCC 0 0 -0.254   0 0 1.524   0.889
    */
    {
        G4double r = 0.889*cm;
        G4double h = 1.524*cm;

        G4Tubs* solid15 =
            new G4Tubs("TARGET_15",0,r,h/2,0,360*deg);

        G4LogicalVolume* logic15 =
            new G4LogicalVolume(solid15, matCu, "TARGET_15");

        new G4PVPlacement(
            0, G4ThreeVector(0,0,-0.254*cm + h/2),
            logic15,"TARGET_15",logicWorld,false,0,true);
    }


    // --- Superficie 16 ---
    /*
        16 RCC 0 0 -0.254   0 0 1.524   0.301
    */
    {
        G4double r = 0.301*cm;
        G4double h = 1.524*cm;

        G4Tubs* solid16 =
            new G4Tubs("TARGET_16",0,r,h/2,0,360*deg);

        G4LogicalVolume* logic16 =
            new G4LogicalVolume(solid16, matCu, "TARGET_16");

        new G4PVPlacement(
            0, G4ThreeVector(0,0,-0.254*cm + h/2),
            logic16,"TARGET_16",logicWorld,false,0,true);
    }


    // ==========================================================
    // VISUALIZACIÓN
    // ========================================================== 
    logicPrimColl->SetVisAttributes(new G4VisAttributes(G4Color(1,0,0)));
    logicPlatePbSb->SetVisAttributes(new G4VisAttributes(G4Color(0.3,0.3,0.3)));












 
    
G4RotationMatrix* rotXJawMas = new G4RotationMatrix();
rotXJawMas->rotateY(2.8624*deg);
G4RotationMatrix* rotXJawMenos = new G4RotationMatrix();
rotXJawMenos->rotateY(-2.8624*deg);
// -------------------- X JAWS --------------------
// half-lengths y centro calculados de las RPP 69/70
G4double jawX_halfX = 11.938*cm/2;
G4double jawX_halfY = 10.795*cm;
G4double jawX_halfZ = 3.9*cm;
G4double jawX_centerZ = (36.703+44.503)/2*cm;
G4double jawX_shift=8.0066*cm;
G4Box* jawXBox = new G4Box("JawXBox", jawX_halfX, jawX_halfY, jawX_halfZ);

// Jaw X izquierda (RPP 69): centerX = -5.969 cm
G4LogicalVolume* jawXnegLog = new G4LogicalVolume(jawXBox, tungsten, "JawXneg_log");
new G4PVPlacement(rotXJawMas,
                  G4ThreeVector(- jawX_shift, 0.0, jawX_centerZ),
                  jawXnegLog,
                  "JawXneg_phys",
                  motherVolume,
                  false,
                  0,
                  true);

// Jaw X derecha (RPP 70): centerX = +5.969 cm
G4LogicalVolume* jawXposLog = new G4LogicalVolume(jawXBox, tungsten, "JawXpos_log");
new G4PVPlacement(rotXJawMenos,
                  G4ThreeVector(+ jawX_shift, 0.0, jawX_centerZ),
                  jawXposLog,
                  "JawXpos_phys",
                  motherVolume,
                  false,
                  0,
                  true);




                  
    
G4RotationMatrix* rotYJawMas = new G4RotationMatrix();
rotYJawMas->rotateX(2.8624*deg);
G4RotationMatrix* rotYJawMenos = new G4RotationMatrix();
rotYJawMenos->rotateX(-2.8624*deg);
// -------------------- Y JAWS --------------------
// half-lengths y centro calculados de las RPP 67/68
G4double jawY_halfX = 9.398*cm;
G4double jawY_halfY =11.938/2*cm;
G4double jawY_halfZ = 3.9*cm;
G4double jawY_centerZ =  (28.067+35.867)/2*cm;
G4double jawY_shift=7.5748*cm;

G4Box* jawYBox = new G4Box("JawYBox", jawY_halfX, jawY_halfY, jawY_halfZ);

// Jaw Y negativa (RPP 67): centerY = -5.969 cm
G4LogicalVolume* jawYnegLog = new G4LogicalVolume(jawYBox, tungsten, "JawYneg_log");
new G4PVPlacement(rotYJawMenos,
                  G4ThreeVector(0.0, -jawY_shift, jawY_centerZ),
                  jawYnegLog,
                  "JawYneg_phys",
                  motherVolume,
                  false,
                  0,
                  true);

// Jaw Y positiva (RPP 68): centerY = +5.969 cm
G4LogicalVolume* jawYposLog = new G4LogicalVolume(jawYBox, tungsten, "JawYpos_log");
new G4PVPlacement(rotYJawMas,
                  G4ThreeVector(0.0, +jawY_shift, jawY_centerZ),
                  jawYposLog,
                  "JawYpos_phys",
                  motherVolume,
                  false,
                  0,
                  true);

// Visual
G4VisAttributes* visJaw = new G4VisAttributes(G4Colour(0.2,0.2,0.2));
visJaw->SetForceSolid(true);


jawYnegLog->SetVisAttributes(visJaw);
jawYposLog->SetVisAttributes(visJaw);

jawXnegLog->SetVisAttributes(visJaw);
jawXposLog->SetVisAttributes(visJaw);

    
    
    
     

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
G4int numDetectors =  61; // 23;  
  
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

    if (logicVoxel != nullptr) {
       logicVoxel->SetSensitiveDetector(sensDet);
    }
}
