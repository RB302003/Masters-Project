/**
 * @file MyDetectorConstruction.cc
 * @brief Defines the physical geometry and materials of the simulation.
 * 
 * WHAT THIS CLASS DOES:
 * This class is the architect of the simulation. It builds the detector universe:
 * 1. Defines the materials (e.g., Lead Tungstate PbWO4, Air, Vacuum).
 * 2. Assigns OPTICAL properties to these materials (refractive index, scintillation 
 *    yield, absorption length) so photons behave realistically.
 * 3. Constructs the actual 3D shapes (a 7x7x24 grid of PbWO4 crystals).
 * 4. Applies a reflective skin surface (wrapping) to the crystals to bounce 
 *    photons back inside, simulating real-world detector wrapping.
 */

#include "MyDetectorConstruction.hh"
#include "MySensitiveDetector.hh"
#include "G4SDManager.hh"

MyDetectorConstruction::MyDetectorConstruction()
  : logicDetector(nullptr),
    logicAbsorber(nullptr),
    fMessenger(nullptr),
    fMode("homogeneous")
{
  nCols = 7;
  nRows = 7;

  // Homogeneous mode parameters
  nLayers = 24;

  // Sampling mode parameters
  absorberLayers = 20;
  activeLayers   = 20;

  DefineMaterial();

  xWorld = 2.0 * m;
  yWorld = 2.05 * m;
  zWorld = 4.0 * m;

  // Create messenger for runtime mode selection
  fMessenger = new G4GenericMessenger(this, "/detector/", "Detector configuration");
  fMessenger->DeclareProperty("mode", fMode,
      "Calorimeter mode: 'homogeneous' (PbWO4 24-layer) or 'sampling' (Pb+Scint 20-layer)")
      .SetDefaultValue("homogeneous");
}

MyDetectorConstruction::~MyDetectorConstruction() {
  delete fMessenger;
}

void MyDetectorConstruction::DefineMaterial() {
  G4NistManager* nist = G4NistManager::Instance();

  Aluminium = nist->FindOrBuildMaterial("G4_Al");
  Air = nist->FindOrBuildMaterial("G4_AIR");

  PbWO4 = nist->FindOrBuildMaterial("G4_PbWO4");

  Lead = nist->FindOrBuildMaterial("G4_Pb");
  Scintillator = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

  // Define optical properties for materials
  DefinePbWO4OpticalProperties();
  DefineAirOpticalProperties();
}

void MyDetectorConstruction::DefinePbWO4OpticalProperties() {

  G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

  // Energy points covering 300 nm – 700 nm
  // E = hc/lambda, where hc = 1239.84 eV·nm
  const G4int nEntries = 8;
  G4double photonEnergy[nEntries] = {
    1.771*eV,  // 700 nm
    2.066*eV,  // 600 nm
    2.480*eV,  // 500 nm
    2.755*eV,  // 450 nm
    2.952*eV,  // 420 nm  (scintillation peak)
    3.100*eV,  // 400 nm
    3.542*eV,  // 350 nm
    4.133*eV   // 300 nm
  };

  G4double rindex[nEntries] = {
    2.160, 2.170, 2.190, 2.200, 2.210, 2.220, 2.250, 2.300
  };
  mpt->AddProperty("RINDEX", photonEnergy, rindex, nEntries);

  G4double absLength[nEntries] = {
    200.0*cm,  // 700 nm - very transparent
    180.0*cm,  // 600 nm
    150.0*cm,  // 500 nm
    120.0*cm,  // 450 nm
    100.0*cm,  // 420 nm - around emission peak
     80.0*cm,  // 400 nm
     40.0*cm,  // 350 nm - UV starts absorbing
     10.0*cm   // 300 nm - strong UV absorption
  };
  mpt->AddProperty("ABSLENGTH", photonEnergy, absLength, nEntries);

  // ------- Scintillation Yield -------
  // PbWO4: ~200 photons/MeV at 18°C (low compared to NaI)
  // Using reduced yield of 10/MeV for practical simulation times
  // Scale up to 200/MeV for production runs
  // CMS ECAL nominal yield: ~100-200 ph/MeV at 18°C
  mpt->AddConstProperty("SCINTILLATIONYIELD", 200.0 / MeV);
  mpt->AddConstProperty("RESOLUTIONSCALE", 1.0);

  G4double scintComp1[nEntries] = {
    0.00,  // 700 nm
    0.02,  // 600 nm
    0.10,  // 500 nm
    0.40,  // 450 nm
    1.00,  // 420 nm  (peak)
    0.60,  // 400 nm
    0.10,  // 350 nm
    0.01   // 300 nm
  };
  mpt->AddProperty("SCINTILLATIONCOMPONENT1", photonEnergy, scintComp1, nEntries);
  mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 10.0*ns);
  mpt->AddConstProperty("SCINTILLATIONYIELD1", 0.80);

  G4double scintComp2[nEntries] = {
    0.01,  // 700 nm
    0.10,  // 600 nm
    1.00,  // 500 nm  (peak)
    0.50,  // 450 nm
    0.20,  // 420 nm
    0.05,  // 400 nm
    0.01,  // 350 nm
    0.00   // 300 nm
  };
  mpt->AddProperty("SCINTILLATIONCOMPONENT2", photonEnergy, scintComp2, nEntries);
  mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 30.0*ns);
  mpt->AddConstProperty("SCINTILLATIONYIELD2", 0.20);

  PbWO4->SetMaterialPropertiesTable(mpt);

  G4cout << "========================================" << G4endl;
  G4cout << "  PbWO4 optical properties defined" << G4endl;
  G4cout << "  Scintillation yield: 10 photons/MeV" << G4endl;
  G4cout << "  Refractive index: 2.16 - 2.30" << G4endl;
  G4cout << "========================================" << G4endl;
}

void MyDetectorConstruction::DefineAirOpticalProperties() {

  G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

  const G4int nEntries = 2;
  G4double photonEnergy[nEntries] = { 1.771*eV, 4.133*eV };
  G4double rindex[nEntries]      = { 1.000293, 1.000293 };

  mpt->AddProperty("RINDEX", photonEnergy, rindex, nEntries);

  Air->SetMaterialPropertiesTable(mpt);
}

void MyDetectorConstruction::DefineCrystalSurface() {
  if (!logicDetector) return;

  // Create a polished reflective wrapping (simulating Tyvek/ESR foil)
  G4OpticalSurface* crystalSurface = new G4OpticalSurface("CrystalWrapping");
  crystalSurface->SetModel(unified);
  crystalSurface->SetFinish(polished);
  crystalSurface->SetType(dielectric_metal);

  // Surface properties: 95% reflectivity across all wavelengths
  const G4int nEntries = 2;
  G4double photonEnergy[nEntries] = { 1.771*eV, 4.133*eV };
  G4double reflectivity[nEntries] = { 0.95, 0.95 };
  G4double efficiency[nEntries]   = { 0.0,  0.0 };  // no photon detection at surface

  G4MaterialPropertiesTable* surfMPT = new G4MaterialPropertiesTable();
  surfMPT->AddProperty("REFLECTIVITY", photonEnergy, reflectivity, nEntries);
  surfMPT->AddProperty("EFFICIENCY",   photonEnergy, efficiency,   nEntries);
  crystalSurface->SetMaterialPropertiesTable(surfMPT);

  // Apply as a skin surface on all PbWO4 crystals
  new G4LogicalSkinSurface("CrystalSkin", logicDetector, crystalSurface);

  G4cout << "  Crystal wrapping applied: 95% reflective skin surface" << G4endl;
}

G4VPhysicalVolume* MyDetectorConstruction::Construct() {
  solidWorld = new G4Box("World", 0.5 * xWorld, 0.5 * yWorld, 0.5 * zWorld);
  logicWorld = new G4LogicalVolume(solidWorld, Air, "World");
  physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", nullptr, false, 0, false);

  G4cout << "============================================" << G4endl;
  G4cout << "  Detector mode: " << fMode << G4endl;
  G4cout << "============================================" << G4endl;

  if (fMode == "sampling") {

    G4double absorberThickness = 5.0 * mm;
    G4double activeThickness   = 2.0 * mm;
    G4double layerThickness    = absorberThickness + activeThickness;

    G4double xsize = 2.2 * cm;
    G4double ysize = 2.2 * cm;

    G4double totalThickness = absorberLayers * layerThickness;

    // Aluminium enclosure sized for sampling geometry
    G4double alHalfX = 0.5 * xsize * nCols + 0.2 * cm;
    G4double alHalfY = 0.5 * ysize * nRows + 0.2 * cm;
    G4double alHalfZ = 0.5 * totalThickness + 0.2 * cm;

    G4VSolid* solidAlBox = new G4Box("AluminiumBox", alHalfX, alHalfY, alHalfZ);
    G4LogicalVolume* logicAlBox = new G4LogicalVolume(solidAlBox, Aluminium, "AluminiumBox");
    new G4PVPlacement(0, G4ThreeVector(), logicAlBox, "AluminiumBox", logicWorld, false, 0, false);

    G4VSolid* solidInside = new G4Box("InsideBox", alHalfX - 0.2*cm, alHalfY - 0.2*cm, alHalfZ - 0.2*cm);
    G4LogicalVolume* logicInside = new G4LogicalVolume(solidInside, Air, "InsideBox");
    new G4PVPlacement(0, G4ThreeVector(), logicInside, "InsideBox", logicAlBox, false, 0, false);

    // Absorber tile (Lead)
    G4VSolid* solidAbsorber = new G4Box("Absorber",
                                        0.5 * xsize,
                                        0.5 * ysize,
                                        0.5 * absorberThickness);
    logicAbsorber = new G4LogicalVolume(solidAbsorber, Lead, "Absorber");

    // Active tile (Scintillator)
    G4VSolid* solidDetector = new G4Box("Detector",
                                        0.5 * xsize,
                                        0.5 * ysize,
                                        0.5 * activeThickness);
    logicDetector = new G4LogicalVolume(solidDetector, Scintillator, "Detector");

    for (G4int layer = 0; layer < absorberLayers; ++layer) {

      G4double zAbs = -0.5 * totalThickness +
                      layer * layerThickness +
                      0.5 * absorberThickness;

      G4double zAct = -0.5 * totalThickness +
                      layer * layerThickness +
                      absorberThickness +
                      0.5 * activeThickness;

      for (G4int col = 0; col < nCols; ++col) {
        G4double ypos = -0.5 * ysize * (nCols - 1) + col * ysize;

        for (G4int row = 0; row < nRows; ++row) {
          G4double xpos = -0.5 * xsize * (nRows - 1) + row * xsize;

          G4int absorberCopyNo = layer * 10000 + col * 100 + row;
          G4int activeCopyNo   = (layer + 100) * 10000 + col * 100 + row;

          // Place absorber
          new G4PVPlacement(
              0,
              G4ThreeVector(xpos, ypos, zAbs),
              logicAbsorber,
              "Absorber",
              logicInside,
              false,
              absorberCopyNo,
              false);

          // Place active detector
          new G4PVPlacement(
              0,
              G4ThreeVector(xpos, ypos, zAct),
              logicDetector,
              "Detector",
              logicInside,
              false,
              activeCopyNo,
              false);
        }
      }
    }

  } else {

    G4double activeThickness = 10.0 * mm;
    G4double xsize = 2.2 * cm;
    G4double ysize = 2.2 * cm;

    G4VSolid* solidAlBox = new G4Box("AluminiumBox", 25.2 * cm, 15.2 * cm, 60.95 * cm);
    G4LogicalVolume* logicAlBox = new G4LogicalVolume(solidAlBox, Aluminium, "AluminiumBox");
    new G4PVPlacement(0, G4ThreeVector(), logicAlBox, "AluminiumBox", logicWorld, false, 0, false);

    G4VSolid* solidInside = new G4Box("InsideBox", 25.0 * cm, 15.0 * cm, 60.75 * cm);
    G4LogicalVolume* logicInside = new G4LogicalVolume(solidInside, Air, "InsideBox");
    new G4PVPlacement(0, G4ThreeVector(), logicInside, "InsideBox", logicAlBox, false, 0, false);

    G4VSolid* solidDetector = new G4Box("Detector",
                                        0.5 * xsize,
                                        0.5 * ysize,
                                        0.5 * activeThickness);

    logicDetector =
        new G4LogicalVolume(solidDetector, PbWO4, "Detector");

    G4double totalThickness = (nLayers * activeThickness);

    for (G4int layer = 0; layer < nLayers; ++layer) {

      G4double zAct = -0.5 * totalThickness +
                      layer * activeThickness +
                      0.5 * activeThickness;

      for (G4int col = 0; col < nCols; ++col) {
        G4double ypos = -0.5 * ysize * (nCols - 1) + col * ysize;

        for (G4int row = 0; row < nRows; ++row) {
          G4double xpos = -0.5 * xsize * (nRows - 1) + row * xsize;

          G4int activeCopyNo   = (layer + 100) * 10000 + col * 100 + row;

          new G4PVPlacement(
              0,
              G4ThreeVector(xpos, ypos, zAct),
              logicDetector,
              "Detector",
              logicInside,
              false,
              activeCopyNo,
              false);
        }
      }
    }

    logicAbsorber = nullptr;

    DefineCrystalSurface();
  }

  fScoringVolume = logicDetector;

  return physWorld;
}

void MyDetectorConstruction::ConstructSDandField() {
  G4SDManager* sdManager = G4SDManager::GetSDMpointer();

  auto* sd = new MySensitiveDetector("MySD");
  sdManager->AddNewDetector(sd);

  if (logicDetector) {
    logicDetector->SetSensitiveDetector(sd);
  }
}
