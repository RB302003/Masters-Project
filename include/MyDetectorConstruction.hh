#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH

#include "G4SystemOfUnits.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4GenericMessenger.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"

class MyDetectorConstruction : public G4VUserDetectorConstruction {
public:
  MyDetectorConstruction();
  ~MyDetectorConstruction();

  G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

  // Returns current mode string so other classes can query it
  const G4String& GetMode() const { return fMode; }

  virtual G4VPhysicalVolume *Construct();
  virtual void ConstructSDandField();

private:
  void DefineMaterial();
  void DefinePbWO4OpticalProperties();
  void DefineAirOpticalProperties();
  void DefineCrystalSurface();

  G4LogicalVolume *logicDetector;
  G4LogicalVolume *logicAbsorber;  // only used in sampling mode

  G4int nCols, nRows;
  G4int nLayers;            // used in homogeneous mode (24)
  G4int absorberLayers;     // used in sampling mode (20)
  G4int activeLayers;       // used in sampling mode (20)

  G4Box *solidWorld;

  G4LogicalVolume *logicWorld;

  G4VPhysicalVolume *physWorld;

  // Materials — all three always built
  G4Material *Air;
  G4Material *Aluminium;
  G4Material *PbWO4;         // homogeneous mode
  G4Material *Lead;          // sampling mode absorber
  G4Material *Scintillator;  // sampling mode active

  G4double xWorld, yWorld, zWorld;

  // Runtime mode switching
  G4GenericMessenger *fMessenger;
  G4String fMode;  // "homogeneous" or "sampling"

protected:
  G4LogicalVolume* fScoringVolume = nullptr;
};

#endif
