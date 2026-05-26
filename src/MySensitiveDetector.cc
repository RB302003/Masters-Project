#include "MySensitiveDetector.hh"
#include "MyHit.hh"
#include "G4UnitsTable.hh"

#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4VTouchable.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"

thread_local MySensitiveDetector* MySensitiveDetector::sInstance = nullptr;

MySensitiveDetector::MySensitiveDetector(const G4String &name)
  : G4VSensitiveDetector(name),
    hitsCollection(nullptr),
    myName(name)
{
  collectionName.insert("MyHitCollection");
  sInstance = this;
}

MySensitiveDetector::~MySensitiveDetector() {
  if (sInstance == this) sInstance = nullptr;
}

void MySensitiveDetector::Initialize(G4HCofThisEvent* hce)
{
  hitsCollection = new MyHitCollection(SensitiveDetectorName, collectionName[0]);

  if (HCID < 0) {
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  }

  if (hce && HCID >= 0) {
    hce->AddHitsCollection(HCID, hitsCollection);
  }

  fHitMap.clear();
  fTrackInitKEMap.clear();
}

G4bool MySensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  G4double edep = step->GetTotalEnergyDeposit();
  G4double stepLength = step->GetStepLength();

  if (stepLength <= 0.0) return false;

  G4Track* track = step->GetTrack();

  // Scintillation signal: proportional to energy deposited
  G4double scintYield = 10000.0 / MeV;
  G4double scintSignal = scintYield * edep;

  // Cherenkov signal: proportional to path length above threshold
  G4double cherenkovSignal = 0.0;
  G4double beta = track->GetVelocity() / CLHEP::c_light;

  // Determine refractive index from material
  G4String matName = step->GetPreStepPoint()->GetMaterial()->GetName();
  G4double n = 2.2; // default: PbWO4
  if (matName == "G4_PLASTIC_SC_VINYLTOLUENE") {
    n = 1.8; // Scintillator (sampling mode)
  }

  if (beta * n > 1.0) {
    G4double factor = 1.0 - 1.0 / (beta * beta * n * n);
    cherenkovSignal = stepLength * factor;
  }

  const G4VTouchable* touch = step->GetPreStepPoint()->GetTouchable();
  G4int detID = touch->GetCopyNumber();
  G4int trackID = track->GetTrackID();

  if (fTrackInitKEMap.find(trackID) == fTrackInitKEMap.end()) {
    G4double initKE = track->GetVertexKineticEnergy();
    fTrackInitKEMap.emplace(trackID, initKE);
  }

  auto it = fHitMap.find(detID);

  if (it == fHitMap.end()) {

    MyHit* hit = new MyHit(detID, -1, 0.0);

    hit->SetTrackID(trackID);
    hit->SetTime(step->GetPreStepPoint()->GetGlobalTime());
    hit->SetPosition(step->GetPreStepPoint()->GetPosition());
    hit->SetMomentum(step->GetPreStepPoint()->GetMomentum());

    hit->AddTrackEdep(trackID, edep);
    hit->AddScint(scintSignal);
    hit->AddCher(cherenkovSignal);

    hitsCollection->insert(hit);
    fHitMap.emplace(detID, hit);

  } else {

    MyHit* hit = it->second;

    hit->AddTrackEdep(trackID, edep);
    hit->AddScint(scintSignal);
    hit->AddCher(cherenkovSignal);
  }

  return true;
}

void MySensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{
}
