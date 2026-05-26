#ifndef DETECTOR_HH
#define DETECTOR_HH

#include "G4VSensitiveDetector.hh"
#include "G4SystemOfUnits.hh"
#include "MyHit.hh"
#include <unordered_map>

/*
  MySensitiveDetector:
  - Maintains a hits collection per event (MyHitCollection).
  - Also records a per-event map of trackID -> initial kinetic energy (vertex KE).
  - Provides a thread-local Instance() accessor so EventAction can read the per-event
    track initial energy map and write it to the new Tracks tree.
*/

class MySensitiveDetector : public G4VSensitiveDetector {
public:
  explicit MySensitiveDetector(const G4String &name);
  ~MySensitiveDetector() override;

  void Initialize(G4HCofThisEvent* hce) override;
  G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;
  void EndOfEvent(G4HCofThisEvent* hce) override;

  // Return the per-event map of track initial kinetic energies
  const std::unordered_map<G4int, G4double>& GetTrackInitKEMap() const { return fTrackInitKEMap; }

  // Thread-local instance accessor (set in constructor). EventAction uses this.
  static MySensitiveDetector* Instance() { return sInstance; }

private:
  MyHitCollection* hitsCollection;
  G4int HCID = -1;
  G4String myName;

  // Fast lookup for same-event accumulation: detectorID -> pointer to MyHit inside hitsCollection
  std::unordered_map<G4int, MyHit*> fHitMap;

  // Per-event map: trackID -> initial kinetic energy (vertex KE) recorded first time the track appears.
  std::unordered_map<G4int, G4double> fTrackInitKEMap;

  // thread-local pointer to the current SD instance (one per thread). Not owning.
  static thread_local MySensitiveDetector* sInstance;
};

#endif
