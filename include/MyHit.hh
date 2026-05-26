#ifndef MYHIT_HH
#define MYHIT_HH

#include "G4VHit.hh"
#include "G4ThreeVector.hh"
#include "G4THitsCollection.hh"
#include <unordered_map>
#include <vector>
#include <utility>

class MyHit;

using MyHitCollection = G4THitsCollection<MyHit>;

class MyHit : public G4VHit {
public:
  MyHit(G4int id = -1, G4int layer = -1, G4double edep = 0.0);
  MyHit(const MyHit& right);
  const MyHit& operator=(const MyHit& right);
  virtual ~MyHit();

  G4bool operator==(const MyHit& right) const;

  // --- Setters ---
  void SetID(G4int id);
  void SetLayer(G4int layer);

  void AddTrackEdep(G4int trackID, G4double edep, G4double initKE = -1.0);
  void AddEdep(G4double edep);

  // Signal adders
  void AddScint(G4double val);
  void AddCher(G4double val);

  void SetTime(G4double time);
  void SetPosition(const G4ThreeVector& pos);
  void SetMomentum(const G4ThreeVector& mom);
  void SetTrackID(G4int tid);

  // --- Getters ---
  G4int GetID() const;
  G4int GetLayer() const;
  G4double GetEdep() const;

  // Signal getters
  G4double GetScint() const;
  G4double GetCher() const;

  G4double GetTime() const;
  G4ThreeVector GetPosition() const;
  G4ThreeVector GetMomentum() const;
  G4int GetTrackID() const;

  const std::unordered_map<G4int, G4double>& GetTrackEdepMap() const;
  const std::unordered_map<G4int, G4double>& GetTrackInitKEMap() const;

  size_t GetTrackCount() const;

  std::vector<std::pair<G4int, G4double>> GetTopNTracks(size_t N) const;
  std::pair<G4int, G4double> GetMaxContributor() const;

private:
  G4int fDetectorID;
  G4int fLayerID;

  G4double fEdep;

  // Signal storage
  G4double fScint;
  G4double fCher;

  G4double fTime;
  G4ThreeVector fPosition;
  G4ThreeVector fMomentum;
  G4int fTrackID;

  std::unordered_map<G4int, G4double> fTrackEdepMap;
  std::unordered_map<G4int, G4double> fTrackInitKE;
};

#endif
