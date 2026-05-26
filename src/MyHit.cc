#include "MyHit.hh"
#include "G4THitsCollection.hh"
#include "G4ios.hh"
#include <algorithm>

G4ThreadLocal G4Allocator<MyHit>* MyHitAllocator = 0;

// Constructor
MyHit::MyHit(G4int id, G4int layer, G4double edep)
  : G4VHit(),
    fDetectorID(id),
    fLayerID(layer >= 0 ? layer : (id/10000)),
    fEdep(edep),
    fScint(0.0),
    fCher(0.0),
    fTime(0.0),
    fPosition(),
    fMomentum(),
    fTrackID(-1),
    fTrackEdepMap(),
    fTrackInitKE()
{}

MyHit::~MyHit() = default;

// copy ctor
MyHit::MyHit(const MyHit& right) : G4VHit(right) {
  fDetectorID    = right.fDetectorID;
  fLayerID       = right.fLayerID;
  fEdep          = right.fEdep;
  fScint         = right.fScint;
  fCher          = right.fCher;
  fTime          = right.fTime;
  fPosition      = right.fPosition;
  fMomentum      = right.fMomentum;
  fTrackID       = right.fTrackID;
  fTrackEdepMap  = right.fTrackEdepMap;
  fTrackInitKE   = right.fTrackInitKE;
}

const MyHit& MyHit::operator=(const MyHit& right) {
  if (this != &right) {
    fDetectorID    = right.fDetectorID;
    fLayerID       = right.fLayerID;
    fEdep          = right.fEdep;
    fScint         = right.fScint;
    fCher          = right.fCher;
    fTime          = right.fTime;
    fPosition      = right.fPosition;
    fMomentum      = right.fMomentum;
    fTrackID       = right.fTrackID;
    fTrackEdepMap  = right.fTrackEdepMap;
    fTrackInitKE   = right.fTrackInitKE;
  }
  return *this;
}

// equality
G4bool MyHit::operator==(const MyHit& right) const {
  return (this->fDetectorID == right.fDetectorID &&
          this->fLayerID == right.fLayerID &&
          this->fTrackID == right.fTrackID);
}

// --- Setters ---
void MyHit::SetID(G4int id) { fDetectorID = id; }
void MyHit::SetLayer(G4int layer) { fLayerID = layer; }
void MyHit::AddEdep(G4double edep) { fEdep += edep; }

void MyHit::AddScint(G4double val) { fScint += val; }
void MyHit::AddCher(G4double val) { fCher += val; }

void MyHit::SetTime(G4double time) { fTime = time; }
void MyHit::SetPosition(const G4ThreeVector& pos) { fPosition = pos; }
void MyHit::SetMomentum(const G4ThreeVector& mom) { fMomentum = mom; }
void MyHit::SetTrackID(G4int tid) { fTrackID = tid; }

// AddTrackEdep
void MyHit::AddTrackEdep(G4int trackID, G4double edep, G4double initKE) {
  if (edep <= 0.0) return;

  fEdep += edep;

  auto it = fTrackEdepMap.find(trackID);
  if (it == fTrackEdepMap.end()) {
    fTrackEdepMap.emplace(trackID, edep);
  } else {
    it->second += edep;
  }

  if (initKE >= 0.0) {
    auto itk = fTrackInitKE.find(trackID);
    if (itk == fTrackInitKE.end()) {
      fTrackInitKE.emplace(trackID, initKE);
    }
  }
}

// --- Getters ---
G4int MyHit::GetID() const { return fDetectorID; }
G4int MyHit::GetLayer() const { return fLayerID; }
G4double MyHit::GetEdep() const { return fEdep; }

G4double MyHit::GetScint() const { return fScint; }
G4double MyHit::GetCher() const { return fCher; }

G4double MyHit::GetTime() const { return fTime; }
G4ThreeVector MyHit::GetPosition() const { return fPosition; }
G4ThreeVector MyHit::GetMomentum() const { return fMomentum; }
G4int MyHit::GetTrackID() const { return fTrackID; }

const std::unordered_map<G4int, G4double>& MyHit::GetTrackEdepMap() const {
  return fTrackEdepMap;
}
const std::unordered_map<G4int, G4double>& MyHit::GetTrackInitKEMap() const {
  return fTrackInitKE;
}

size_t MyHit::GetTrackCount() const {
  return fTrackEdepMap.size();
}

std::vector<std::pair<G4int, G4double>> MyHit::GetTopNTracks(size_t N) const {
  std::vector<std::pair<G4int, G4double>> vec;
  vec.reserve(fTrackEdepMap.size());
  for (const auto &kv : fTrackEdepMap) vec.emplace_back(kv.first, kv.second);

  std::sort(vec.begin(), vec.end(),
            [](const std::pair<G4int, G4double>& a,
               const std::pair<G4int, G4double>& b){
              return a.second > b.second;
            });

  if (vec.size() > N) vec.resize(N);
  return vec;
}

std::pair<G4int, G4double> MyHit::GetMaxContributor() const {
  if (fTrackEdepMap.empty()) return std::make_pair(-1, 0.0);

  auto it = std::max_element(
      fTrackEdepMap.begin(), fTrackEdepMap.end(),
      [](const std::pair<G4int,G4double>& a,
         const std::pair<G4int,G4double>& b){
          return a.second < b.second;
      });

  return *it;
}
