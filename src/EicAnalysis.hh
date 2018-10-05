/** @file EicAnalysis.hh
    @author Kolja Kauder
    @version Revision 0.1
    @brief Class for A<SUB>J</SUB> analysis
    @details Uses JetAnalyzer objects to perform z<SUB>g</SUB> analysis.
    @date Mar 02, 2015
*/

#ifndef __EICANALYSIS_HH
#define __EICANALYSIS_HH

#include "EicEnums.hh"
#include "AnalysisParameters.hh"
#include "JetAnalyzer.hh"

#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TString.h"
#include "TRandom.h"
#include "TChain.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TFile.h"
#include "TSystem.h"
#include "TParameter.h"
#include "TClonesArray.h"

// #include "fastjet/contrib/Recluster.hh"
#include "fastjet/contrib/SoftDrop.hh"

#include "eicsmear/erhic/EventBase.h"
#include "eicsmear/erhic/EventPythia.h"
#include "eicsmear/erhic/Particle.h"

#include <assert.h>
#include <iostream>
#include <cmath>
#include <climits>
#include <sstream>

using namespace std;
using namespace fastjet;
using namespace contrib;

#include <random>
#include <algorithm>

/** 
    A helper for geant data
 */
double LookupXsec( TString filename );

/**
   For sorting with a different key
*/
typedef pair<PseudoJet,double> PseudoJetPt;
struct PseudoJetPtGreater {
  bool operator()( PseudoJetPt const & a, PseudoJetPt const & b) { 
    return a.second > b.second ;
  }
};

/**
  To keep original and groomed jets connected
 */
class GroomingResultStruct{
public:
  PseudoJet orig;
  PseudoJet groomed;
  double zg;
  GroomingResultStruct ( PseudoJet orig, PseudoJet groomed, double zg )
    : orig(orig), groomed(groomed),zg(zg)
  {};
  
  static bool origptgreater( GroomingResultStruct const & a, GroomingResultStruct const & b) { 
    return a.orig.pt()>b.orig.pt();
  };
  
  static bool groomedptgreater( GroomingResultStruct const & a, GroomingResultStruct const & b) { 
    return a.groomed.pt()>b.groomed.pt();
  };
  
};


/** 
    convenient output
*/
ostream & operator<<(ostream & ostr, const PseudoJet& jet);

/**
   Selectors
*/
static const Selector NotGhost = !fastjet::SelectorIsPureGhost ();     ///< Helper useful outside the class as well
static const Selector OnlyCharged = NotGhost && ( SelectorChargeRange( -3, -1) || SelectorChargeRange( 1, 3) );     ///< Helper useful outside the class as well
static const Selector OnlyNeutral = NotGhost && SelectorChargeRange( 0, 0);     ///< Helper useful outside the class as well


/**
   The main class
 */
class EicAnalysis {

private :

  // These need to be initialized
  // ----------------------------
  EicParameters pars;   ///< container to have all analysis parameters in one place

  // Storage for sorting out the various MC particles
  // ------------------------------------------------
  // Standard pythia6 event record from Brian's files:
  // 1     21         11        0        3        4
  // 2     21       2212        0        5        0
  // 3     21         11        1        0        0
  // 4     21         22        1        0        0
  // 5     21       2212        2        0        0
  // -- 1 and 2 are the incoming e+P
  // -- 3 and 4 are the (still incoming) e+gamma (virt)
  // -- 5 is a copy of proton 2, can supposedly sometimes differ but I don't see it
  // -- then comes a number of E=0 (mostly) gluons with status 21.
  // -- the first final particle (typically 12) seems to always be "the" outgoing (final) electron
  // -- then comes the actual event, with statuses 1 (final), 11 and 12 (intermediate, but what's the difference?)
 
  vector <Particle> InitialBeam; //< lines 1-2
  vector <Particle> InitialEandGamma; //< lines 3-4
  vector <Particle> InitialRest; //< lines with status 21
  vector <Particle> FinalElectrons; //< expect the first one to be "the" electron (~line 12), but of course there may be more...
  vector <Particle> FinalRest; //< lines 13+ with status 1
  vector <Particle> Remainder; //< Everything from statuses 11, 12

  // Internal
  // --------  
  float EtaJetCut;       ///< jet eta
  float EtaGhostCut;     ///< ghost eta
  
  fastjet::JetDefinition JetDef;       ///< jet definition

  // Relevant jet candidates
  fastjet::Selector select_jet_eta;        ///< jet rapidity selector
  fastjet::Selector select_jet_pt;         ///< jet p<SUB>T</SUB> selector
  fastjet::Selector select_jet;            ///< compound jet selector

  fastjet::GhostedAreaSpec AreaSpec;      ///< ghosted area specification
  fastjet::AreaDefinition AreaDef;        ///< jet area definition

  // Data
  // ----
  Long64_t NEvents=-1;
  TChain* Events = 0;
  TClonesArray* pFullEvent = 0;             ///< Constituents
  TClonesArray* pHardPartons = 0;           ///< For pythia data, the original hard scatter

  double Q2;
  double trueQ2;
  double X;
  double trueX;

  vector<PseudoJet> particles;
  vector<PseudoJet> partons;
  double rho=0;                             ///< background density
  
  Long64_t evi=0;
  Long64_t NEmbEvents=-1;
  TChain* EmbEvents = 0;
  TClonesArray* pEmbEvent = 0;
  
  bool Embedding=false;
  Long64_t Embevi =0; // Starting point will be changed!  


  // Note that the following are NOT unique across MC trees
  // My fault, will try to fix with a quick and dirty trick:
  // runid (and eventid) are below 1M (itself not exactly optimal)
  // So we'll just hash the data file name and add it to the runid,
  // making sure that it stays in the right range
  int eventid;
  int runid;
  int Embeventid;
  int Embrunid;

  double refmult=0;
  double Embrefmult=0;

  double weight=0;    

  JetAnalyzer *pJA=0;
  JetAnalyzer *pEmbJA=0;
  Subtractor* pBackgroundSubtractor =  0;
  Subtractor* pEmbBackgroundSubtractor =  0;
  ConstituentSubtractor* pConstituentBackgroundSubtractor = 0;

  /// Resolution Smearing width for tracks
  /// This function should return sigma (pT - pT_truth) for a given pT
  TF1 * SigmaPt;
  /// Actual smearing factor
  /// Gaussian around 1 with the width from SigmaPt
  TF1 * SmearPt;

  // For matching
  // ------------
  fastjet::Selector SelectClose;

  JetAnalyzer* pJAhi;                      ///< JetAnalyzer object for high pT
  JetAnalyzer* pJAlo;                      ///< JetAnalyzer object for low pT
  
  // --------------------------------------------------------------
  // JetAlgorithm ReclusterJetAlgorithm = fastjet::cambridge_algorithm;
  fastjet::JetAlgorithm ReclusterJetAlgorithm;
  fastjet::JetDefinition ReclusterJetDef;
  Recluster * recluster=0; 

  // Results to be passed around
  TClonesArray *pSavedHardPartons=0;           ///< original hard scatter in PYTHIA
  TClonesArray *pSavedHardPartonNames=0;       ///< original hard scatter PID in PYTHIA
  vector<GroomingResultStruct> GroomingResult; ///< grooming result in a nice structured package

  erhic::EventPythia* inEvent;
    
public:  

  /** Standard constructor. Parse vectorized argv.
      \param argc: number of arguments
      \param argv: string array of command line options
   */
  EicAnalysis ( const int argc, const char** const  );

  /** Destructor. Clean things up
   */
  virtual ~EicAnalysis();

  /** Decoupled chain initialization for readability
   */
  bool InitChains ( );

  /** Main routine for one event.
      \return false if at the end of the chain      
   */
  EVENTRESULT RunEvent ( );
 
  // Quick and dirty QA histos - keep them public
  // --------------------------------------------
  TH2D* QA_TowerEt;
  
  // Getters and Setters
  // -------------------
  inline EicParameters& GetPars()         { return pars; };

  /// Get jet radius
  inline double GetR ( )                   { return pars.R; };
  /// Set jet radius
  inline void   SetR ( const double newv ) { pars.R=newv;   };

  /// Get the weight of the current event (mainly for PYTHIA)
  inline double GetEventWeight()           { return weight; };

  /// Get the refmult of the current event
  inline double GetRefmult()           { return refmult; };

  /// Get the runid of the current event
  inline double GetRunid()           { return runid; };

  /// Get the eventid of the current event
  inline double GetEventid()           { return eventid; };

  /// The main result of the analysis
  inline const vector<GroomingResultStruct>& GetGroomingResult()    { return GroomingResult; }

  inline TClonesArray* GetHardPartons()              { return pSavedHardPartons; }
  inline TClonesArray* GetHardPartonNames()          { return pSavedHardPartonNames; }
  
  inline double GetRho()                             { return rho; }
  inline const vector<PseudoJet>& GetParticles() const     { return particles; }
  
  /// Handle to selector for jet candidates
  inline fastjet::Selector& GetJetSelector () { return select_jet; }

  /// Handle to ghosted area specification
  inline fastjet::GhostedAreaSpec& GetAreaSpec () { return AreaSpec; }
  /// Handle to jet area definition
  inline fastjet::AreaDefinition& GetAreaDef () { return AreaDef; }

  /// Handle to constituents
  inline std::vector<fastjet::PseudoJet> GetConstituents() {return particles; };

  /// Current Q2
  inline double GetQ2()                             { return Q2; }
  /// Current true Q2
  inline double GetTrueQ2()                         { return trueQ2; }
  /// Current X
  inline double GetX()                              { return X; }
  /// Current true X
  inline double GetTrueX()                          { return trueX; }
   
  
};  

#endif // __EICANALYSIS_HH
