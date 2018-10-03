/** @file EicAnalysis.cxx
    @author Kolja Kauder
    @version Revision 0.1
    @brief Class for Jet substructure analysis
    @details Uses JetAnalyzer objects to perform A<SUB>J</SUB> analysis.
    @date Mar 02, 2015
*/

#include "EicAnalysis.hh"
#include <stdlib.h>     // for getenv, atof, atoi
#include<map>

#include <TParticlePDG.h>
#include "eicsmear/erhic/Pid.h"

using std::vector;
using std::map;
using std::cout;
using std::cerr;
using std::endl;

// ===================================================================
/** 
    For PseudoJet arithmetic
 */
PseudoJet operator-(const PseudoJet &p)
{
  return PseudoJet(-p.px(), -p.py(), -p.pz(), p.E());
}

// ===================================================================
// Rotations (or Breit Frame)
//-------------------------------------------------------------------
// rotation around the x axis
PseudoJet rotateX(const PseudoJet p, const double& psi){
  double sp = sin(psi);
  double cp = cos(psi);
  return PseudoJet(p.px(),
		   cp*p.py()-sp*p.pz(),
		   sp*p.py()+cp*p.pz(),
		   p.E());
}

// rotation around the y axis
PseudoJet rotateY(const PseudoJet p, const double& psi){
  double sp = sin(psi);
  double cp = cos(psi);
  return PseudoJet(cp*p.px()+sp*p.pz(),
		   p.py(),
		   cp*p.pz()-sp*p.px(),
		   p.E());
}

// rotation around the z axis
PseudoJet rotateZ(const PseudoJet p, const double& psi){
  double sp = sin(psi);
  double cp = cos(psi);
  return PseudoJet(cp*p.px()-sp*p.py(),
		   sp*p.px()+cp*p.py(),
		   p.pz(),
		   p.E());
}

// ===================================================================

// Standard ctor
EicAnalysis::EicAnalysis ( const int argc, const char** const argv )
{
  // Parse arguments
  // ---------------
  vector<string> arguments(argv + 1, argv + argc);
  bool argsokay=true;
  NEvents=-1;
  for ( auto parg = arguments.begin() ; parg!=arguments.end() ; ++parg){
    string arg=*parg;
    if ( arg == "-R" ){      
      if ( ++parg == arguments.end() ){ argsokay=false; break; }
      pars.R = atof( parg->data());
    } else if ( arg == "-breit" ){      
      if ( ++parg == arguments.end() ){ argsokay=false; break; }
      if ( parg->size() != 1 || (*parg)[0] < '0' || (*parg)[0] > '1' ) { argsokay=false; break; }
      pars.BreitFrame = ( (*parg)[0] == '1' );
    } else if ( arg == "-lja" ){     
      if (++parg==arguments.end() ){ argsokay=false; break; }
      pars.LargeJetAlgorithm = AlgoFromString ( *parg);
    } else if ( arg == "-rcja" ){
      if (++parg==arguments.end() ){ argsokay=false; break; }
      pars.CustomRecluster=true;
      pars.CustomReclusterJetAlgorithm=AlgoFromString ( *parg);
      ReclusterJetAlgorithm = AlgoFromString ( *parg);
      ReclusterJetDef = JetDefinition( ReclusterJetAlgorithm, 2.0*pars.R );
      recluster = new Recluster( ReclusterJetDef, 2.0*pars.R );
    } else if ( arg == "-bg" ){
      if (++parg==arguments.end() ){ argsokay=false; break; }
      pars.SubtractBg=BGTYPE(atoi(parg->data()));
      if ( pars.SubtractBg<0 || pars.SubtractBg>2 ){ argsokay=false; break; }
    } else if ( arg == "-Embbg" ){      
      if (++parg==arguments.end() ){ argsokay=false; break; }
      pars.EmbSubtractBg=BGTYPE(atoi(parg->data()));
      if ( pars.EmbSubtractBg<0 || pars.EmbSubtractBg>2 ){ argsokay=false; break; }
    } else if ( arg == "-b" ){
      if (++parg==arguments.end() ){ argsokay=false; break; }
      pars.beta = atof((parg)->data());
    } else if ( arg == "-z" ){
      if (++parg==arguments.end() ){ argsokay=false; break; }
      pars.z_cut = atof((parg)->data());
    } else if ( arg == "-pj" ){
      if (++parg ==arguments.end() ){ argsokay=false; break; }
      pars.PtJetMin = atof((parg)->data());
      if (++parg ==arguments.end() ){ argsokay=false; break; }
      pars.PtJetMax = atof((parg)->data());
    } else if ( arg == "-ec" ){      
      if (++parg ==arguments.end() ){ argsokay=false; break; }
      pars.EtaConsCut = atof((parg)->data());      
    } else if ( arg == "-pc" ){      
      if (++parg ==arguments.end() ){ argsokay=false; break; }
      pars.PtConsMin = atof((parg)->data());      
      if (++parg ==arguments.end() ){ argsokay=false; break; }
      pars.PtConsMax = atof((parg)->data());
    } else if ( arg == "-hadcorr" ){
      if ( ++parg == arguments.end() ){ argsokay=false; break; }
      pars.HadronicCorr = atof( parg->data());
    } else if ( arg == "-o" ){     
      if (++parg ==arguments.end() ){ argsokay=false; break; }
      pars.OutFileName=*parg;
    } else if ( arg == "-i" ){
      if (++parg ==arguments.end() ){ argsokay=false; break; }
      pars.InputName=*parg;
    } else if ( arg == "-embi" ){
      if (++parg ==arguments.end() ){ argsokay=false; break; }
      pars.EmbInputName=*parg;
      // Add a shortcut
      if (pars.EmbInputName == "FAKERHIC" ) {
	pars.EmbInputName = "Data/FakeAuAu20_*root";
	pars.Embintype = MCTREE;           ///< Embedding input type (can be a pico dst, a result tree, an MC tree)
	pars.EmbChainName = "tree";        ///< Name of the embedding input chain
	// pars.EmbTriggerName = "All";       ///< Embedding trigger type (All, MB, HT, pp, ppHT, ppJP)
      }
    } else if ( arg == "-c" ){
      if (++parg ==arguments.end() ){ argsokay=false; break; }
      pars.ChainName=*parg;
    } else if ( arg == "-trig" ){
      if (++parg ==arguments.end() ){ argsokay=false; break; }
      pars.TriggerName=*parg;
    } else if ( arg == "-embtrig" ){
      if (++parg ==arguments.end() ){ argsokay=false; break; }
      pars.EmbTriggerName=*parg;
    } else if ( arg == "-embc" ){
      if (++parg ==arguments.end() ){ argsokay=false; break; }
      pars.EmbChainName=*parg;
    } else if ( arg == "-nmix" ){      
      if (++parg==arguments.end() ){ argsokay=false; break; }
      pars.nMix=atoi(parg->data());
    } else if ( arg == "-intype" ){
      if (++parg ==arguments.end() ){ argsokay=false; break; }
      if ( *parg == "eptree" ){
	pars.intype = EPTREE;
	continue;
      }
      if ( *parg == "mctree" ){
	pars.intype = MCTREE;
	continue;
      }
      if ( *parg == "tree" ){
	pars.intype = INTREE;
	continue;
      }
      if ( *parg == "herwigtree" ){
	pars.intype = HERWIGTREE;
	continue;
      }
      argsokay=false;
      break;
    } else if ( arg == "-embintype" ){
      if (++parg ==arguments.end() ){ argsokay=false; break; }
      if ( *parg == "eptree" ){
	pars.Embintype = EPTREE;
	continue;
      }
      if ( *parg == "mctree" ){
	pars.Embintype = MCTREE;
	continue;
      }
      if ( *parg == "tree" ){
	pars.Embintype = INTREE;
	continue;
      }
      argsokay=false;
      break;
    } else if ( arg == "-N" ){      
      if (++parg==arguments.end() ){ argsokay=false; break; }
      NEvents=atoi(parg->data());
    } else if ( arg == "-tracksmear" ){
      if (++parg==arguments.end() ){ argsokay=false; break; }
      switch ( atoi(parg->data()) ){
      case 0:
	cout << "Track Smearing disabled." << endl;
	break;
      case 1:
	cout << "Track Smearing set to pp primaries. DeltapT/pT = 0.01+0.005 pT " << endl;
	SigmaPt = new TF1("SigmaPt","[0] + [1]*x",0,100);
	SigmaPt->FixParameter ( 0, 0.01);
	SigmaPt->FixParameter ( 1, 0.005);
	break;
      case 2:
	cout << "Track Smearing set to AA primaries. DeltapT/pT = 0.005+0.0025 pT " << endl;
	SigmaPt = new TF1("SigmaPt","[0] + [1]*x",0,100);
	SigmaPt->FixParameter ( 0, 0.005);
	SigmaPt->FixParameter ( 1, 0.0025);

	break;
      case 3:
	cout << "Track Smearing set to AA globals not implemented. DeltapT/pT = 0.01 * pT^2" << endl;
	SigmaPt = new TF1("SigmaPt","[0]*x*x",0,100);
	SigmaPt->FixParameter ( 0, 0.01);
	break;
      case 4:
	cout << "Unrecognized pT smearing option. " << endl;
	argsokay=false;
	break;	
      }
    } else if ( arg == "-fakeeff" ){
      if ( ++parg == arguments.end() ){ argsokay=false; break; }
      pars.FakeEff = atof( parg->data());
      cout << "Setting fake efficiency to " << pars.FakeEff << endl;
      if ( pars.FakeEff<0 || pars.FakeEff>1 ){ argsokay=false; break; }
    } else if ( arg == "-towunc" ){
      if ( ++parg == arguments.end() ){ argsokay=false; break; }
      pars.IntTowScale = atoi( parg->data());
      pars.fTowScale = 1.0 + pars.IntTowScale*pars.fTowUnc;
      cout << "Setting tower scale to " << pars.fTowScale << endl;
      if ( pars.IntTowScale<-1 || pars.FakeEff>1 ){ argsokay=false; break; }
    } else if ( arg == "-jetnef" ){
      if ( ++parg == arguments.end() ){ argsokay=false; break; }
      pars.MaxJetNEF = atof( parg->data());
      cout << "Setting Max Jet NEF to " << pars.MaxJetNEF << endl;
      if ( pars.MaxJetNEF<0 || pars.MaxJetNEF>1 ){ argsokay=false; break; }      
    } else {
      argsokay=false;
      break;
    }
  }
  
  if ( !argsokay ) {
    cerr << "usage: " << argv[0] << endl
	 << " [-o OutFileName]"  << endl
      	 << " [-N Nevents (<0 for all)]" << endl
      	 << " [-breit Use Breit frame (0,1)]" << endl
      	 << " [-lja LargeJetAlgorithm]" << endl
	 << " [-rcja ReclusteringJetAlgorithm]" << endl
	 << " [-bg BackgroundSubtraction (NONE=0, AREA=1, CONSTSUBPRE=2, CONSTSUBPOST=3)]" << endl
      	 << " [-Embbg EmbBackgroundSubtraction (NONE=0, AREA=1, CONSTSUBPRE=2, CONSTSUBPOST=3)]" << endl
      	 << " [-i infilepattern]" << endl
      	 << " [-c chainname]" << endl
      	 << " [-intype eptree|tree|mctree|herwigtree]" << endl
	 << " [-trig trigger name (e.g. HT)]" << endl
      	 << " [-ht offline high tower cut]" << endl
	 << " [-ht snap to high tower true|false ]" << endl
      	 << " [-embi embedding infilepattern]" << endl
      	 << " [-embc embedding chainname]" << endl
      	 << " [-embintype embedding eptree|tree|mctree]" << endl
	 << " [-nmix nMix]" << endl
	 << " [-R radius]" << endl
	 << " [-b beta]" << endl
	 << " [-z z_cut]" << endl
	 << " [-pj PtJetMin PtJetMax]" << endl
	 << " [-ec EtaConsCut]" << endl
	 << " [-pc PtConsMin PtConsMax]" << endl
	 << " [-hadcorr HadronicCorrection]  -- Set to a negative value for MIP correction." << endl
      	 << " [-psc PtSubConsMin PtSubConsMax]" << endl
	 << " [-tracksmear number] -- enable track pT smearing. " << endl
	 << "                      -- 1: pp primaries, 2: AuAu primaries, 3: AuAu (maybe also pp) globals." << endl
      	 << " [-fakeeff (0..1)] -- enable fake efficiency for systematics. 0.95 is a reasonable example." << endl
      	 << " [-towunc -1|0|1 ] -- Shift tower energy by this times " << pars.fTowUnc << endl
	 << endl << endl
	 << "NOTE: Wildcarded file patterns should be in single quotes." << endl
	 << endl;      
    throw std::runtime_error("Not a valid list of options");
  }

  // Consistency checks
  // ------------------
  if ( pars.SubtractBg==CONSTSUBPOST || pars.EmbSubtractBg==CONSTSUBPOST ){
    throw std::runtime_error("Constituent Subtraction after jetfinding not yet implemented, sorry.");
  }

  if ( pars.PtJetMin<=0 ){
    throw std::runtime_error("PtJetMin needs to be positive (0.001 will work).");
  }


  // Derived rapidity cuts
  // ---------------------
  EtaJetCut     = pars.EtaConsCut - pars.R;
  EtaGhostCut   = EtaJetCut  + 2.0*pars.R;
  
  // Jet candidate selectors
  // -----------------------
  // select_jet_eta     = SelectorAbsRapMax( EtaJetCut );
  select_jet_eta     = SelectorAbsEtaMax( EtaJetCut );
  select_jet_pt      = SelectorPtRange( pars.PtJetMin, pars.PtJetMax );
  select_jet         = select_jet_eta * select_jet_pt;     
  
  // Repeat on subjets?
  // ------------------
  pars.Recursive = pars.InputName.Contains("Pythia") && false;

  // Initialize jet finding
  // ----------------------
  AreaSpec = GhostedAreaSpec ( EtaGhostCut, pars.GhostRepeat, pars.GhostArea );
  AreaDef = AreaDefinition (fastjet::active_area_explicit_ghosts, AreaSpec);
  JetDef    = JetDefinition( pars.LargeJetAlgorithm, pars.R );

  SelectClose = fastjet::SelectorCircle( pars.R );

  cout << " R = " << pars.R << endl;
  cout << " beta = " << pars.beta << endl;
  cout << " z_cut = " << pars.z_cut << endl;
  cout << " Original jet algorithm : "<< pars.LargeJetAlgorithm << endl;
  if ( pars.CustomRecluster )
    cout << " Reclustering with a different algorithm "<< ReclusterJetAlgorithm << endl;
  cout << " PtJetMin = " << pars.PtJetMin << endl;
  cout << " PtJetMax = " << pars.PtJetMax << endl;
  cout << " PtConsMin = " << pars.PtConsMin << endl;
  cout << " PtConsMax = " << pars.PtConsMax << endl;
  cout << " Constituent eta cut = " << pars.EtaConsCut << endl;
  cout << " Jet eta cut = " << EtaJetCut << endl;
  cout << " Ghosts out to eta = " << EtaGhostCut << endl;
  cout << " Reading tree named \""<< pars.ChainName << "\" from " << pars.InputName << endl;
  cout << " intype = " << pars.intype << endl;
  cout << " Writing to " << pars.OutFileName << endl;
  if ( Embedding ) {
    cout << " Embedding into tree named \""<< pars.EmbChainName << "\" from " << pars.EmbInputName << endl;
    cout << " Embintype = " << pars.Embintype << endl;
    cout << " nMix = " << pars.nMix << endl;
  }
  cout << " ----" << endl;
  cout << " Writing to " << pars.OutFileName << endl;
  cout << " ----------------------------" << endl;

  // Quick and dirty QA histos
  // -------------------------
  QA_TowerEt = new TH2D( "QA_TowerEt","", 4800, 0.5, 4800.5, 80, 0, 80);

  // Provide a gaussian for track pt smearing
  // ----------------------------------------
  SmearPt = new TF1( "SmearPt","gaus(0)",-1,1);
}
//----------------------------------------------------------------------
EicAnalysis::~EicAnalysis(){
  if (pJA){
    delete pJA; pJA=0;
  }  
  if ( pConstituentBackgroundSubtractor ){  
    // delete pConstituentBackgroundSubtractor; // TAKEN CARE OFF IN pJA
    pConstituentBackgroundSubtractor =  0;
  }
  if ( pBackgroundSubtractor ){
    // delete pBackgroundSubtractor; // TAKEN CARE OFF IN pJA
    pBackgroundSubtractor =  0;
  }
}

//----------------------------------------------------------------------
bool EicAnalysis::InitChains(){

  Events = new TChain(pars.ChainName);
  Events->Add(pars.InputName);
  if ( NEvents<0 ) NEvents = INT_MAX;

  // Setup Input Event Buffer
  if ( pars.intype==EPTREE ){
    assert ( Events->GetEntries()>0 && "Something went wrong loading events.");
    NEvents=min(NEvents,Events->GetEntries() );
    
    // pFullEvent = new TClonesArray("TLorentzVector");
    inEvent = new erhic::EventPythia;
    Events->SetBranchAddress("event",&inEvent);
  }

  cout << "N = " << NEvents << endl;

  // For picoDSTs
  // -------------
  Embedding = false && (pars.EmbInputName !="" && pars.EmbInputName!="NONE");  

  // if ( Embedding ){
  //   EmbEvents = new TChain(pars.EmbChainName);
  //   EmbEvents->Add(pars.EmbInputName);    
  //   assert ( EmbEvents->GetEntries()>0 && "Something went wrong loading the embedding data.");
  //   if ( NEmbEvents<0 ) NEmbEvents = EmbEvents->GetEntries();
  //   NEmbEvents=min(NEmbEvents,EmbEvents->GetEntries() );
  //   gRandom->SetSeed(0);
  //   Embevi = gRandom->Integer(NEmbEvents); // Start at a random point

  //   if ( pars.Embintype==MCTREE ){

  //     pEmbEvent = new TClonesArray("TLorentzVector");
  //     EmbEvents->GetBranch("McParticles")->SetAutoDelete(kFALSE);
  //     EmbEvents->SetBranchAddress("McParticles", &pEmbEvent);
  //   } else if ( pars.Embintype==INTREE ){      

  //     pEmbEvent = new TClonesArray("TLorentzVector");
  //     EmbEvents->GetBranch("FullEvent")->SetAutoDelete(kFALSE);
  //     EmbEvents->SetBranchAddress("FullEvent", &pEmbEvent);
  //   } else {
  //     throw std::runtime_error("Unknown embedding type.");
  //   }

  //   cout << "Starting Embedding with event number " << Embevi << endl;
    
  // }  
    
  cout << "Done initializing chains. " << endl;
  return true;
}
//----------------------------------------------------------------------
// Main routine for one event.
EVENTRESULT EicAnalysis::RunEvent (){
  // cout << "-----------------------" << endl;
  // cout << "Entering EicAnalysis::RunEvent " << endl;
  TLorentzVector* sv;

  UInt_t filehash = 0;
  TString cname = "";
  
  // Reset results (from last event)
  // -------------------------------
  GroomingResult.clear();
  weight=1;
  partons.clear();
  particles.clear();

  if ( pHardPartons ) pHardPartons->Clear();
  if ( pHardPartonNames ) pHardPartonNames->Clear();

  switch (pars.intype) {
    // =====================================================
  case EPTREE:
    if ( evi>= NEvents ) {
      return EVENTRESULT::ENDOFINPUT;
      break;
    }
    if ( !(evi%500) ) cout << "Working on " << evi << " / " << NEvents << endl;
    inEvent->Clear();
    Events->GetEntry(evi);
    eventid = 0;
    runid = 0;

    ++evi;
    break;
  case INTREE :
  case MCTREE :
  case HERWIGTREE :
  default:
    cerr << "Unknown/unsupported intype " << pars.intype << endl;
    return EVENTRESULT::PROBLEM;
  }
    
  // // FIXME: May (will) not work as intended unless both inputs are picoDSTs!
  // runevent = ULong64_t(runid)*10000000LL + eventid;
    
  // Fill particle container
  // -----------------------

  int nparticles = 0;
  if ( pars.intype==EPTREE ){
    nparticles = inEvent->GetNTracks(); 
  }

  const Particle* inParticle =  nullptr;
  InitialBeam.clear();
  InitialEandGamma.clear();
  InitialRest.clear();
  FinalElectrons.clear();
  FinalRest.clear();
  Remainder.clear();

     
  if ( pars.intype == EPTREE )  {
    if ( nparticles < 5 ) {
      throw std::runtime_error("Not enough particles in the event");
    }
    
    for ( int i=0 ; i<nparticles ; ++i ){
      inParticle = inEvent->GetTrack(i);
      auto info = inParticle->Id().Info();

      // Let's keep track of unusual particles, which lead to problems with Id().Info()
      // check QCD effective states in http://home.fnal.gov/~mrenna/lutp0613man2/node44.html
      if ( !info ){
	if ( abs( inParticle->id )  > 9900000 && abs( inParticle->id )  < 9910000  // diffractive states (and exotics)
	     || abs( inParticle->id ) == 110 // reggeon
	     || abs( inParticle->id ) == 990 // pomeron
	     ) {
	  // cout << "Found PID " << inParticle->id
	  //      << " pt = " << inParticle->GetPt()
	  //      << " eta = " << inParticle->GetEta()
	  //      << " phi = " << inParticle->GetPhi()
	  //      << endl;
	} else {
	  // There are others that we're not catching
	  cout << "Found unexpected PID " << inParticle->id << endl;
	  throw std::runtime_error("Found unexpected PID");
	}
      }


      switch ( inParticle->GetStatus() ){
	// more or less pre-collision
	// --------------------------
      case 21 :
	if ( i==0 || i==1 ) {
	  InitialBeam.push_back ( *inParticle );
	  break;
	}
	
	if ( i==2 || i==3 ) {
	  InitialEandGamma.push_back ( *inParticle );
	  break;
	}

	InitialRest.push_back ( *inParticle );
	break;

	// Final state. This is what we could detect
	// -----------------------------------------
      case 1 : // status
	if ( inParticle->Id() == 11 ) FinalElectrons.push_back ( *inParticle );
	else                             FinalRest.push_back ( *inParticle );
	break;

	// Some intermediate state.
	// ------------------------
      case 11:
      case 12:
	Remainder.push_back ( *inParticle );
	break;
	
	// That should be all
	// ------------------
      default:
	cerr << " Unknown particle status " << inParticle->GetStatus() << endl;
	throw std::runtime_error("Unknown particle status.");
	break;
      }
		            
    } // particle loop

    // Sanity checks
    // -------------
    if ( InitialBeam.size() !=2 )  	 	throw std::runtime_error("InitialBeam.size != 2");
    if ( InitialBeam.at(0).Id() != 11 ) 	throw std::runtime_error("InitialBeam.at(0).Id() != 11");
    if ( InitialBeam.at(1).Id() != 2212 ) 	throw std::runtime_error("InitialBeam.at(1).Id() != 2212");
    
    if ( InitialEandGamma.size () !=2 )   	throw std::runtime_error("InitialEandGamma.size != 2");
    if ( InitialEandGamma.at(0).Id() != 11 ) 	throw std::runtime_error("InitialEandGamma.at(0).Id() != 11");
    if ( InitialEandGamma.at(1).Id() != 22 ) 	throw std::runtime_error("InitialEandGamma.at(1).Id() != 22");
    
    if ( InitialRest.at(0).Id() != 2212 )      	throw std::runtime_error("InitialRest.at(0).Id() != 2212"); 
   
    if ( FinalElectrons.size () ==0 )  	 	throw std::runtime_error("FinalElectrons is empty");
    for ( auto& e : FinalElectrons ) {
      if ( e.Id() != 11 ) 		throw std::runtime_error("FinalElectrons contains non-electron");
      if ( e.GetStatus() != 1 ) 	throw std::runtime_error("FinalElectrons contains non-final particle");
    }
    
    if ( FinalRest.size () ==0 )  	 	throw std::runtime_error("Empty final state");
    for ( auto& r : FinalRest ) {
      if ( r.GetStatus() != 1 ) 	throw std::runtime_error("FinalRest contains non-final particle");
    }
    
    for ( auto& r : Remainder ) {
      if ( r.GetStatus() != 11 && r.GetStatus() != 12 ) 	throw std::runtime_error("Unknown status in Remainder");
    }
    
  } else { // different input format    			    
    cerr << "Anything other than Brian's ep trees is currently not accepted." << endl;
    throw std::runtime_error("Found intype");
  }

  // Determine kinematics
  Q2 = -999;
  trueQ2 = inEvent->GetTrueQ2();
  X = -999;
  trueX = inEvent->GetTrueX();
  switch ( pars.KinematicsMethod ){
  case ELECTRON :
    X = inEvent->GetX();
    Q2 = inEvent->GetQ2();
  case JACQUETBLONDEL :
    X = inEvent->GetXJacquetBlondel();
    break;
  case DOUBLEANGLE : 
    X = inEvent->GetXDoubleAngle();
    break;
  case TRUEKINEMATICS :
    X = inEvent->GetTrueX();
  default :
    throw std::runtime_error("Unknown KinematicsMethod");    
  }



  // // Some experiments about kinematics
  // // ---------------------------------
  // auto& e1 = InitialBeam.at(0);
  // auto& p1 = InitialBeam.at(1);
  // auto& e2 = FinalElectrons.at(0);

  // auto s = 4 * e1.GetE() * p1.GetE();
  // // Following https://arxiv.org/pdf/0911.0884.pdf
  // // Electron method, eq. (10)
  // auto SigmaE = e2.GetE() * (1.-std::cos(e2.GetTheta() ));
  // auto yE = 1. -  SigmaE/2./e1.GetE();
  // auto Q2E = pow(e2.GetPt(),2) / (1-yE);
  // auto xE = Q2E / yE / s ;

  // Works similarly using the sum of all hadrons (eq 11)
  // SigmaH = Sum (Ei − pz,i) etc.
  // --> Jacquet-Blondel
  
  // cout << "root (s) " << sqrt(s) << endl;
  // cout << "y = " << y << endl;  
  // cout << Events->GetLeaf("trueY")->GetValue() << endl;
  // cout << "Q^2 = " << Q2 << endl;
  
  
  
  if (trueQ2>pow(2,2) ) {
    // cout << evi << "  " << Events->GetLeaf("trueQ2")->GetValue() << endl;
  }
  else return EVENTRESULT::NOTACCEPTED; // Let's focus on interesting events.

  // hard partons?
  //cout << " ==================================================" << endl;
  for ( auto& parton : Remainder ){
    // ignore hadrons 
    auto pdg = abs( parton.Id() );
    if ( pdg > 100 &&  pdg < 1000 ) continue; // mesons
    if ( pdg > 2000 &&  pdg < 6000 ) continue; // baryons
    
    //parton.Print();
    //cout << "  --> pt = " << parton.GetPt() << "  y = " << parton.GetRapidity() << "  eta = " << parton.GetEta() << "  phi = " << parton.GetPhi() << endl;
  }
  
  // cout << "x = " << x << endl;  
  // cout << Events->GetLeaf("trueX")->GetValue() << endl;
  // cout  << endl;

  // Set Up Boost to Breit Frame
  // ---------------------------
  auto gammaOrig = InitialEandGamma.at(1);
  auto protonOrig = InitialBeam.at(1); // line 2, what Brian uses
  // auto proton = InitialRest.at(0); // line 5
  if ( gammaOrig.Id() != 22 )      	throw std::runtime_error("gamma.Id() != 22");
  if ( protonOrig.Id() != 2212 )      	throw std::runtime_error("proton.Id() != 2212");

  // // e = e' + gamma ?
  // auto eglv = e2.Get4Vector()+gammaOrig.Get4Vector();
  // auto esum = e2.Get4Vector() + eglv;
  // if ( fabs ( eglv.E() - e1.GetE() ) > 1e-4 ) {
  //   cerr << "Virtual photon doesn't preserve energy " << endl;
  //   throw std::runtime_error("Virtual photon doesn't preserve energy");
  // }
  
  auto gamma = PseudoJet ( gammaOrig.Get4Vector() ) ;
  auto proton = PseudoJet ( protonOrig.Get4Vector() ) ;

  // TODO: Check this math!
  map <KINEMATICSMETHOD,PseudoJet> boosts;
  map <KINEMATICSMETHOD,PseudoJet> boostedprotons;
  // electron method
  PseudoJet boost_vectorE  = -(gamma + 2.0*inEvent->GetX()*proton);
  PseudoJet boosted_protonE = proton;
  boosted_protonE.boost(boost_vectorE);
  boosts[ELECTRON]=boost_vectorE;
  boostedprotons[ELECTRON]=boosted_protonE;

  // Jacquet-Blondel method
  PseudoJet boost_vectorJB = -(gamma + 2.0*inEvent->GetXJacquetBlondel()*proton);
  PseudoJet boosted_protonJB = proton;
  boosted_protonJB.boost(boost_vectorJB);
  boosts[JACQUETBLONDEL]=boost_vectorJB;
  boostedprotons[JACQUETBLONDEL]=boosted_protonJB;

  // DA method
  PseudoJet boost_vectorDA = -(gamma + 2.0*inEvent->GetXDoubleAngle()*proton);
  PseudoJet boosted_protonDA = proton;
  boosted_protonDA.boost(boost_vectorDA);
  boosts[DOUBLEANGLE]=boost_vectorDA;
  boostedprotons[DOUBLEANGLE]=boosted_protonDA;

  // True kinematics
  PseudoJet boost_vectorT  = -(gamma + 2.0*inEvent->GetTrueX()*proton);
  PseudoJet boosted_protonT = proton;
  boosted_protonT.boost(boost_vectorT);
  boosts[TRUEKINEMATICS]=boost_vectorT;
  boostedprotons[TRUEKINEMATICS]=boosted_protonT;

  // Works?
  bool checkboosts=true;
  if ( checkboosts ){
    for ( auto& toCheckPair : boosts ){
      auto & toCheck = toCheckPair.second;
      double bxL = toCheck.px()/toCheck.E();
      double byL = toCheck.py()/toCheck.E();
      double bzL = toCheck.pz()/toCheck.E();    
      if(bxL*bxL + byL*byL + bzL*bzL >= 1.0) {
	cout << "Bad boost type " << toCheckPair.first << " in Event = " << evi << endl;
	cout << bxL << " " << byL << " " << bzL << endl;
	std::runtime_error("Bad boost.");
      }
    }
  }

  // Let's put particles together!
  // -----------------------------
  vector<PseudoJet> lab_particles;
  vector<PseudoJet> boosted_particles;

  // Choose method used to determine Breit boost
  auto boosttype = pars.KinematicsMethod;
  auto breitboost = boosts[boosttype];
  auto breitproton = boostedprotons[boosttype];
  double phi_proton = breitproton.phi();
  double theta_proton = TMath::ATan2(breitproton.perp(),breitproton.pz());   // not provided by PseudoJet

  // For now, only exclude the correct electron
  // FIXME: Repeat for all
  for ( auto i = 1; i<FinalElectrons.size (); ++i ){
    auto& e = FinalElectrons.at(i);

    // CUTS
    if ( e.GetPt()< pars.PtConsMin )             continue;
    if ( fabs( e.GetEta() )>pars.EtaConsCut )    continue;

    PseudoJet pj =  PseudoJet ( e.Get4Vector() );
    lab_particles.push_back( pj );

    // boost to Breit frame
    // TODO: Check this math!    
    PseudoJet boostedpj =  pj;
    boostedpj.boost ( breitboost );
    boostedpj = rotateZ(boostedpj, -phi_proton);
    boostedpj = rotateY(boostedpj, -theta_proton);
    boosted_particles.push_back( boostedpj);
  }

  // Remaining hadrons
  for ( auto& r : FinalRest ) {
    // CUTS
    if ( r.GetPt()< pars.PtConsMin )             continue;
    if ( fabs( r.GetEta() )>pars.EtaConsCut )    continue;

    PseudoJet pj =  PseudoJet ( r.Get4Vector() );
    lab_particles.push_back( pj );

    // boost to Breit frame
    // TODO: Check this math!    
    PseudoJet boostedpj =  pj;
    boostedpj.boost ( breitboost );
    boostedpj = rotateZ(boostedpj, -phi_proton);
    boostedpj = rotateY(boostedpj, -theta_proton);
    boosted_particles.push_back( boostedpj);

  }

  // Which ones to use?
  if ( pars.BreitFrame ) {
    particles = boosted_particles;
  } else {
    particles = lab_particles;
  }
  
  if ( particles.size()==0 ) return EVENTRESULT::NOCONSTS;

  // For pythia, use cross section as weight
  // ---------------------------------------
  if ( TParameter<double>* sigmaGen=(TParameter<double>*) Events->GetCurrentFile()->Get("sigmaGen") ){
    weight=sigmaGen->GetVal();
  }
  
  // For HERWIG, have e-by-e cross section
  // -------------------------------------
  if ( Events->GetLeaf("weight") ){
    weight=Events->GetLeaf("weight")->GetValue();
  }
  
  // For pythia, fill leading parton container
  // -----------------------------------------  
  if ( pHardPartons ){
    for ( int i=0 ; i<pHardPartons->GetEntries() ; ++i ){
      sv = (TLorentzVector*) pHardPartons->At(i);
      
      PseudoJet pj = PseudoJet (*sv );
      
      // // flavor info
      // TString& s = ((TObjString*)(pHardPartonNames->At(i)))->String();
      // int qcharge=-999;
      // if ( s=="g" ) qcharge = 0;
      
      // if ( s(0)=='u' || s(0)=='c' || s(0)=='t' ) qcharge  = 2;
      // if ( s(0)=='d' || s(0)=='s' || s(0)=='b' ) qcharge = -1;
      // if ( s.Contains("bar") ) qcharge*=-1;
      
      // if ( abs ( qcharge ) >3 ) cout<< s << endl;
      // pj.set_user_info ( new JetAnalysisUserInfo( qcharge ) );

      partons.push_back( pj );
      
      // Save them too
      TLorentzVector* ppart = (TLorentzVector*) pSavedHardPartons->ConstructedAt( i );
      *ppart = *sv;

      // TObjString* saves = (TObjString*) pSavedHardPartonNames->ConstructedAt( i );
      // *saves=s.Data();
    }
  }

  
  // Run analysis
  // ------------
  if (pJA){
    delete pJA; pJA=0;
  }
  if ( pConstituentBackgroundSubtractor ){
    // delete pConstituentBackgroundSubtractor; // TAKEN CARE OFF IN pJA
    pConstituentBackgroundSubtractor =  0;
  }
  if ( pBackgroundSubtractor ){
    // delete pBackgroundSubtractor;// TAKEN CARE OFF IN pJA
    pBackgroundSubtractor =  0;
  }
  rho=0;
  
  // cout << particles.size() << endl; //return true;
  // if ( particles.size()>40 ) {
  //   for (auto&p : particles)    cout << p << endl;
  // }
    
  switch ( pars.SubtractBg ){
  case AREA:
    pJA = new JetAnalyzer( particles, JetDef, AreaDef, select_jet_eta * (!fastjet::SelectorNHardest(2)) ) ;
    pBackgroundSubtractor = pJA->GetBackgroundSubtractor();
    rho = pJA->GetBackgroundEstimator()->rho();
  break;
  case CONSTSUBPRE :
    // estimate the background    
    pJA = new JetAnalyzer( particles, JetDef, AreaDef, select_jet_eta * (!fastjet::SelectorNHardest(2)) ) ;
    pConstituentBackgroundSubtractor =  pJA->GetConstituentBackgroundSubtractor();
    // cout << pConstituentBackgroundSubtractor->description() << endl;

    // if ( ForceRho >=0 ) {
    //   std::cerr << "Can't currently handle ForceRho" << std::endl;
    //   throw(-1);
    //   // BackgroundSubtractor = new fastjet::contrib::ConstituentSubtractor( ForceRho );
    // }

    // Now correct the full event and create a new JetAnalyzer with these
    rho = pJA->GetBackgroundEstimator()->rho();     // save rho first
    particles = pConstituentBackgroundSubtractor->subtract_event(particles, pars.EtaConsCut);
    delete pJA; pJA=0;
    pJA = new JetAnalyzer( particles, JetDef, AreaDef ) ;
    
    break;
  case NONE:
    pJA = new JetAnalyzer( particles, JetDef ) ;
    // return EVENTRESULT::JETSFOUND;
    break;
  default :
    cerr << "Unsupported Background subtraction method" << endl;
    return EVENTRESULT::PROBLEM;
    break;      

  }

  // cout << particles.size() << endl;
  
  JetAnalyzer& JA = *pJA;
  vector<PseudoJet> JAResult = sorted_by_pt( select_jet ( JA.inclusive_jets() ) );
  // for ( auto& jet : JAResult ){
  //   if ( fabs(jet.eta())>0.6 ){
  //     cout << jet << endl;
  //   }
  // }
  
  if ( JAResult.size()==0 ) {
    // cout << "Nothing found" << endl;
    return EVENTRESULT::NOJETS;
  } 

  contrib::SoftDrop sd( pars.beta, pars.z_cut);
  if ( pars.CustomRecluster ) {
    sd.set_reclustering(true, recluster);
  }
  
  if ( pars.SubtractBg == AREA) {
    sd.set_subtractor(pBackgroundSubtractor);
    sd.set_input_jet_is_subtracted( false );
  }
  // cout << pBackgroundSubtractor->description() << endl;
  // sd.set_input_jet_is_subtracted( true );
  // sd.set_subtractor( 0 );
  // contrib::SoftDrop::_verbose=true;

  int njets = JAResult.size();
  // cout << "-----------------------" << endl;
  // cout << "We have " << njets << " jets. " << endl;

  for (unsigned ijet = 0; ijet < JAResult.size(); ijet++) {
    PseudoJet& CurrentJet = JAResult[ijet];

    // charged jets are problematic
    // // Check whether it fulfills neutral energy fraction criteria
    // // This could be done earlier in the jet selector
    // PseudoJet NeutralPart = join ( OnlyNeutral( CurrentJet.constituents() ) );
    // PseudoJet ChargedPart = join ( OnlyCharged( CurrentJet.constituents() ) );
    // double q = 0;
    // for ( PseudoJet& c : ChargedPart.constituents() ){
    //   q+= c.user_info<JetAnalysisUserInfo>().GetQuarkCharge();
    // }
    
    // JetAnalysisUserInfo* userinfo = new JetAnalysisUserInfo( q );
    // // Save neutral energy fraction in multi-purpose field
    // userinfo->SetNumber(NeutralPart.pt()  / CurrentJet.pt());
    // CurrentJet.set_user_info ( userinfo );

    // if ( pars.MaxJetNEF<1.0 &&  NeutralPart.pt()  / CurrentJet.pt() > pars.MaxJetNEF ) continue;

    // DEBUG - Recluster now
    // Well, good-ish. That doesn't make a difference in pt
    // cout << CurrentJet.pt() << " ---> ";
    // JetAnalysisUserInfo userinfo = CurrentJet.user_info<JetAnalysisUserInfo> ();
    // JetAnalysisUserInfo* newuserinfo = new JetAnalysisUserInfo( userinfo.GetQuarkCharge(), userinfo.GetTag(), userinfo.GetNumber() );
    // PseudoJet newjet = contrib::Recluster(cambridge_algorithm, JetDef.max_allowable_R)( CurrentJet );
    // if ( fabs(newjet.pt() - CurrentJet.pt()) > 1e-7 ){
    //   cout << CurrentJet.pt() << " ---> " << newjet.pt() << endl;
    // }
    // CurrentJet = newjet;
    // CurrentJet.set_user_info ( newuserinfo );
    
    // Run SoftDrop and examine the output
    PseudoJet sd_jet = sd( CurrentJet );

    // // DEBUG
    // vector<PseudoJet> temp=reshuffle (CurrentJet.constituents());
    // if ( CurrentJet.constituents().size() != temp.size() ){
    //   cerr << " !!!" << CurrentJet.constituents().size() << "  " << temp.size() << endl;
    //   throw -1;
    // }
    // JetDefinition TempJetDef    = JetDefinition( fastjet::cambridge_algorithm, 3*pars.R ); // gotta catch them all
    // JetAnalyzer TempJA (temp, TempJetDef);
    // vector<PseudoJet> shuffled = sorted_by_pt( select_jet ( TempJA.inclusive_jets() ) );
    // // if ( shuffled.size() !=1 ) cerr << shuffled.size() << endl;
    // // if ( shuffled.size() >1 ) cerr << " --> " << CurrentJet.pt() << "  " << shuffled.at(0).pt() << "  " << shuffled.at(1).pt()  << endl;
    // // if ( shuffled.size() >1 ) cerr << " --> " << CurrentJet.eta() << "  " << shuffled.at(0).eta() << "  " << shuffled.at(1).eta()  << endl;
    // // if ( shuffled.size() >1 ) cerr << " --> " << CurrentJet.constituents().size() << "  " << shuffled.at(0).constituents().size() << "  " << shuffled.at(1).constituents().size() << endl;
    // if ( shuffled.size() ==0 ) continue;
    // PseudoJet newjet = shuffled.at(0);
    // PseudoJet sd_jet = sd( newjet );
    // // END DEBUG
    
    // cout << CurrentJet.constituents().size() << endl;
    // cout << " Grooming Result: " << CurrentJet.pt() << "  --> " << sd_jet.pt() << endl << endl;
    if ( sd_jet == 0){
      cout <<  " FOREGROUND Original Jet:   " << CurrentJet << endl;
      if ( pBackgroundSubtractor ) cout <<  " FOREGROUND rho A: " << JA.GetBackgroundEstimator()->rho() * CurrentJet.area() << endl;	  
      if ( pBackgroundSubtractor ) cout <<  " FOREGROUND Subtracted Jet: " << (*pBackgroundSubtractor)( CurrentJet ) << endl;	  
      cout << " --- Skipped. Something caused SoftDrop to return 0 ---" << endl;
      continue;
    }

    if ( pars.SubtractBg == AREA) {
      // cout << "hello subtractor " << ijet << endl;
      CurrentJet = (*pBackgroundSubtractor)( CurrentJet );	
    }
    double zg = sd_jet.structure_of<contrib::SoftDrop>().symmetry();

    // if ( zg==0 ){
    //   // cout << CurrentJet.constituents().size() << endl;
    //   if ( CurrentJet.constituents().size() >5 ){
    // 	cout << " === " << endl;
    // 	for ( auto c : CurrentJet.constituents() ) cout << c << endl;
    // 	cout << " ====== " << endl;
    // 	cout << CurrentJet << endl;
    // 	cout << sd_jet << endl;
    //   }
    // }

    GroomingResult.push_back ( GroomingResultStruct ( CurrentJet, sd_jet, zg ) );

    // cout << CurrentJet.pt() - sd_jet.pt() << endl;
    // Shouldn't be negative
    if ( CurrentJet.pt() - sd_jet.pt() < -1e-7  && false )
      cout << "CurrentJet.pt() is smaller than sd_jet.pt()"
	   << CurrentJet.pt() - sd_jet.pt() << endl;
    
    
  }  
  // By default, sort for original jet pt
  sort ( GroomingResult.begin(), GroomingResult.end(), GroomingResultStruct::origptgreater);
  
  return EVENTRESULT::JETSFOUND;
}
  
//----------------------------------------------------------------------
double LookupXsec( TString filename ){
  // Some data for geant
  // -------------------
  //cross-sections for simulated GEANT data sample
  // From Renee.
  // also available via
  // http://people.physics.tamu.edu/sakuma/star/jets/c101121_event_selection/s0150_mclist_001/web.php
  // Double_t MinbXsec=28.12;
  // Double_t Xsec[12];
  // Xsec[0]=28.11;//2
  // Xsec[1]=1.287;//3
  // Xsec[2]=0.3117;//4
  // Xsec[3]=0.1360;//5
  // Xsec[4]=0.02305;//7
  // Xsec[5]=0.005494;//9
  // Xsec[6]=0.002228;//11
  // Xsec[7]=0.0003895;//15
  // Xsec[8]=0.00001016;//25
  // Xsec[9]=0.0000005010;//35
  // Xsec[10]=0.0000000283;//45
  // Xsec[11]=0.000000001443;//55
  
  static const Double_t MinbXsec=28.12;
  // static const Double_t Xsec[12] = {
  //   28.11,		// 2-3
  //   1.287,		// 3-4
  //   0.3117,		// 4-5
  //   0.1360,		// 5-7
  //   0.02305,		// 7-9
  //   0.005494,		// 9-11
  //   0.002228,		// 11-15
  //   0.0003895,		// 15-25
  //   0.00001016,		// 25-35
  //   0.0000005010,	// 35-45
  //   0.0000000283,	// 45-55
  //   0.000000001443	// 55-65
  // };

  static const Double_t Xsec[12] = {
    1.0,        // Placeholder for 2-3
    1.30E+09,	// 3-4
    3.15E+08,	// 4-5
    1.37E+08,	// 5-7
    2.30E+07,	// 7-9
    5.53E+06,	// 9-11
    2.22E+06,	// 11-15
    3.90E+05,	// 15-25
    1.02E+04,	// 25-35
    5.01E+02,	// 35-45
    2.86E+01,	// 45-55
    1.46E+00	// 55-65
  };

  static const Double_t Nmc[12] = {
    1,			// 2-3
    672518,		// 3-4
    672447,		// 4-5
    393498,		// 5-7
    417659,		// 7-9
    412652,		// 9-11
    419030,		// 11-15
    396744,		// 15-25
    399919,		// 25-35
    119995,		// 35-45
    117999,		// 45-55
    119999		// 55-65
  };

  Double_t w[12];
  for ( int i=0; i<12 ; ++i ){
    w[i] = Xsec[i] / Nmc[i];
    // w[i] = Nmc[i] / Xsec[i] ;
  }

  // static const Double_t w[12] = {
  //   1,			// Placeholder
  //   1.90E+03,
  //   6.30E+02,
  //   3.43E+02,
  //   5.49E+01,
  //   1.33E+01,
  //   5.30E+00,
  //   9.81E-01,
  //   2.56E-02,
  //   4.56E-03,
  //   2.43E-04,
  //   1.20E-05
  // };
    
  if ( filename.Contains("picoDst_3_4") ) return w[1];
  if ( filename.Contains("picoDst_4_5") ) return w[2];
  if ( filename.Contains("picoDst_5_7") ) return w[3];
  if ( filename.Contains("picoDst_7_9") ) return w[4];
  if ( filename.Contains("picoDst_9_11") ) return w[5];
  if ( filename.Contains("picoDst_11_15") ) return w[6];
  if ( filename.Contains("picoDst_15_25") ) return w[7];
  if ( filename.Contains("picoDst_25_35") ) return w[8];
  if ( filename.Contains("picoDst_35_45") ) return w[9];
  if ( filename.Contains("picoDst_45_55") ) return w[10];
  if ( filename.Contains("picoDst_55_65") ) return w[11];

  return 1;

}

//----------------------------------------------------------------------
/** 
    convenient output
*/
ostream & operator<<(ostream & ostr, const PseudoJet& jet) {
  if (jet == 0) {
    ostr << " 0 ";
  } else {
    ostr << " pt = " << jet.pt()
         << " m = " << jet.m()
         << " y = " << jet.rap()
         << " phi = " << jet.phi()
         << " ClusSeq = " << (jet.has_associated_cs() ? "yes" : "no");
  }
  return ostr;
}
//----------------------------------------------------------------------


