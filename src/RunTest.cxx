/** @file RunTest.cxx
    @author Kolja Kauder
    @version Revision 0.1
    @brief Subjet Analysis for pythia and/or real events.
    @details Perform Subjet analysis on a chain of LorentzVectors
    @date Mar 23, 2017
*/

// #include "DevSubjetAnalysis.hh"
#include "AnalysisParameters.hh"
#include "EicAnalysis.hh"

#include "StEpSimuJet.h"
#include "StEpSimuJetParticle.h"

#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TChain.h>
#include <TFile.h>
#include <TBranch.h>
#include <TMath.h>
#include <TRandom.h>
#include <TParameter.h>
#include "TString.h"
#include "TObjString.h"

#include <set>
#include <vector>
#include <algorithm>

#include <cmath>
#include <climits>

// #include "fastjet/contrib/Recluster.hh"
#include "fastjet/contrib/SoftDrop.hh"

#include <exception>

using namespace std;
using namespace fastjet;
using namespace contrib;

// DEBUG
void decluster (PseudoJet j);

/** 
    - Parse parameters
    - Set up input tree
    - Set up output histos and tree
    - Initialize SubjetAnalysis object
    - If needed, combine input from two sources
    - Loop through events
    \arg argv: flags.
    Display options with
    <BR><tt>% RunTest -h </tt> 
    <BR>Note that wildcarded file patterns should be in single quotes.
*/

int main( int argc, const char** argv ){

  // shared_ptr<EicAnalysis> eic = nullptr;
  EicAnalysis* eic = nullptr; 
  try {
    // eic = make_shared<EicAnalysis>(argc, argv );
    eic = new EicAnalysis(argc, argv );
  } catch ( std::exception& e ){
    cerr << "Initialization failed with exception " << e.what() << endl;
    return -1;
  }

  // TF1* SigmaPt = new TF1("SigmaPt","[0] + [1]*x",0,100);
  
  // // // Kolja, primary in p+p:
  // // // sigma ( pT /pT )  ∼ 0.01 + 0.005*pT /(GeV/c)  //  --> 1% @ 0 GeV, 1.5% @ 1 GeV, 4% @ 6 GeV
  // // SigmaPt->FixParameter ( 0, 0.01);
  // // SigmaPt->FixParameter ( 1, 0.005);

  // // // Steven, primary in Au+Au:
  // // // sigma pT/pT  =  (0.5+0.25 pT)%. //  --> 0.5% @ 0 GeV, 1% @ 1 GeV, 2% @ 6 GeV
  // // SigmaPt->FixParameter ( 0, 0.005);
  // // SigmaPt->FixParameter ( 1, 0.0025);

  // // Nick: relative:
  // // 1% @ 1 GeV, 3% @ 5 GeV
  // // Nick: in quadrature --> 1.5% @ 1 GeV, 4% @ 5 GeV
  // SigmaPt->FixParameter ( 0, 0.00875);
  // SigmaPt->FixParameter ( 1, 0.00625);
  
  // // Jan, globals:
  // // Can't be right....
  // // sigma ( pT /pT ) ~ 0.01 * pT^2 // --> 0.0% @ 0 GeV, 1% @ 1 GeV, 3.6% @ 6 GeV
  // // TF1* SigmaPt = new TF1("SigmaPt","[0]*x*x",0,100);
  // // SigmaPt->FixParameter ( 0, 0.01);
    

  // eic->SetSigmaPt ( SigmaPt );
  
  if ( eic->InitChains() == false ){
    cerr << "Chain initialization failed" << endl;
    return -1;
  }
  
  // Get parameters we used
  // ----------------------
  const EicParameters pars  = eic->GetPars();

  // Files and histograms
  // --------------------

  TFile* fout = new TFile( pars.OutFileName, "RECREATE");
  
  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);
  TH3::SetDefaultSumw2(true);
    
  // TH1D* hzg      = new TH1D( "hzg", "z_{g}"              , 100*(pars.R+0.1),  0.0, (pars.R+0.1)  );
  // TH1D* hEmbzg   = new TH1D( "hEmbzg", "z_{g}, embedded" , 100*(pars.R+0.1),  0.0, (pars.R+0.1)  );
  TH3D* cptphieta = new TH3D("cptphieta","",500, 0.2, 50.2, 100, 0, TMath::TwoPi(), 100, -1, 1);
  TH3D* nptphieta = new TH3D("nptphieta","",500, 0.2, 50.2, 100, 0, TMath::TwoPi(), 100, -1, 1);
  
  // List of miscellaneous info
  // --------------------------
  TTree* info = new TTree("info", "Information");
  info->Branch("InputName"         , (void*)pars.InputName.Data()     , "InputName/C" );
  info->Branch("ChainName"         , (void*)pars.ChainName.Data()     , "ChainName/C" );
  bool Embedding = (pars.EmbInputName != "");
  if ( Embedding ){
    info->Branch("EmbInputName"         , (void*)pars.EmbInputName.Data()     , "EmbInputName/C" );
    info->Branch("EmbChainName"         , (void*)pars.EmbChainName.Data()     , "EmbChainName/C" );
  }
  info->Branch("R"                 , (void*) &pars.R                          , "R/D" );
  info->Branch("z_cut"             , (void*) &pars.z_cut                      , "z_cut/D" );
  info->Branch("beta"              , (void*) &pars.beta                       , "beta/D" );
  info->Branch("LargeJetAlgorithm" , (UInt_t*)&pars.LargeJetAlgorithm , "LargeJetAlgorithm/i" );

  if ( pars.CustomRecluster )
    info->Branch("ReclusterJetAlgorithm" , (UInt_t*)&pars.CustomReclusterJetAlgorithm , "ReclusterJetAlgorithm/i" );

  info->Branch("PtJetMin"          , (void*) &pars.PtJetMin                   , "PtJetMin/D" );
  info->Branch("PtJetMax"          , (void*) &pars.PtJetMax                   , "PtJetMax/D" );
  info->Branch("EtaConsCut"        , (void*) &pars.EtaConsCut                 , "EtaConsCut/D" );
  info->Branch("PtConsMin"         , (void*) &pars.PtConsMin                  , "PtConsMin/D" );
  info->Branch("PtConsMax"         , (void*) &pars.PtConsMax                  , "PtConsMax/D" );


  // Save results
  // ------------
  TTree* ResultTree=new TTree("ResultTree","Result Jets");
  
  TClonesArray HardPartons( "TLorentzVector" );
  ResultTree->Branch("HardPartons", &HardPartons );

  TClonesArray Jets( "StEpSimuJet" );
  ResultTree->Branch("Jets", &Jets );
  TClonesArray GroomedJets( "StEpSimuJet" );
  ResultTree->Branch("GroomedJets", &GroomedJets );

  // TClonesArray Jets( "TLorentzVector" );
  // ResultTree->Branch("Jets", &Jets );
  // TClonesArray GroomedJets( "TLorentzVector" ); 
  // ResultTree->Branch("GroomedJets", &GroomedJets );

  TClonesArray sj1( "TLorentzVector" );
  ResultTree->Branch("sj1", &sj1 );
  TClonesArray sj2( "TLorentzVector" );
  ResultTree->Branch("sj2", &sj2 );
  double refmult;
  ResultTree->Branch("refmult",&refmult, "refmult/d");  
  
  int njets=0;
  ResultTree->Branch("njets",   &njets, "njets/I" );
  double zg[1000];
  ResultTree->Branch("zg",       zg, "zg[njets]/D" );
  double delta_R[1000];
  ResultTree->Branch("delta_R",  delta_R, "delta_R[njets]/D" );
  double nef[1000];
  ResultTree->Branch("nef",  nef, "nef[njets]/D" );
  double mu[1000];
  ResultTree->Branch("mu",       mu, "mu[njets]/D" );
  double rho=-1;
  ResultTree->Branch("rho",      &rho, "rho/D" );

  double Q2=-999;
  ResultTree->Branch("Q2",      &Q2, "Q2/D" );
  double trueQ2=-999;
  ResultTree->Branch("trueQ2",  &trueQ2, "trueQ2/D" );
  double X=-999;
  ResultTree->Branch("X",       &X, "X/D" );
  double trueX=-999;
  ResultTree->Branch("trueX",   &trueX, "trueX/D" );



  double zg1[1000];
  double zg2[1000];
  if ( pars.Recursive){
    ResultTree->Branch("zg1",       zg1, "zg1[njets]/D" );
    ResultTree->Branch("zg2",       zg2, "zg2[njets]/D" );
  }

  double weight=1;
  ResultTree->Branch("weight",      &weight, "weight/D" );

  double Embrho=-1;
  ResultTree->Branch("Embrho",      &Embrho, "Embrho/D" );

  // Give each event a unique ID to compare event by event with different runs
  int runid;
  int eventid;
  ResultTree->Branch("eventid", &eventid, "eventid/I");
  ResultTree->Branch("runid",   &runid, "runid/I");

  // Helpers
  TLorentzVector* sv;
  TObjString* tobjs;
  
  // Go through events
  // -----------------
  // Long64_t Ntot=0;
  Long64_t Naccepted=0;
  // Long64_t Nrejected=0;
  cout << "Running analysis" << endl;
  try {
    bool ContinueReading = true;

    while ( ContinueReading ){

      HardPartons.Clear();
      Jets.Clear();
      GroomedJets.Clear();

      //   sj1.Clear();
      //   sj2.Clear();
      rho=-1;
      refmult=0;
      runid   =-(INT_MAX-1);
      eventid =-(INT_MAX-1);

      // for ( auto i=0; i<sizeof(zg) / sizeof(zg[0]); ++i ) zg[i]=0;
      std::fill_n( zg, sizeof(zg)/ sizeof(zg[0]), 0);
      EVENTRESULT ret=eic->RunEvent();

      // Understand what happened in the event
      switch (ret){
      case  EVENTRESULT::PROBLEM :
	cerr << "Encountered a serious issue" << endl;
	return -1;
	break;	
      case  EVENTRESULT::ENDOFINPUT :
	// cout << "End of Input" << endl;
	ContinueReading=false;
	continue;
	break;
      case  EVENTRESULT::NOTACCEPTED :
	// cout << "Event rejected" << endl;
	continue;
	break;
      case  EVENTRESULT::NOCONSTS :
	// cout << "Event empty." << endl;
	continue;
	break;
      case  EVENTRESULT::NOJETS :
	// cout << "No jets found." << endl;
	continue;
	break;
      case  EVENTRESULT::JETSFOUND:
	// The only way not to break out or go back to the top
	// Do Something
	Naccepted++;
	break;
      default :
	cerr << "Unknown return value." << endl;
	return -1;
	// I understand that a break after continue or return is silly...
	// But it's necessary in nested switches in root and I don't want to lose the habit    
	break;
      }
	  
      // cout << ++Ntot << endl;
      
      // Now we can pull out details and results
      // ---------------------------------------
      weight = eic->GetEventWeight();
      refmult = eic->GetRefmult();
      rho = eic->GetRho();
      runid = eic->GetRunid();
      eventid = eic->GetEventid();
      
      Q2 = eic->GetQ2();
      trueQ2 = eic->GetTrueQ2();
      X = eic->GetX();
      trueX = eic->GetTrueX();


      TClonesArray* pHardPartons =  eic->GetHardPartons();
      if ( pHardPartons ){
	sv = (TLorentzVector*) HardPartons.ConstructedAt( 0 ); *sv = *( (TLorentzVector*) pHardPartons->At(0));
	sv = (TLorentzVector*) HardPartons.ConstructedAt( 1 ); *sv = *( (TLorentzVector*) pHardPartons->At(1));
      }
      
      vector<GroomingResultStruct> GroomingResult = eic->GetGroomingResult(); 
      // sort ( GroomingResult.begin(), GroomingResult.end(), GroomingResultStruct::groomedptgreater);
      // sort ( GroomingResult.begin(), GroomingResult.end(), GroomingResultStruct::origptgreater);
            
      njets=GroomingResult.size();
      int ijet=0;
      for ( auto& gr : GroomingResult ){	
	auto& origjet = gr.orig;
	auto& groomedjet = gr.groomed;
	
	zg[ijet] = gr.zg;
	delta_R[ijet]=groomedjet.structure_of<contrib::SoftDrop>().delta_R();
	// nef[ijet] = gr.orig.user_info<JetAnalysisUserInfo>().GetNumber() ;

	// Old - using TLorentzVector
	// Also note that using constructed_at is better anyway
	// TLorentzVector sv = TLorentzVector( MakeTLorentzVector( gr.orig) );
	// sv.SetCharge( gr.orig.user_info<JetAnalysisUserInfo>().GetQuarkCharge() / 3 );      
 	// new ( Jets[ijet] )               TLorentzVector ( sv );
	// new ( GroomedJets[ijet] )        TLorentzVector ( TLorentzVector( MakeTLorentzVector( gr.groomed) ) );


	StEpSimuJet* Jet = (StEpSimuJet*) Jets.ConstructedAt( ijet );
	*Jet = StEpSimuJet ( origjet.pt(),origjet.eta(),origjet.phi(),origjet.E() );
	Jet->SetZg(gr.zg);
	Jet->SetRg( groomedjet.structure_of<contrib::SoftDrop>().delta_R() );
	Jet->SetNconsts(origjet.constituents().size());

	StEpSimuJet* GroomedJet = (StEpSimuJet*) GroomedJets.ConstructedAt( ijet );
	*GroomedJet = StEpSimuJet ( groomedjet.pt(),groomedjet.eta(),groomedjet.phi(),groomedjet.E() );
	GroomedJet->SetZg(gr.zg);
	GroomedJet->SetRg( groomedjet.structure_of<contrib::SoftDrop>().delta_R() );
	GroomedJet->SetNconsts(groomedjet.constituents().size());

	ijet++;
      }
      
      ResultTree->Fill();

      const vector<PseudoJet>& particles = eic->GetParticles();
      for (auto& p : particles ){
	if ( p.has_user_info<JetAnalysisConstituentInfo>() ){
	  if ( abs( p.user_info<JetAnalysisConstituentInfo>().GetQuarkCharge() )>0 ){
	    cptphieta->Fill(p.pt(),p.phi(),p.eta());
	  } else {
	    nptphieta->Fill(p.pt(),p.phi(),p.eta());
	  }
	}else{
	  // cout << "constituent: no charge info" << endl;
	}
      }	     
    }
  } catch (std::string& s) {
    cerr << "RunEvent failed with string " << s << endl;
    return -1;
  } catch ( std::exception& e ){
    cerr << "RunEvent failed with exception " << e.what() << endl;
    return -1;
  }

  
  // info->Branch("Ntot"      , &Ntot      , "Ntot/L" );
  info->Branch("Naccepted" , &Naccepted , "Naccepted/L" );
  // info->Branch("Nrejected" , &Nrejected , "Nrejected/L" );
  info->Fill();
  
  fout->Write();

  if (eic->QA_TowerEt ) eic->QA_TowerEt->Write();

  cout << "Done." << endl;


  delete eic;
  return 0;
}
//----------------------------------------------------------------------

// DEBUG
void decluster (PseudoJet j){
  cout << " ### Declustering ### " << endl;
  cout << j.pt() << endl;

  PseudoJet piece1, piece2;

  static int level=0;
  if ( j.has_parents(piece1, piece2) ){
    cout << "level = " << level++ << endl;
    cout << "piece1.pt() = " << piece1.pt() << endl;
    cout << "piece2.pt() = " << piece2.pt() << endl;

    if (! piece1.is_pure_ghost() ) decluster(piece1);
    if (! piece2.is_pure_ghost() ) decluster(piece2);

  } else cout << " Done with this branch" << endl;
  return;
}

//----------------------------------------------------------------------
bool readinbadrunlist(vector<int> & badrun, TString csvfile) {
	
  // open infile
  std::string line;
  std::ifstream inFile (csvfile );
	
  std::cout<<"Loading bad run id from "<< csvfile.Data()<<std::endl;;
	        
  if ( !inFile.good() ) {
    std::cout<<"Can't open "<<csvfile.Data()<<std::endl;
    return false;
  }
	
  while (std::getline (inFile, line) ){
    if ( line.size()==0 ) continue; // skip empty lines
    if ( line[0] == '#' ) continue; // skip comments
	
    std::istringstream ss( line );
    while( ss ){
      std::string entry;
      std::getline( ss, entry, ',' );
      int ientry = atoi(entry.c_str());
      if (ientry) {
	badrun.push_back( ientry );
	std::cout<<"Added bad runid "<<ientry<<std::endl;
      }
    }
  }
	
  return true;
}

