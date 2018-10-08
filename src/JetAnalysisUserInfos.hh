#ifndef JETANALYSISUSERINFOS_HH
#define JETANALYSISUSERINFOS_HH

/** @author Kauder:Kolja
    @version Revision 0.1
    @description Some UserInfo UserInfo derivatives.
    Only one can be attached per PseudoJet.
    Splitting into constituent info and jet info (to save space).
    @date Mar 03, 2015
*/


// =============================================================================

/** Constituent information
    For quick prototyping, you can use the string Tag and the float Number.
 */
class JetAnalysisConstituentInfo: public fastjet::PseudoJet::UserInfoBase {
public:
  /// Standard Constructor
  JetAnalysisConstituentInfo(int quarkcharge=-999, int pdg=0,
			     float number=-1, std::string tag="" )
    :  mQuarkcharge( quarkcharge ), mPdg (pdg),
       mNumber(number), mTag(tag)
  {};
  
  /// Charge in units of e/3
  int QuarkCharge() const { return mQuarkcharge; };
  /// Set charge in units of e/3
  void SetQuarkCharge(const int Quarkcharge )  { mQuarkcharge = Quarkcharge; };

  /// Get PDG code, quick reference at http://home.fnal.gov/~mrenna/lutp0613man2/node44.html
  int Pdg() const { return mPdg; };
  /// Set PDG code
  void SetPdg(const int pdg)  { mPdg=pdg; };

  /// Get multi-purpose description
  std::string Tag() const { return mTag; };
  /// Set multi-purpose description
  void SetTag( const std::string tag ) { mTag=tag; };

  /// Get multi-purpose numerical value
  float Number() const { return mNumber; };
  /// Set multi-purpose numerical value
  void SetNumber(const float f)  { mNumber=f; };

protected:
  int mQuarkcharge;    ///< Charge in units of e/3
  int mPdg;            ///< PDG code

  std::string mTag;    ///< Multi-purpose
  float mNumber;       ///< Multi-purpose
};

// =============================================================================

/** Jet information
    For quick prototyping, you can use the string Tag and the float Number.
 */
class JetAnalysisJetInfo: public fastjet::PseudoJet::UserInfoBase {
public:
  /// Standard Constructor
  JetAnalysisJetInfo(int quarkcharge=-999, 
		     double nef = -1,
		     double zg=-1,
		     double rg=-1,
		     float number=-1, std::string tag="" )
    :  mQuarkcharge( quarkcharge ), 
       mNef(nef),
       mZg(zg),
       mRg(rg),
       mNumber(number), mTag(tag)
  {};
  
  /// Charge in units of e/3
  int QuarkCharge() const { return mQuarkcharge; };
  /// Set charge in units of e/3
  void SetQuarkCharge(const int Quarkcharge )  { mQuarkcharge = Quarkcharge; };


  /// Get groomed momentum sharing z_g
  double Zg() const {return mZg;};
  /// Set groomed momentum sharing z_g
  void SetZg( const double Zg ) { mZg=Zg; };

  /// Get groomed radius R_g
  double Rg() const {return mRg;};
  /// Set groomed momentum sharing z_g
  void SetRg( const double Rg ) { mRg=Rg; };

  /// Get neutral energy fraction
  double Nef() const {return mNef;};
  /// Set groomed momentum sharing z_g
  void SetNef( const double Nef ) { mNef=Nef; };

  /// Get multi-purpose description
  std::string Tag() const { return mTag; };
  /// Set multi-purpose description
  void SetTag( const std::string tag ) { mTag=tag; };

  /// Get multi-purpose numerical value
  float Number() const { return mNumber; };
  /// Set multi-purpose numerical value
  void SetNumber(const float f)  { mNumber=f; };

protected:
  int mQuarkcharge;    ///< Charge in units of e/3
  double mZg;          ///< groomed momentum sharing z_g
  double mRg;          ///< groomed radius R_g
  double mNef;         ///< neutral energy fraction
 
  std::string mTag;    ///< Multi-purpose string
  float mNumber;       ///< Multi-purpose float
};

#endif // JETANALYSISUSERINFOS_HH
