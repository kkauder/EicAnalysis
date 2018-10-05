#ifndef JETANALYSISUSERINFOS_HH
#define JETANALYSISUSERINFOS_HH

/** @author Kauder:Kolja
    @version Revision 0.1
    @brief Some UserInfo UserInfo derivatives
    @date Mar 03, 2015
*/


// =============================================================================

/** Only one can be attached per PseudoJet.
    Splitting into constituent info and jet info (to save space).
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
  int GetQuarkCharge() const { return mQuarkcharge; };
  /// Set charge in units of e/3
  void SetQuarkCharge(const int Quarkcharge )  { mQuarkcharge = Quarkcharge; };

  /// Get PDG code, quick reference at http://home.fnal.gov/~mrenna/lutp0613man2/node44.html
  int GetPdg() const { return mPdg; };
  /// Set PDG code
  void SetPdg(const int pdg)  { mPdg=pdg; };

  /// Get multi-purpose description
  std::string GetTag() const { return mTag; };
  /// Set multi-purpose description
  void SetTag( const std::string tag ) { mTag=tag; };

  /// Get multi-purpose numerical value
  float GetNumber() const { return mNumber; };
  /// Set multi-purpose numerical value
  void SetNumber(const float f)  { mNumber=f; };

protected:
  int mQuarkcharge;    ///< Charge in units of e/3
  int mPdg;            ///< PDG code

  std::string mTag;    ///< Multi-purpose
  float mNumber;       ///< Multi-purpose
};

#endif // JETANALYSISUSERINFOS_HH
