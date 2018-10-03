/** @file PpZgParameters.hh
    @author Kolja Kauder
    @brief Common parameters
    @details Used to quickly include the same parameters into different macros.
    @date Mar 23, 2017
   
 */

#ifndef EICENUMS_HH
#define EICENUMS_HH

/// method to use for x, Q^2, Breit frame, ...
enum KINEMATICSMETHOD {
  ELECTRON=0,
  JACQUETBLONDEL=1,
  DOUBLEANGLE=2,
  TRUEKINEMATICS=3
};

/// Input file format
enum INTYPE{ MCTREE, INTREE, EPTREE, HERWIGTREE };

/// Background Subtraction method
enum BGTYPE{ NONE=0, AREA=1, CONSTSUBPRE=2, CONSTSUBPOST=3 };

/// Return values for the main routine
enum class EVENTRESULT{
  PROBLEM,
  ENDOFINPUT,
  NOTACCEPTED,
  NOCONSTS,
  NOJETS,
  JETSFOUND
};

#endif // EICENUMS_HH
