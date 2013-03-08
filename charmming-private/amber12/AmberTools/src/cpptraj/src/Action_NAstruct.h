#ifndef INC_ACTION_NASTRUCT_H
#define INC_ACTION_NASTRUCT_H
// NAstruct
#include <vector>
#include <map>
#include "Action.h"
#include "AxisType.h"
#include "Range.h"
/// Basic Nucleic acid structure analysis. 
/** Calculate nucleic acid base/base pair structural parameters.
  * Algorithms for calculation of base/base pair structural parameters
  * adapted from:
  *   Babcock MS, Pednault EPD, Olson WK, "Nucleic Acid Structure Analysis: 
  *   Mathematics for Local Cartesian and Helical Structure Parameters That
  *   Are Truly Comparable Between Structures", J. Mol. Biol. (1994) 237,
  *   125-156.
  * NA base reference frame coordinates taken from:
  *   Olson WK, Bansal M, Burley SK, Dickerson RE, Gerstein M, Harvey SC,
  *   Heinemann U, Lu XJ, Neidle S, Shekked Z, Sklenar H, Suzuki M, Tung CS,
  *   Westhof E, Wolberger C, Berman H, "A Standard Reference Frame for the 
  *   Description of Nucleic Acid Base-pair Geometry", J. Mol. Biol. (2001)
  *   313, 229-237.
  */
class NAstruct: public Action {
    // Variables
    std::vector<AxisType> RefCoords;    ///< Hold reference frame coordinates for each base
    std::vector<AxisType> BaseAxes;     ///< Hold axis coordinates for each base
    std::vector<AxisType> BasePairAxes; ///< Hold axis coordinates for each base pair
    std::vector<AtomMask> ExpMasks;     ///< Hold atom masks of each base for input coordinates
    std::vector<AtomMask> FitMasks;     ///< Hold atom masks used to RMS-fit each base
    std::vector<int> BasePair;          ///< [Base1-0,Base2-0,IsAnti-0], [Base1-1...
    int Nbp;                            ///< Total number of base pairs
    int Nbases;                         ///< Total number of NA bases
    double HBdistCut2;                  ///< distance Cutoff^2 for determining hydrogen bonds
    double HBangleCut2;                 ///< Angle Cutoff^2 for determining if bases can h bond
    double originCut2;                  ///< Cutoff^2 for determining base-pairing vi origins
    int Nframe;                         ///< Keep track of # frames for print() function
    bool useReference_;
    std::map<std::string,AxisType::NAbaseType> CustomMap;
    std::map<std::string,AxisType::NAbaseType>::iterator customRes;
    CpptrajFile BPOut;
    CpptrajFile BPstepOut;
    CpptrajFile HelixOut;
 /*   // Datasets
    DataSetList SHEAR;
    DataSetList STRETCH;
    DataSetList STAGGER;
    DataSetList BUCKLE;
    DataSetList PROPELLER;
    DataSetList OPENING;
    DataSetList SHIFT;
    DataSetList SLIDE;
    DataSetList RISE;
    DataSetList TILT;
    DataSetList ROLL;
    DataSetList TWIST;*/
    // Init Args
    Range resRange;
    //char *outFilename;
    //char *naoutFilename;
    bool noheader;
    // Functions
    void ClearLists();
    bool GCpair(AxisType *, AxisType *);
    bool ATpair(AxisType *, AxisType *);
    bool basesArePaired(AxisType *, AxisType *);
    int determineBasePairing();
    int setupBaseAxes(Frame *);
    int calculateParameters(AxisType &, AxisType &, AxisType*, double*);
    int helicalParameters(AxisType &, AxisType &, double *);
    int determineBaseParameters();
    int determineBasepairParameters();
#   ifdef NASTRUCTDEBUG
    // DEBUG - Class to hold PDB output
    class AxisPDBwriter {
        int pdbatom;
        CpptrajFile pdbfile;
      public:
        AxisPDBwriter() { pdbatom=0; }
        ~AxisPDBwriter() { pdbfile.CloseFile(); }
        void Open(const char *fname) {
          pdbfile.SetupFile((char*)fname, WRITE, 0);
          pdbfile.OpenFile();
        }
        void Write(AxisType &axis, int res, char *resname) {
          axis.WritePDB(pdbfile, res, resname, &pdbatom);
        }
        void WriteAxes(AxisType &axis, int res, char *resname) {
          axis.WriteAxesPDB(pdbfile, res, resname, &pdbatom);
        }
    };
    // DEBUG - used to trigger AxisPDBwriter for first call of calculateParameters
    bool calcparam;
#   endif
  public:
    NAstruct();
    ~NAstruct();

    int init();
    int setup();
    int action();
    void print();
};
#endif
