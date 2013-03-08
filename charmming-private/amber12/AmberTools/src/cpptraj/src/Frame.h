#ifndef INC_FRAME_H
#define INC_FRAME_H
#include "AtomMask.h"
// Class: Frame
/// Hold coordinates, perform various operations/transformations on them.
/** Intended to hold coordinates e.g. from a trajectory or reference frame,
  * along with box coordinates (used in imaging calculations) and optionally 
  * with mass information and/or velocity information. Mass is stored since 
  * several functions (like COM, RADGYR etc) have the option to factor in 
  * the mass of the atoms involved, and this avoids having to pass a mass 
  * pointer in, which takes the burden of keeping track of mass away from 
  * actions etc. Mass is stored when the frame is initially created, and is 
  * modified if necessary by SetFrameFromMask (which is the case when e.g. 
  * calculating RMSD).
  */
class Frame {
  protected:
    static const size_t COORDSIZE;
    static const size_t BOXSIZE;
  public:
    double *X;     ///< Coord array, X0 Y0 Z0 X1 Y1 Z1 ...
    int natom;     ///< Number of atoms
    int maxnatom;  ///< Number of atoms for which space has been allocated
    int N;         ///< Number of coords, natom*3
    double box[6]; ///< Box coords, 3xlengths, 3xangles
    double T;      ///< Temperature
    double *V;     ///< Velocities
    double *Mass;  ///< Mass

    Frame();
    virtual ~Frame();             // Destructor is virtual since this class can be inherited
    Frame & operator=(const Frame&);
    Frame & operator+=(const Frame&);
    Frame & operator-=(const Frame&);
    int Divide(Frame&, double); 
    Frame(const Frame&);

    int SetupFrame(int, double*);
    int SetupFrameV(int,double*,bool);

    int SetupFrameFromMask(AtomMask *, double *);
    int SetupFrameFromCoords(float *, int);
    int SetupFrameFromCoords(float *);
    Frame *FrameCopy();
    //int Resize(int,bool,bool);

    // Coordinate manipulation
    void ZeroCoords();
    //void AddCoord(Frame*);
    void Divide(double);
    void Translate(double *);
    void Translate(double *, int,int);
    void Trans_Rot_Trans(double *, double *);
    void Rotate(double *);
    void InverseRotate(double *);
    void Center(AtomMask&, bool,bool);
    void CenterReference(double *, bool);
    void ShiftToGeometricCenter();
    void SetupImageTruncoct(double*, AtomMask*,bool,bool);
    void ImageNonortho(bool, double*, double*, double*, bool, bool, bool, std::vector<int> &);
    void ImageNonortho(double*, double*, bool, bool, double*, double*, double*);
    int SetupImageOrtho(double*, double*, bool);
    void ImageOrtho(double*,double*, bool, bool, std::vector<int> &);
    void ImageOrtho(double*, double*, double*, double*);
    // Coordinate assignment/extraction
    void printAtomCoord(int);
    void GetCoord(double *, int);
    void SetCoord(int, double *);
    double *Coord(int);
    void SetFrameFromMask(Frame*, AtomMask *);
    int SetFrameCoordsFromMask(double *, AtomMask *);
    int SetFrameCoords(double *);
    // Center of mass
    double CenterOfMass(AtomMask*, double *);
    double GeometricCenter(AtomMask*, double *);
    double CenterOfMass(double*,int,int);
    double GeometricCenter(double*,int,int);
    // Coordinate calculation
    double BoxToRecip(double *, double *);
    double DIST2(AtomMask*, AtomMask*, bool, int, double *, double *);
    double DIST2(int, int, int, double *, double *);
    double DIST2(double*, int, int, double *, double *);
    double DIST(int, int);
    double DIST2(int, int);
    double COORDDIST(int, int);
    double COORDDIST2(int, int);
    double ANGLE(AtomMask*, AtomMask*, AtomMask*,bool);
    double ANGLE(int, int, int);
    double DIHEDRAL(AtomMask *, AtomMask *, AtomMask *, AtomMask *,bool);
    double DIHEDRAL(int,int,int,int);
    double PUCKER(AtomMask*,AtomMask*,AtomMask*,AtomMask*,AtomMask*,int,bool,bool);
    double RADGYR(AtomMask *, bool, double *);
    double RMSD(Frame*, double*, double*,bool);
    double RMSD_CenteredRef( Frame &, double[9], double[6], bool);
    double RMSD(Frame*,bool);
    double DISTRMSD( Frame * );

    void SetAxisOfRotation(double *, int, int);
    void RotateAroundAxis(double *, double, AtomMask &);
};
#endif
