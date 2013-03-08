#ifndef INC_TRAJ_AMBERRESTART_H
#define INC_TRAJ_AMBERRESTART_H
#include "TrajectoryIO.h"
// Class: AmberRestart.h
/// Reads and writes formatted (ASCII text) amber
class AmberRestart : public TrajectoryIO {
    int restartAtoms;     ///< Number of atoms in restart file
    int natom3;           ///< Number of coords
    int frameSize;        ///< Size of 1 coord frame in bytes, inc box & velo if present
    size_t coordSize_;    ///< Size of 1 coord frame in bytes, used for blank reads.
    char *frameBuffer;    ///< Used to read in restart coord
    int numBoxCoords;     ///< Number of box coords (3 or 6)
    double restartTime;   ///< Time in restart file, read in
    double restartTemp;   ///< (Optional) replica temperature, read in.
    double time0;         ///< For writes, restart time offset
    double dt;            ///< For writes, restart timestep (scaling)
    bool singleWrite;     ///< If false, frame # will be appended to output filename

    // Inherited functions
    int setupRead(AmberParm*);
    int setupWrite(AmberParm*,int);
    int openTraj();
    void closeTraj();
    int readFrame(int,double*,double*,double*,double*);
    int writeFrame(int,double*,double*,double*,double);
    int processWriteArgs(ArgList*);
    void info();

    int getBoxAngles(char *, int);

  public:

    AmberRestart();
    ~AmberRestart();
    // AmberRestart-specific functions
    void SetNoVelocity();
};
#endif
