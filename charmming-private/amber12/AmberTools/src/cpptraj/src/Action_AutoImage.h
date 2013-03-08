#ifndef INC_ACTION_AUTOIMAGE_H
#define INC_ACTION_AUTOIMAGE_H
#include "Action.h"
class Action_AutoImage : public Action {
  public:
    Action_AutoImage();

  private:
    int init();
    int setup();
    int action();

    AtomMask anchorMask_; ///< Used to center anchor region.
    char* anchor_;  ///< Mask expression for anchor region.
    char* fixed_;   ///< Mask expression for fixed region.
    char* mobile_;  ///< Mask expression for mobile region.

    bool origin_;         ///< If true imaging occurs w.r.t. coordinate origin.
    bool ortho_;          ///< If true imaging is orthogonal.
    bool center_;         ///< If true imaging of mobile region uses molecule center.
    bool truncoct_;
    enum TriclinicArg {OFF, FORCE, FAMILIAR};
    TriclinicArg triclinic_; ///< Determine whether triclinic code should be used.

    typedef std::vector<int> pairList;
    pairList anchorList_;
    pairList fixedList_;
    pairList mobileList_;

    pairList SetupAtomRanges(const char*);
};
#endif
