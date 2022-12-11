#ifndef _INCL_FAST_PROPAGATE_
#define _INCL_FAST_PROPAGATE_

#include "base.h"
#include "atomic_data.h"
#include "spec.h"


int propagate_lc_los(int numCells, cellProp cell[], outArrayBlock *specBlockIn, double zuni, int numOut, double tOutList[], double bkg[][6], char outputDir[], char prefix[], int trimSpec, int outCellGap);

#endif
