#ifndef MASCARADE_CONSTANTS_H
#define MASCARADE_CONSTANTS_H

// Weight of the viewport-overflow penalty (how far a label box sticks out of the panel),
// relative to leader length. Shared by the discrete effective length (effective.cpp) and the
// polish energy (polish.cpp) so the two stages agree.
static const double OVERFLOW_WEIGHT = 10.0;

#endif
