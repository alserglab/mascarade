#ifndef MASCARADE_EFFECTIVE_H
#define MASCARADE_EFFECTIVE_H

// The effective length shared by the candidate ranking (effectiveLength()) and the continuous
// polish energy (forcePolish()): base leader length + soft viewport-overflow penalty +
// foreign-cluster arc. Defining it in one place keeps the two stages optimising exactly the
// same quantity.

#include "geometry.h"    // ClusterArcs -- the foreign-cluster arc field
#include <algorithm>

// Weight of the viewport-overflow penalty (how far a label box sticks out of the panel),
// relative to leader length. Both users of the effective length -- the discrete candidate
// ranking (effectiveLength()) and the continuous polish energy (forcePolish()) -- go through
// effectiveLengthImpl() below, so defining it here keeps the two stages in agreement.
static const double OVERFLOW_WEIGHT = 10.0;

// Effective length of a single leader.
//   arcs            prebuilt foreign-cluster arc field
//   own             0-based own cluster (skipped in the arc term)
//   ax,ay -> bx,by  leader segment (anchor -> pole)
//   base            base leader length (label centre -> pole)
//   box             padded label box (for the overflow term)
//   view            viewport bounds
static inline double effectiveLengthImpl(const ClusterArcs& arcs, int own,
                                         double ax, double ay, double bx, double by,
                                         double base, const Rect& box, const Rect& view) {
  return base + OVERFLOW_WEIGHT * box.overflow(view) + arcs.foreignArc(own, ax, ay, bx, by);
}

#endif
