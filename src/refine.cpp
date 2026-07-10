#include "geometry.h"
#include <vector>
#include <algorithm>
using namespace Rcpp;

// Per-candidate box + leader geometry shared by the two refine kernels: the pairwise conflict
// predicates over candidate indices and the per-label conflict counts against a current
// selection. (`len` is not held here: oneMoveSweep ranks by raw length, twoMoveBnB by squared.)
struct CandidateSet {
  NumericVector cxmin, cxmax, cymin, cymax;   // padded box extents, one per candidate
  NumericVector ex, ey, tx, ty;               // leader per candidate: anchor -> pole
  // per-candidate leader bounding box, cached by prepare() for the broad-phase rejects below
  std::vector<double> lxmin, lxmax, lymin, lymax;

  // cache each leader's axis-aligned bounding box (call once after construction). The pairwise
  // predicates are dominated by the Boost segment tests, so a cheap AABB reject in front of them
  // skips the (common) far-apart pairs without a full segment/segment or segment/box test.
  void prepare() {
    int n = ex.size();
    lxmin.resize(n);
    lxmax.resize(n);
    lymin.resize(n);
    lymax.resize(n);
    for (int i = 0; i < n; ++i) {
      lxmin[i] = std::min(ex[i], tx[i]);
      lxmax[i] = std::max(ex[i], tx[i]);
      lymin[i] = std::min(ey[i], ty[i]);
      lymax[i] = std::max(ey[i], ty[i]);
    }
  }
  // do the padded boxes of candidates a and b overlap?
  bool boxesOverlap(int a, int b) const {
    return cxmin[a] < cxmax[b] && cxmin[b] < cxmax[a]
        && cymin[a] < cymax[b] && cymin[b] < cymax[a];
  }
  // does leader `s` (its cached AABB) come near box `b`? (broad-phase for leaderBoxConflict)
  bool leaderNearBox(int s, int b) const {
    return lxmin[s] <= cxmax[b] && cxmin[b] <= lxmax[s]
        && lymin[s] <= cymax[b] && cymin[b] <= lymax[s];
  }
  // does either candidate's leader pass through the other's box?
  bool leaderBoxConflict(int a, int b) const {
    if (leaderNearBox(a, b)
        && segbox(ex[a], ey[a], tx[a], ty[a], cxmin[b], cxmax[b], cymin[b], cymax[b])) {
      return true;
    }
    if (leaderNearBox(b, a)
        && segbox(ex[b], ey[b], tx[b], ty[b], cxmin[a], cxmax[a], cymin[a], cymax[a])) {
      return true;
    }
    return false;
  }
  // do the two leaders cross?
  bool leadersCross(int a, int b) const {
    if (lxmin[a] > lxmax[b] || lxmax[a] < lxmin[b]
        || lymin[a] > lymax[b] || lymax[a] < lymin[b]) {
      return false;
    }
    return segcross(ex[a], ey[a], tx[a], ty[a], ex[b], ey[b], tx[b], ty[b]);
  }
  // any of the three conflict types between candidates a and b?
  bool anyConflict(int a, int b) const {
    return boxesOverlap(a, b) || leaderBoxConflict(a, b) || leadersCross(a, b);
  }
  // conflict counts (box-box, leader-leader, leader-box) of `cand` placed for `label`, versus
  // the current placements `cur`, excluding `label` and `exclude` (pass -1 for no exclusion).
  void countConflicts(int cand, int label, int exclude, const std::vector<int>& cur,
                      int& boxBox, int& leaderLeader, int& leaderBox) const {
    boxBox = 0;
    leaderLeader = 0;
    leaderBox = 0;
    int K = (int) cur.size();
    for (int m = 0; m < K; ++m) {
      if (m == label || m == exclude) {
        continue;
      }
      if (boxesOverlap(cand, cur[m])) {
        ++boxBox;
      }
      if (leadersCross(cand, cur[m])) {
        ++leaderLeader;
      }
      if (leaderBoxConflict(cand, cur[m])) {
        ++leaderBox;
      }
    }
  }
};

//' Single-move conflict sweep
//'
//' Multi-pass Gauss-Seidel refinement. Each pass reorders labels by their current conflict
//' vector (box-box, leader-leader, leader-box, then leader length) descending, then moves each
//' to its lexicographically best candidate versus the current others (including staying put).
//' Stops when a pass makes no change.
//'
//' @param cxmin,cxmax,cymin,cymax Numeric padded box extents, one per candidate.
//' @param ex,ey Numeric leader start (anchor) per candidate.
//' @param tx,ty Numeric leader target (pole) per candidate.
//' @param len Numeric ranking length per candidate.
//' @param rows List of integer vectors: the 0-indexed candidate rows available to each label.
//' @param init Integer vector: the starting candidate index per label.
//' @param maxpass Integer cap on the number of sweeps.
//' @return Integer vector: the chosen candidate index per label.
//' @keywords internal
// [[Rcpp::export]]
IntegerVector oneMoveSweep(NumericVector cxmin, NumericVector cxmax,
                           NumericVector cymin, NumericVector cymax,
                           NumericVector ex, NumericVector ey,
                           NumericVector tx, NumericVector ty,
                           NumericVector len, List rows, IntegerVector init, int maxpass) {
  int K = rows.size();
  std::vector<std::vector<int> > candidatesFor(K);
  for (int k = 0; k < K; ++k) {
    IntegerVector v = rows[k];
    candidatesFor[k].assign(v.begin(), v.end());
  }
  std::vector<int> cur(init.begin(), init.end());
  CandidateSet cs{cxmin, cxmax, cymin, cymax, ex, ey, tx, ty};
  cs.prepare();

  int applied = 0;
  int passes = 0;
  std::vector<int> order(K), curBox(K), curCross(K), curLeaderBox(K);
  for (int pass = 0; pass < maxpass; ++pass) {
    ++passes;
    for (int label = 0; label < K; ++label) {
      cs.countConflicts(cur[label], label, -1, cur, curBox[label], curCross[label], curLeaderBox[label]);
    }
    for (int i = 0; i < K; ++i) {
      order[i] = i;
    }
    // visit the most-conflicted (then longest-leader) labels first
    std::sort(order.begin(), order.end(), [&](int a, int b) {
      if (curBox[a] != curBox[b]) {
        return curBox[a] > curBox[b];
      }
      if (curCross[a] != curCross[b]) {
        return curCross[a] > curCross[b];
      }
      if (curLeaderBox[a] != curLeaderBox[b]) {
        return curLeaderBox[a] > curLeaderBox[b];
      }
      return len[cur[a]] > len[cur[b]];
    });

    bool changed = false;
    for (int oi = 0; oi < K; ++oi) {
      int label = order[oi];
      int bestBox = 1 << 30;
      int bestCross = 1 << 30;
      int bestLeaderBox = 1 << 30;
      double bestLen = 1e18;
      int bestCand = cur[label];
      for (int cand : candidatesFor[label]) {
        int boxBox, leaderLeader, leaderBox;
        cs.countConflicts(cand, label, -1, cur, boxBox, leaderLeader, leaderBox);
        double candLen = len[cand];
        bool better;
        if (boxBox != bestBox) {
          better = boxBox < bestBox;
        } else if (leaderLeader != bestCross) {
          better = leaderLeader < bestCross;
        } else if (leaderBox != bestLeaderBox) {
          better = leaderBox < bestLeaderBox;
        } else {
          better = candLen < bestLen - 1e-9;
        }
        if (better) {
          bestBox = boxBox;
          bestCross = leaderLeader;
          bestLeaderBox = leaderBox;
          bestLen = candLen;
          bestCand = cand;
        }
      }
      if (bestCand != cur[label]) {
        cur[label] = bestCand;
        changed = true;
        ++applied;
      }
    }
    if (!changed) {
      break;
    }
  }

  IntegerVector out(K);
  for (int k = 0; k < K; ++k) {
    out[k] = cur[k];
  }
  out.attr("passes") = passes;
  out.attr("applied") = applied;
  return out;
}

//' Two-move conflict / length refinement
//'
//' Per label (heaviest leader first), picks the lexicographically best move (by change in
//' box-box, leader-leader, leader-box conflicts, then change in length) and applies it if it
//' improves. A CONFLICT-FREE label uses a length branch-and-bound: candidates are
//' length-ascending and a shorter c1 needs a partner only when exactly one other label
//' currently sits in that slot (that label is the partner -- it need not itself be in
//' conflict), with the length bound pruning the length-increasing tail -- the exact per-step
//' optimum. A CONFLICTED label drops the pruning and searches ALL pairs of candidates (this
//' label's against every other label's), driving out conflicts present in the input; only
//' conflicted labels pay this cost.
//'
//' @param cxmin,cxmax,cymin,cymax Numeric padded box extents, one per candidate.
//' @param ex,ey Numeric leader start (anchor) per candidate.
//' @param tx,ty Numeric leader target (pole) per candidate.
//' @param len Numeric ranking length per candidate.
//' @param rows List of integer vectors: the 0-indexed candidate rows available to each label.
//' @param init Integer vector: the starting candidate index per label.
//' @param maxpass Integer cap on the number of passes.
//' @param sq Logical; if `TRUE` the length objective uses the squared length.
//' @return Integer vector: the chosen candidate index per label.
//' @keywords internal
// [[Rcpp::export]]
IntegerVector twoMoveBnB(NumericVector cxmin, NumericVector cxmax,
                         NumericVector cymin, NumericVector cymax,
                         NumericVector ex, NumericVector ey,
                         NumericVector tx, NumericVector ty,
                         NumericVector len, List rows, IntegerVector init, int maxpass, bool sq) {
  int K = rows.size();
  std::vector<std::vector<int> > candidatesFor(K);
  for (int k = 0; k < K; ++k) {
    IntegerVector v = rows[k];
    candidatesFor[k].assign(v.begin(), v.end());
  }
  std::vector<int> cur(init.begin(), init.end());
  CandidateSet cs{cxmin, cxmax, cymin, cymax, ex, ey, tx, ty};
  cs.prepare();
  // ranking length: squared (sq) or raw
  auto lengthCost = [&](int cand) -> double {
    return sq ? len[cand] * len[cand] : len[cand];
  };

  int passes = 0;
  std::vector<int> order(K);
  for (int pass = 0; pass < maxpass; ++pass) {
    ++passes;
    bool changed = false;
    for (int i = 0; i < K; ++i) {
      order[i] = i;
    }
    std::sort(order.begin(), order.end(), [&](int a, int b) {
      return len[cur[a]] > len[cur[b]];
    });
    // smallest length reduction available across all labels (for the conflict-free bound)
    double minPartner = 1e18;
    for (int m = 0; m < K; ++m) {
      double d = lengthCost(candidatesFor[m][0]) - lengthCost(cur[m]);
      if (d < minPartner) {
        minPartner = d;
      }
    }

    for (int oi = 0; oi < K; ++oi) {
      int label1 = order[oi];
      double curLen1 = lengthCost(cur[label1]);
      bool haveMove = false;
      int bestDBox = 0;
      int bestDCross = 0;
      int bestDLeaderBox = 0;
      double bestDLen = 0;
      int bestL2 = -1;
      int bestC1 = -1;
      int bestC2 = -1;

      // is (dBox, dCross, dLeaderBox, dLen) lexicographically better than the best so far?
      // With no move yet, "better" means an improvement over doing nothing (all deltas < 0).
      auto isBetter = [&](int dBox, int dCross, int dLeaderBox, double dLen) -> bool {
        if (!haveMove) {
          if (dBox != 0) {
            return dBox < 0;
          }
          if (dCross != 0) {
            return dCross < 0;
          }
          if (dLeaderBox != 0) {
            return dLeaderBox < 0;
          }
          return dLen < -1e-9;
        }
        if (dBox != bestDBox) {
          return dBox < bestDBox;
        }
        if (dCross != bestDCross) {
          return dCross < bestDCross;
        }
        if (dLeaderBox != bestDLeaderBox) {
          return dLeaderBox < bestDLeaderBox;
        }
        return dLen < bestDLen - 1e-9;
      };
      auto consider = [&](int dBox, int dCross, int dLeaderBox, double dLen,
                          int l2, int c1, int c2) {
        if (isBetter(dBox, dCross, dLeaderBox, dLen)) {
          haveMove = true;
          bestDBox = dBox;
          bestDCross = dCross;
          bestDLeaderBox = dLeaderBox;
          bestDLen = dLen;
          bestL2 = l2;
          bestC1 = c1;
          bestC2 = c2;
        }
      };

      bool label1Conflicted = false;
      for (int m = 0; m < K; ++m) {
        if (m != label1 && cs.anyConflict(cur[label1], cur[m])) {
          label1Conflicted = true;
          break;
        }
      }

      if (label1Conflicted) {
        // exhaustive: every c1 for label1, paired with every c2 of every other label
        int curBox1, curCross1, curLeaderBox1;
        cs.countConflicts(cur[label1], label1, -1, cur, curBox1, curCross1, curLeaderBox1);
        for (int c1 : candidatesFor[label1]) {
          // single move: label1 -> c1
          int box1, cross1, leaderBox1;
          cs.countConflicts(c1, label1, -1, cur, box1, cross1, leaderBox1);
          consider(box1 - curBox1, cross1 - curCross1, leaderBox1 - curLeaderBox1,
                   lengthCost(c1) - curLen1, -1, c1, -1);
          // two-move: label1 -> c1 together with label2 -> c2
          for (int label2 = 0; label2 < K; ++label2) {
            if (label2 == label1) {
              continue;
            }
            int old1Box, old1Cross, old1LeaderBox;
            cs.countConflicts(cur[label1], label1, label2, cur, old1Box, old1Cross, old1LeaderBox);
            int old2Box, old2Cross, old2LeaderBox;
            cs.countConflicts(cur[label2], label2, label1, cur, old2Box, old2Cross, old2LeaderBox);
            // the current label1<->label2 conflict is excluded from old* above; the new c1<->c2
            // conflict is 0 (conflicting pairs are skipped), so subtract the current one.
            int mutualBox = cs.boxesOverlap(cur[label1], cur[label2]) ? 1 : 0;
            int mutualCross = cs.leadersCross(cur[label1], cur[label2]) ? 1 : 0;
            int mutualLeaderBox = cs.leaderBoxConflict(cur[label1], cur[label2]) ? 1 : 0;
            int new1Box, new1Cross, new1LeaderBox;
            cs.countConflicts(c1, label1, label2, cur, new1Box, new1Cross, new1LeaderBox);
            double curLen2 = lengthCost(cur[label2]);
            for (int c2 : candidatesFor[label2]) {
              if (cs.anyConflict(c1, c2)) {
                continue;
              }
              int new2Box, new2Cross, new2LeaderBox;
              cs.countConflicts(c2, label2, label1, cur, new2Box, new2Cross, new2LeaderBox);
              int dBox = (new1Box + new2Box) - (old1Box + old2Box) - mutualBox;
              int dCross = (new1Cross + new2Cross) - (old1Cross + old2Cross) - mutualCross;
              int dLeaderBox = (new1LeaderBox + new2LeaderBox)
                             - (old1LeaderBox + old2LeaderBox) - mutualLeaderBox;
              double dLen = (lengthCost(c1) - curLen1) + (lengthCost(c2) - curLen2);
              consider(dBox, dCross, dLeaderBox, dLen, label2, c1, c2);
            }
          }
        }
      } else {
        // conflict-free: length branch-and-bound with the single-blocker partner rule
        for (int c1 : candidatesFor[label1]) {
          double d1 = lengthCost(c1) - curLen1;
          bool lengthPruned = haveMove ? (d1 + minPartner >= bestDLen) : (d1 + minPartner >= 0);
          if (lengthPruned) {
            break;
          }
          int blockers = 0;
          int blocker = -1;
          for (int m = 0; m < K; ++m) {
            if (m == label1) {
              continue;
            }
            if (cs.anyConflict(c1, cur[m])) {
              ++blockers;
              if (blockers > 1) {
                break;
              }
              blocker = m;
            }
          }
          if (blockers == 0) {
            // single move: label1 -> c1 (c1 is conflict-free here)
            int box1, cross1, leaderBox1;
            cs.countConflicts(c1, label1, -1, cur, box1, cross1, leaderBox1);
            int curBox1, curCross1, curLeaderBox1;
            cs.countConflicts(cur[label1], label1, -1, cur, curBox1, curCross1, curLeaderBox1);
            consider(box1 - curBox1, cross1 - curCross1, leaderBox1 - curLeaderBox1,
                     d1, -1, c1, -1);
            continue;
          }
          if (blockers > 1) {
            continue;
          }
          int label2 = blocker;
          double curLen2 = lengthCost(cur[label2]);
          int old1Box, old1Cross, old1LeaderBox;
          cs.countConflicts(cur[label1], label1, label2, cur, old1Box, old1Cross, old1LeaderBox);
          int old2Box, old2Cross, old2LeaderBox;
          cs.countConflicts(cur[label2], label2, label1, cur, old2Box, old2Cross, old2LeaderBox);
          for (int c2 : candidatesFor[label2]) {
            double dLen = d1 + (lengthCost(c2) - curLen2);
            bool bound = haveMove ? (dLen >= bestDLen) : (dLen >= 0);
            if (bound) {
              break;
            }
            if (cs.anyConflict(c1, c2)) {
              continue;
            }
            int new1Box, new1Cross, new1LeaderBox;
            cs.countConflicts(c1, label1, label2, cur, new1Box, new1Cross, new1LeaderBox);
            int new2Box, new2Cross, new2LeaderBox;
            cs.countConflicts(c2, label2, label1, cur, new2Box, new2Cross, new2LeaderBox);
            int dBox = (new1Box + new2Box) - (old1Box + old2Box);
            int dCross = (new1Cross + new2Cross) - (old1Cross + old2Cross);
            int dLeaderBox = (new1LeaderBox + new2LeaderBox) - (old1LeaderBox + old2LeaderBox);
            consider(dBox, dCross, dLeaderBox, dLen, label2, c1, c2);
          }
        }
      }

      if (haveMove
          && (bestDBox < 0 || bestDCross < 0 || bestDLeaderBox < 0 || bestDLen < -1e-9)) {
        cur[label1] = bestC1;
        if (bestL2 >= 0) {
          cur[bestL2] = bestC2;
        }
        changed = true;
      }
    }
    if (!changed) {
      break;
    }
  }

  IntegerVector out(K);
  for (int k = 0; k < K; ++k) {
    out[k] = cur[k];
  }
  out.attr("passes") = passes;
  return out;
}
