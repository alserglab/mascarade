#include "geometry.h"
#include <vector>
#include <algorithm>
using namespace Rcpp;

// Multi-pass single-move sweep (Gauss-Seidel). Each pass: reorder labels by their
// current conflict vector (bb, ll, lb, leader-length) descending, then visit each once
// and move it to its lexicographically best candidate (bb, ll, lb, len) vs the current
// others, including staying. Stop when a pass makes no change.
// [[Rcpp::export]]
IntegerVector oneMoveSweep(NumericVector cxmin, NumericVector cxmax, NumericVector cymin, NumericVector cymax,
                           NumericVector ex, NumericVector ey, NumericVector tx, NumericVector ty,
                           NumericVector len, List rows, IntegerVector init, int maxpass) {
  int K = rows.size();
  std::vector<std::vector<int> > byl(K);
  for (int k = 0; k < K; ++k) { IntegerVector v = rows[k]; byl[k].assign(v.begin(), v.end()); }
  std::vector<int> cur(init.begin(), init.end());
  auto obb = [&](int a, int b)->bool { return cxmin[a] < cxmax[b] && cxmin[b] < cxmax[a] && cymin[a] < cymax[b] && cymin[b] < cymax[a]; };
  auto olb = [&](int a, int b)->bool { return segbox(ex[a], ey[a], tx[a], ty[a], cxmin[b], cxmax[b], cymin[b], cymax[b])
                                            || segbox(ex[b], ey[b], tx[b], ty[b], cxmin[a], cxmax[a], cymin[a], cymax[a]); };
  auto oll = [&](int a, int b)->bool { return segcross(ex[a], ey[a], tx[a], ty[a], ex[b], ey[b], tx[b], ty[b]); };
  auto ev = [&](int c, int L, int& b, int& x, int& l) { b = 0; x = 0; l = 0;
    for (int M = 0; M < K; ++M) { if (M == L) continue; if (obb(c, cur[M])) ++b; if (oll(c, cur[M])) ++x; if (olb(c, cur[M])) ++l; } };
  int applied = 0, passes = 0; std::vector<int> ord(K), cbb(K), cll(K), clb(K);
  for (int pass = 0; pass < maxpass; ++pass) { ++passes;
    for (int L = 0; L < K; ++L) { int b, x, l; ev(cur[L], L, b, x, l); cbb[L] = b; cll[L] = x; clb[L] = l; }
    for (int i = 0; i < K; ++i) ord[i] = i;
    std::sort(ord.begin(), ord.end(), [&](int a, int b) {
      if (cbb[a] != cbb[b]) return cbb[a] > cbb[b];
      if (cll[a] != cll[b]) return cll[a] > cll[b];
      if (clb[a] != clb[b]) return clb[a] > clb[b];
      return len[cur[a]] > len[cur[b]]; });
    bool changed = false;
    for (int oi = 0; oi < K; ++oi) { int L = ord[oi];
      int bb_ = 1 << 30, ll_ = 1 << 30, lb_ = 1 << 30; double le_ = 1e18; int bc = cur[L];
      for (int c : byl[L]) { int b, x, l; ev(c, L, b, x, l); double le = len[c];
        if (b < bb_ || (b == bb_ && (x < ll_ || (x == ll_ && (l < lb_ || (l == lb_ && le < le_ - 1e-9)))))) { bb_ = b; ll_ = x; lb_ = l; le_ = le; bc = c; } }
      if (bc != cur[L]) { cur[L] = bc; changed = true; ++applied; }
    }
    if (!changed) break;
  }
  IntegerVector out(K); for (int k = 0; k < K; ++k) out[k] = cur[k];
  out.attr("passes") = passes; out.attr("applied") = applied;
  return out;
}

// Two-move refinement, per label (heaviest leader first), picking the lexicographically best
// (dbb, dll, dlb, dlen) move among those it visits and applying it if it is an improvement.
//   * CONFLICT-FREE label: length branch-and-bound. Candidates are length-ascending; a shorter
//     c1 needs a partner only if blocked by exactly one label (that blocker is the partner:
//     0 blockers = a single move, 2+ = unfixable by one partner), and the length bound prunes
//     the length-increasing tail. At a conflict-free start this is the exact per-step optimum.
//   * CONFLICTED label: drop the length pruning and search ALL pairs of candidates (this
//     label's c1 x every other label's c2), so conflicts present in the INPUT are driven out.
//     Only conflicted labels pay this cost, so it stays far cheaper than an all-pairs search
//     over every label, and it reaches a conflict-free layout quickly.
// [[Rcpp::export]]
IntegerVector twoMoveBnB(NumericVector cxmin, NumericVector cxmax, NumericVector cymin, NumericVector cymax,
                         NumericVector ex, NumericVector ey, NumericVector tx, NumericVector ty,
                         NumericVector len, List rows, IntegerVector init, int maxpass, bool sq) {
  int K = rows.size();
  std::vector<std::vector<int> > byl(K);
  for (int k = 0; k < K; ++k) {
    IntegerVector v = rows[k];
    byl[k].assign(v.begin(), v.end());
  }
  std::vector<int> cur(init.begin(), init.end());
  auto bbf = [&](int a, int b)->bool { return cxmin[a] < cxmax[b] && cxmin[b] < cxmax[a] && cymin[a] < cymax[b] && cymin[b] < cymax[a]; };
  auto lbf = [&](int a, int b)->bool { return segbox(ex[a], ey[a], tx[a], ty[a], cxmin[b], cxmax[b], cymin[b], cymax[b])
                                            || segbox(ex[b], ey[b], tx[b], ty[b], cxmin[a], cxmax[a], cymin[a], cymax[a]); };
  auto llf = [&](int a, int b)->bool { return segcross(ex[a], ey[a], tx[a], ty[a], ex[b], ey[b], tx[b], ty[b]); };
  auto conf = [&](int a, int b)->bool { return bbf(a, b) || lbf(a, b) || llf(a, b); };
  auto lv = [&](int c)->double { return sq ? len[c] * len[c] : len[c]; };
  // conflicts of candidate c (placed for label L) with the current others, excluding L and exl
  auto cnt3 = [&](int c, int L, int exl, int& b, int& x, int& l) {
    b = 0; x = 0; l = 0;
    for (int M = 0; M < K; ++M) {
      if (M == L || M == exl) { continue; }
      if (bbf(c, cur[M])) { ++b; }
      if (llf(c, cur[M])) { ++x; }
      if (lbf(c, cur[M])) { ++l; }
    }
  };

  int passes = 0;
  std::vector<int> ord(K);
  for (int pass = 0; pass < maxpass; ++pass) {
    ++passes;
    bool changed = false;
    for (int i = 0; i < K; ++i) {
      ord[i] = i;
    }
    std::sort(ord.begin(), ord.end(), [&](int a, int b) { return len[cur[a]] > len[cur[b]]; });
    // smallest length reduction available across all labels (for the conflict-free length bound)
    double minPart = 1e18;
    for (int M = 0; M < K; ++M) {
      double d = lv(byl[M][0]) - lv(cur[M]);
      if (d < minPart) { minPart = d; }
    }

    for (int oi = 0; oi < K; ++oi) {
      int L1 = ord[oi];
      double curL1 = lv(cur[L1]);
      bool have = false;
      int bBB = 0, bXX = 0, bLB = 0; double bDL = 0;
      int bL2 = -1, bc1 = -1, bc2 = -1;
      auto better = [&](int db, int dx, int dl, double dlen)->bool {
        if (!have) { return db < 0 || (db == 0 && (dx < 0 || (dx == 0 && (dl < 0 || (dl == 0 && dlen < -1e-9))))); }
        if (db != bBB) { return db < bBB; }
        if (dx != bXX) { return dx < bXX; }
        if (dl != bLB) { return dl < bLB; }
        return dlen < bDL - 1e-9;
      };
      auto consider = [&](int db, int dx, int dl, double dlen, int L2, int c1, int c2) {
        if (better(db, dx, dl, dlen)) {
          have = true; bBB = db; bXX = dx; bLB = dl; bDL = dlen; bL2 = L2; bc1 = c1; bc2 = c2;
        }
      };

      // is this label currently in conflict with anyone?
      bool l1Conflicted = false;
      for (int M = 0; M < K; ++M) {
        if (M != L1 && conf(cur[L1], cur[M])) { l1Conflicted = true; break; }
      }

      if (l1Conflicted) {
        // exhaustive: every c1 for L1, paired with every c2 of every other label (no pruning)
        int ob1s, ox1s, ol1s;
        cnt3(cur[L1], L1, -1, ob1s, ox1s, ol1s);
        for (int c1 : byl[L1]) {
          int b1, x1, l1;
          cnt3(c1, L1, -1, b1, x1, l1);
          consider(b1 - ob1s, x1 - ox1s, l1 - ol1s, lv(c1) - curL1, -1, c1, -1);   // single move
          for (int L2 = 0; L2 < K; ++L2) {
            if (L2 == L1) { continue; }
            int ob1, ox1, ol1; cnt3(cur[L1], L1, L2, ob1, ox1, ol1);
            int ob2, ox2, ol2; cnt3(cur[L2], L2, L1, ob2, ox2, ol2);
            // the current L1<->L2 conflict leaves ob1/ob2 (cnt3 excludes each other); the new
            // c1<->c2 conflict is 0 (pairs that conflict are skipped), so subtract it explicitly
            int ccb = bbf(cur[L1], cur[L2]) ? 1 : 0;
            int ccx = llf(cur[L1], cur[L2]) ? 1 : 0;
            int ccl = lbf(cur[L1], cur[L2]) ? 1 : 0;
            int c1b, c1x, c1l; cnt3(c1, L1, L2, c1b, c1x, c1l);
            double curL2 = lv(cur[L2]);
            for (int c2 : byl[L2]) {
              if (conf(c1, c2)) { continue; }
              int c2b, c2x, c2l; cnt3(c2, L2, L1, c2b, c2x, c2l);
              int db = (c1b + c2b) - (ob1 + ob2) - ccb;
              int dx = (c1x + c2x) - (ox1 + ox2) - ccx;
              int dl = (c1l + c2l) - (ol1 + ol2) - ccl;
              double dlen = (lv(c1) - curL1) + (lv(c2) - curL2);
              consider(db, dx, dl, dlen, L2, c1, c2);
            }
          }
        }
      } else {
        // conflict-free: length branch-and-bound with the single-blocker partner rule
        for (int c1 : byl[L1]) {
          double d1 = lv(c1) - curL1;
          bool boundBeaten = have ? (d1 + minPart >= bDL) : (d1 + minPart >= 0);
          if (boundBeaten) { break; }
          int nb = 0, blk = -1;
          for (int M = 0; M < K; ++M) {
            if (M == L1) { continue; }
            if (conf(c1, cur[M])) { if (++nb > 1) { break; } blk = M; }
          }
          if (nb == 0) {
            int b1, x1, l1; cnt3(c1, L1, -1, b1, x1, l1);
            int ob1, ox1, ol1; cnt3(cur[L1], L1, -1, ob1, ox1, ol1);
            consider(b1 - ob1, x1 - ox1, l1 - ol1, d1, -1, c1, -1);
            continue;
          }
          if (nb > 1) { continue; }
          int L2 = blk;
          double curL2 = lv(cur[L2]);
          int ob1, ox1, ol1; cnt3(cur[L1], L1, L2, ob1, ox1, ol1);
          int ob2, ox2, ol2; cnt3(cur[L2], L2, L1, ob2, ox2, ol2);
          for (int c2 : byl[L2]) {
            double dd = d1 + (lv(c2) - curL2);
            bool bnd = have ? (dd >= bDL) : (dd >= 0);
            if (bnd) { break; }
            if (conf(c1, c2)) { continue; }
            int c1b, c1x, c1l; cnt3(c1, L1, L2, c1b, c1x, c1l);
            int c2b, c2x, c2l; cnt3(c2, L2, L1, c2b, c2x, c2l);
            int db = (c1b + c2b) - (ob1 + ob2), dx = (c1x + c2x) - (ox1 + ox2), dl = (c1l + c2l) - (ol1 + ol2);
            consider(db, dx, dl, dd, L2, c1, c2);
          }
        }
      }

      if (have && (bBB < 0 || bXX < 0 || bLB < 0 || bDL < -1e-9)) {
        cur[L1] = bc1;
        if (bL2 >= 0) { cur[bL2] = bc2; }
        changed = true;
      }
    }
    if (!changed) { break; }
  }
  IntegerVector out(K);
  for (int k = 0; k < K; ++k) {
    out[k] = cur[k];
  }
  out.attr("passes") = passes;
  return out;
}
