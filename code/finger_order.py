# -*- coding: utf-8 -*-
"""Zinc-finger order in spectral fibres of human C2H2 zinc-finger proteins.

For every reviewed human protein annotated with C2H2-type zinc fingers and with a non-trivial 25-fibre
(N_25 > 1), fibre-mates T (same 25-spectrum, same first 24-mer) are drawn as random Eulerian trails of the
de Bruijn graph.  In S and T, C2H2 fingers are located by the motif C-X(2-4)-C-X(12)-H-X(3-5)-H, and each
finger is summarized by its specificity residues at helix positions -1, 2, 3 and 6 (7, 5, 4 and 1 residues
before the first zinc-coordinating His, which is helix position +7; there is no position 0).  The script checks that the multiset of fingers is the same in S
and T (as it must be for fingers no longer than 25 residues) and measures how their ORDER changes: the
number of array positions whose specificity residues differ.  Under any finger-additive recognition model,
the predicted binding site changes exactly when this ordered list changes.  A randomly initialized
window-25 convolutional encoder is applied to every pair for comparison.

Writes results/finger_order_output.txt, results/finger_order.json and results/finger_pairs.tsv.
"""
import os, sys, re, json, random, collections, time
import torch
import torch.nn as nn
import torch.nn.functional as Fn
import spectral_fibres as SF

torch.set_default_dtype(torch.float64)
K = 25
DRAWS = 300
FINGER = re.compile(r"(?=(C.{2,4}C.{12}(H).{3,5}H))")


def fingers(S):
    """list of (start, end, helix -1..6 as a 7-mer, specificity residues -1,2,3,6)"""
    out = []
    for m in FINGER.finditer(S):
        start = m.start(1)
        h = m.start(2)
        end = m.end(1)
        if out and start < out[-1][1]:
            continue                                    # overlapping match: keep the first
        # helix positions -1,1,...,6 (no position 0) are the 7 residues before the first His at +7,
        # e.g. Zif268 finger 1: R S D E L T R | H
        helix = S[h - 7:h]
        out.append((start, end, helix, S[h - 7] + S[h - 5] + S[h - 4] + S[h - 1]))
    return out


class CNN(nn.Module):
    def __init__(self, R, d=64):
        super().__init__()
        self.conv = nn.Conv1d(21, d, R)

    def forward(self, X):
        return Fn.gelu(self.conv(X))[0].mean(1)


IDX = {a: i for i, a in enumerate("ACDEFGHIKLMNPQRSTVWY")}


def onehot(S):
    X = torch.zeros(1, 21, len(S))
    for i, a in enumerate(S):
        X[0, IDX.get(a, 20), i] = 1.0
    return X


def main():
    t0 = time.time()
    sys.stdout = SF.Tee(os.path.join(SF.RESULTS, "finger_order_output.txt"))
    res = json.load(open(os.path.join(SF.RESULTS, "fibre_results.json")))
    seq = {acc: S for acc, _, _, S in SF.load()}
    cand = sorted([r for r in res if "C2H2-type zinc-finger" in (r["fam"] or "") and (r["logN"].get("25") or 0) > 1e-9],
                  key=lambda r: r["gene"])
    torch.manual_seed(0)
    cnn = CNN(K).eval()
    rnd = random.Random(2026)
    table, pairs = [], []
    tsv = open(os.path.join(SF.RESULTS, "finger_pairs.tsv"), "w", encoding="utf-8")
    tsv.write("gene\tmate\tpositions_changed\tfingers\tarray_positions_changed\tfingerprints_S\tfingerprints_T\n")
    print(f"=== {len(cand)} C2H2 zinc-finger proteins with N_25 > 1; {DRAWS} random trails each ===")
    for r in cand:
        S = seq[r["acc"]]
        fS = fingers(S)
        mates = []
        for _ in range(DRAWS):
            T = SF.random_trail(S, K, rnd)
            if T and T != S and T not in mates:
                mates.append(T)
        changed_orders = 0
        per = []
        for n, T in enumerate(mates):
            assert SF.spectrum(T, K) == SF.spectrum(S, K) and T[:K - 1] == S[:K - 1]
            fT = fingers(T)
            ms = collections.Counter(f[2] for f in fS) == collections.Counter(f[2] for f in fT)
            assert ms and len(fS) == len(fT), "finger multiset must be fibre-invariant"
            ordS, ordT = [f[3] for f in fS], [f[3] for f in fT]
            moved = sum(a != b for a, b in zip(ordS, ordT))
            changed_orders += moved > 0
            with torch.no_grad():
                a, b = cnn(onehot(S)), cnn(onehot(T))
            d = float((a - b).norm() / a.norm())
            assert d < 1e-12
            diff = sum(x != y for x, y in zip(S, T))
            per.append(dict(positions_changed=diff, array_positions_changed=moved))
            tsv.write(f"{r['gene']}\t{n + 1}\t{diff}\t{len(fS)}\t{moved}\t{' '.join(ordS)}\t{' '.join(ordT)}\n")
            pairs.append((r["gene"], S, T, fS, fT))
        row = dict(gene=r["gene"], acc=r["acc"], L=len(S), kappa=r["kappa"], log2N25=r["logN"]["25"],
                   fingers=len(fS), distinct_fingerprints=len({f[3] for f in fS}), mates=len(mates),
                   mates_with_new_order=changed_orders,
                   max_array_positions_changed=max([p["array_positions_changed"] for p in per], default=0),
                   max_positions_changed=max([p["positions_changed"] for p in per], default=0))
        table.append(row)
        print(f"    {r['gene']:8s} L={len(S):5d} kappa={r['kappa']:3d} log2N25={r['logN']['25']:5.1f} fingers {len(fS):3d} "
              f"(distinct specificity residues {row['distinct_fingerprints']:3d}); distinct mates {len(mates):4d}, "
              f"with a new finger order {changed_orders:4d}; max array positions changed "
              f"{row['max_array_positions_changed']:3d}, max residues changed {row['max_positions_changed']:4d}")
    tot = sum(t["mates"] for t in table)
    new = sum(t["mates_with_new_order"] for t in table)
    print(f"\n    {tot} distinct fibre-mates in {len(table)} proteins; {new} change the order of specificity residues; "
          f"finger multisets identical and window-25 CNN outputs identical in all")
    tsv.close()
    json.dump(dict(proteins=table), open(os.path.join(SF.RESULTS, "finger_order.json"), "w"), indent=1)
    # detailed case studies
    for g in ("ZNF91", "ZNF721"):
        cs = [p for p in pairs if p[0] == g]
        if not cs:
            continue
        best = max(cs, key=lambda p: sum(a[3] != b[3] for a, b in zip(p[3], p[4])))
        _, S, T, fS, fT = best
        print(f"\n=== case study {g}: {len(fS)} fingers; the mate with the most reordered fingers ===")
        for i, (a, b) in enumerate(zip(fS, fT)):
            mark = "" if a[3] == b[3] else "   <- moved"
            print(f"      finger {i + 1:2d} at {a[0] + 1:4d}-{a[1]:4d}: helix S {a[2]}  T {b[2]}   "
                  f"(-1,2,3,6) S {a[3]}  T {b[3]}{mark}")
        json.dump(dict(gene=g, S=S, T=T, fS=[list(x) for x in fS], fT=[list(x) for x in fT]),
                  open(os.path.join(SF.RESULTS, f"case_{g}.json"), "w"), indent=1)
    print(f"\n    elapsed {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
