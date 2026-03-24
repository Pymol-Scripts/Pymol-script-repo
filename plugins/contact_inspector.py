"""
Contact Inspector - PyMOL Plugin
================================
Maestro-inspired protein-ligand interaction viewer for PyMOL. Automatically
detects and visualizes all major non-covalent interactions, with ligand
stepping for docking pose review.

H-bonds use PyMOL's built-in polar contact detection (cmd.distance mode=2)
which correctly handles donor/acceptor chemistry. All other interaction
types are detected geometrically.

Interaction categories:
  Non-covalent bonds:  H-bonds, halogen bonds, salt bridges, aromatic H-bonds
  Pi interactions:     Pi-pi stacking (face-to-face & edge-to-face), pi-cation
  Contacts/Clashes:    Good, bad, ugly  (off by default)

Installation:
  1. Plugin > Plugin Manager > Install New Plugin > choose this file, or
  2. run /path/to/contact_inspector.py   then   ci_gui

Authors: Evert J. Homan, PhD; Claude (Anthropic)
Date:    2026-03-24
Version: 1.0
License: MIT
"""

from __future__ import annotations

import math
from collections import defaultdict
from typing import List, Dict, Optional, Set

try:
    import numpy as np
except ImportError:
    np = None

from pymol import cmd, CmdException

# ---------------------------------------------------------------------------
# Colors
# ---------------------------------------------------------------------------

COLORS = {
    "ci_hbond":       (1.00, 0.85, 0.00),
    "ci_halogen":     (0.60, 0.20, 0.90),
    "ci_salt":        (0.90, 0.20, 0.90),
    "ci_arom_hb":     (0.30, 0.85, 0.50),
    "ci_pipi":        (0.30, 0.75, 1.00),
    "ci_pi_cat":      (0.20, 0.80, 0.20),
    "ci_clash_good":  (0.20, 0.80, 0.20),
    "ci_clash_bad":   (1.00, 0.60, 0.00),
    "ci_clash_ugly":  (1.00, 0.15, 0.15),
}

def _register_colors():
    """Register custom colors with PyMOL. Called at runtime to ensure PyMOL is ready."""
    for name, rgb in COLORS.items():
        cmd.set_color(name, list(rgb))

# ---------------------------------------------------------------------------
# Thresholds
# ---------------------------------------------------------------------------

VDW_RADII = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
    "P": 1.80, "S": 1.80, "Cl": 1.75, "Br": 1.85, "I": 1.98,
    "Fe": 1.80, "Zn": 1.39, "Mg": 1.73, "Ca": 1.74, "Mn": 1.80,
    "Se": 1.90, "Si": 2.10, "B": 1.92,
}

HBOND_ELEMENTS = {"N", "O", "S", "F"}
# H-bonds handled entirely by PyMOL's cmd.distance(mode=2)

HALOGEN_DONORS = {"Cl", "Br", "I"}
HALOGEN_ACCEPTORS = {"O", "N", "S"}
HALOGEN_DIST_MAX = 3.5

AROM_HBOND_DIST_MAX = 3.8
AROM_HBOND_ACCEPTORS = {"O", "N", "S"}

SALT_BRIDGE_DIST_MAX = 4.0

PIPI_FTF_DIST_MAX = 4.8
PIPI_FTF_ANGLE_MAX = 40.0
PIPI_ETF_DIST_MAX = 5.5
PIPI_ETF_ANGLE_MIN = 45.0
PIPI_ETF_ANGLE_MAX = 90.0

PI_CATION_DIST_MAX = 6.0

CLASH_GOOD_FRAC = 1.30        # comfortable VDW contact (up to 130% of VDW sum)
CLASH_BAD_FRAC = 0.89         # mild steric overlap
CLASH_UGLY_FRAC = 0.75        # severe steric overlap

SHELL_DIST = 5.0

DASH_RADIUS = 0.06
DASH_GAP = 0.35
DASH_LENGTH = 0.20
LABEL_SIZE = 14

# ---------------------------------------------------------------------------
# Object tracking
# ---------------------------------------------------------------------------

_created_objects: Set[str] = set()

_OBJ_PTS = "_ci_pts"
_OBJ_SHELL = "shell"
_OBJ_SURF  = "_ci_surf"

_INTERACTION_NAMES = {
    "hbonds": "hbonds",
    "halogen": "halogen_bonds",
    "salt": "salt_bridges",
    "arom_hb": "arom_hbonds",
    "pipi": "pi_pi",
    "pi_cation": "pi_cation",
    "clash_good": "clash_good",
    "clash_bad": "clash_bad",
    "clash_ugly": "clash_ugly",
}

def _track(name):
    _created_objects.add(name)

def _clear_all():
    for name in list(_created_objects):
        try: cmd.delete(name)
        except Exception: pass
    _created_objects.clear()

# ---------------------------------------------------------------------------
# Vector helpers
# ---------------------------------------------------------------------------

def _norm(v):
    if np is not None:
        n = np.linalg.norm(v)
        return v / n if n > 1e-9 else v
    mag = math.sqrt(sum(x * x for x in v))
    return [x / mag for x in v] if mag > 1e-9 else v

def _dot(a, b):
    if np is not None:
        return float(np.dot(a, b))
    return sum(x * y for x, y in zip(a, b))

def _dist(a, b):
    if np is not None:
        return float(np.linalg.norm(np.array(a) - np.array(b)))
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))

def _centroid(pts):
    if np is not None:
        return np.mean(pts, axis=0)
    n = len(pts)
    return [sum(p[i] for p in pts) / n for i in range(3)]

def _angle_normals(n1, n2):
    c = abs(_dot(_norm(n1), _norm(n2)))
    return math.degrees(math.acos(min(1.0, max(0.0, c))))

# ---------------------------------------------------------------------------
# Ring detection
# ---------------------------------------------------------------------------

def _ring_planarity_rmsd(coords):
    """Compute RMSD of ring atoms from their best-fit plane.
    Aromatic rings: < 0.1 A.  Cyclohexyl chair: ~ 0.5 A."""
    if np is not None:
        pts = np.array(coords)
        cen = pts.mean(axis=0)
        pts_c = pts - cen
        # SVD to find best-fit plane normal
        _, _, vh = np.linalg.svd(pts_c)
        normal = vh[2]  # smallest singular value = plane normal
        dists = pts_c @ normal
        return float(np.sqrt(np.mean(dists ** 2)))
    else:
        # Pure-python fallback: use ring normal from Newell's method
        n = len(coords)
        cen = [sum(c[i] for c in coords) / n for i in range(3)]
        nm = [0.0, 0.0, 0.0]
        for i in range(n):
            j = (i + 1) % n
            ci, cj = coords[i], coords[j]
            nm[0] += (ci[1] - cj[1]) * (ci[2] + cj[2])
            nm[1] += (ci[2] - cj[2]) * (ci[0] + cj[0])
            nm[2] += (ci[0] - cj[0]) * (ci[1] + cj[1])
        nm = _norm(nm)
        # Distance of each point from plane through centroid
        ss = 0.0
        for c in coords:
            d = sum((c[i] - cen[i]) * nm[i] for i in range(3))
            ss += d * d
        return math.sqrt(ss / n)

def _find_aromatic_rings(model):
    n = len(model.atom)
    adj: Dict[int, Set[int]] = defaultdict(set)
    for bond in model.bond:
        i, j = bond.index
        adj[i].add(j); adj[j].add(i)

    ring_elems = {"C", "N", "O", "S"}
    rings: List[List[int]] = []
    seen: Set[frozenset] = set()

    for start in range(n):
        if model.atom[start].symbol not in ring_elems:
            continue
        _dfs_rings(adj, start, start, [start], set(), rings, seen, model, 6)

    out = []
    for ring in rings:
        if len(ring) not in (5, 6):
            continue
        ok = True
        for idx in ring:
            a = model.atom[idx]
            if a.symbol not in ring_elems:
                ok = False; break
            # sp2 carbon: exactly 3 total bonds (2 ring + 1 substituent/H)
            # sp3 carbon: 4 total bonds → not aromatic
            if a.symbol == "C":
                total_nb = len(adj[idx])
                if total_nb > 3:
                    ok = False; break
        if not ok:
            continue
        # Planarity check: aromatic rings are flat (RMSD < 0.15 A)
        coords = [model.atom[i].coord for i in ring]
        if _ring_planarity_rmsd(coords) > 0.15:
            continue
        out.append(ring)
    return out, adj


def _dfs_rings(adj, start, cur, path, vis, rings, seen, model, mx):
    if len(path) > mx:
        return
    for nb in adj[cur]:
        if nb == start and len(path) >= 5:
            fs = frozenset(path)
            if fs not in seen:
                seen.add(fs); rings.append(list(path))
        elif nb not in vis and nb > start:
            vis.add(nb); path.append(nb)
            _dfs_rings(adj, start, nb, path, vis, rings, seen, model, mx)
            path.pop(); vis.discard(nb)


def _ring_normal(coords):
    n = len(coords)
    nm = [0.0, 0.0, 0.0]
    for i in range(n):
        j = (i + 1) % n
        ci, cj = coords[i], coords[j]
        nm[0] += (ci[1] - cj[1]) * (ci[2] + cj[2])
        nm[1] += (ci[2] - cj[2]) * (ci[0] + cj[0])
        nm[2] += (ci[0] - cj[0]) * (ci[1] + cj[1])
    return _norm(nm)


def _get_rings_info(model):
    rings, adj = _find_aromatic_rings(model)
    out = []
    for r in rings:
        coords = [model.atom[i].coord for i in r]
        out.append((r, _centroid(coords), _ring_normal(coords)))
    return out, adj


def _get_aromatic_ch_atoms(model, rings, adj):
    ring_atoms = set()
    for r in rings:
        ring_atoms.update(r)
    ch_atoms = []
    for idx in ring_atoms:
        a = model.atom[idx]
        if a.symbol != "C":
            continue
        for nb in adj[idx]:
            if model.atom[nb].symbol == "H":
                ch_atoms.append((a, tuple(a.coord), tuple(model.atom[nb].coord)))
                break
    return ch_atoms


# ---------------------------------------------------------------------------
# Interaction result (for non-hbond types detected geometrically)
# ---------------------------------------------------------------------------

class InteractionResult:
    def __init__(self):
        self.hbond_count: int = 0   # set by visualize() from cmd.distance return value
        self.halogen: List[dict] = []
        self.salt_bridges: List[dict] = []
        self.arom_hbonds: List[dict] = []
        self.pipi: List[dict] = []
        self.pi_cation: List[dict] = []
        self.clash_good: List[dict] = []
        self.clash_bad: List[dict] = []
        self.clash_ugly: List[dict] = []


# ---------------------------------------------------------------------------
# Detection (non-hbond interactions only; H-bonds use PyMOL mode=2)
# ---------------------------------------------------------------------------

def detect_interactions(
    lig_sel, prot_sel, state=-1,
    do_halogen=True, do_salt=True, do_arom_hb=True,
    do_pipi=True, do_pi_cation=True,
    do_clash_good=False, do_clash_bad=False, do_clash_ugly=False,
) -> InteractionResult:

    result = InteractionResult()
    search = max(HALOGEN_DIST_MAX, SALT_BRIDGE_DIST_MAX,
                 PI_CATION_DIST_MAX, PIPI_ETF_DIST_MAX) + 2.0

    lig_model = cmd.get_model(lig_sel, state=state)
    prot_model = cmd.get_model(
        f"{prot_sel} within {search} of ({lig_sel})", state=state)
    if not lig_model.atom or not prot_model.atom:
        return result

    lig_atoms = [(a, tuple(a.coord)) for a in lig_model.atom]
    prot_atoms = [(a, tuple(a.coord)) for a in prot_model.atom]

    do_any_clash = do_clash_good or do_clash_bad or do_clash_ugly

    def _el(a): return a.symbol.strip().capitalize()
    def _il(a): return f"{a.resn} {a.name}"
    def _ip(a): return f"{a.chain}/{a.resn}{a.resi}.{a.name}"

    # --- Pairwise ---
    for la, lc in lig_atoms:
        le = _el(la)
        for pa, pc in prot_atoms:
            pe = _el(pa)
            d = _dist(lc, pc)

            if do_halogen and d <= HALOGEN_DIST_MAX:
                if ((le in HALOGEN_DONORS and pe in HALOGEN_ACCEPTORS) or
                    (pe in HALOGEN_DONORS and le in HALOGEN_ACCEPTORS)):
                    result.halogen.append({
                        "p1": lc, "p2": pc, "dist": d,
                        "info1": _il(la), "info2": _ip(pa)})

            if do_salt and d <= SALT_BRIDGE_DIST_MAX:
                lpos = (le == "N" and la.formal_charge > 0)
                lneg = (le == "O" and la.formal_charge < 0)
                ppos = (pa.resn in ("ARG","LYS","HIS","HID","HIE","HIP")
                        and pa.name in ("NH1","NH2","NE","NZ","ND1","NE2"))
                pneg = ((pa.resn == "ASP" and pa.name in ("OD1","OD2"))
                        or (pa.resn == "GLU" and pa.name in ("OE1","OE2")))
                if (lpos and pneg) or (lneg and ppos):
                    result.salt_bridges.append({
                        "p1": lc, "p2": pc, "dist": d,
                        "info1": _il(la), "info2": _ip(pa)})

            if do_any_clash and le != "H" and pe != "H":
                vdw = VDW_RADII.get(le, 1.70) + VDW_RADII.get(pe, 1.70)
                max_d = vdw * CLASH_GOOD_FRAC
                if d <= max_d:
                    frac = d / vdw if vdw > 0 else 1.0
                    # Skip polar pairs for bad/ugly (those are H-bonds)
                    # but allow them for good contacts
                    polar = (le in HBOND_ELEMENTS and pe in HBOND_ELEMENTS)
                    if frac < CLASH_UGLY_FRAC:
                        if do_clash_ugly and not polar:
                            result.clash_ugly.append({
                                "p1": lc, "p2": pc, "dist": d, "vdw": vdw,
                                "quality": "ugly",
                                "info1": _il(la), "info2": _ip(pa)})
                    elif frac < CLASH_BAD_FRAC:
                        if do_clash_bad and not polar:
                            result.clash_bad.append({
                                "p1": lc, "p2": pc, "dist": d, "vdw": vdw,
                                "quality": "bad",
                                "info1": _il(la), "info2": _ip(pa)})
                    elif frac < CLASH_GOOD_FRAC:
                        if do_clash_good:
                            result.clash_good.append({
                                "p1": lc, "p2": pc, "dist": d, "vdw": vdw,
                                "quality": "good",
                                "info1": _il(la), "info2": _ip(pa)})

    # --- Ring-based ---
    lig_rings_info = []
    prot_rings_info = []
    _prm = None
    _prm_adj = None
    _lm = None
    _lm_adj = None

    need_rings = do_pipi or do_arom_hb or do_pi_cation
    if need_rings:
        try:
            _lm = cmd.get_model(lig_sel, state=state)
            lig_rings_info, _lm_adj = _get_rings_info(_lm)
        except Exception:
            pass
        try:
            prs = (f"({prot_sel} within {PI_CATION_DIST_MAX+3} of ({lig_sel}))"
                   f" and (resn PHE+TYR+TRP+HIS+HIE+HID+HIP)")
            _prm = cmd.get_model(prs, state=state)
            prot_rings_info, _prm_adj = _get_rings_info(_prm)
        except Exception:
            pass

    # Pi-pi
    if do_pipi and lig_rings_info and prot_rings_info:
        for _, lrc, lrn in lig_rings_info:
            for pri, prc, prn in prot_rings_info:
                d = _dist(lrc, prc)
                ang = _angle_normals(lrn, prn)
                a0 = _prm.atom[pri[0]]
                info = f"{a0.chain}/{a0.resn}{a0.resi}"
                if d <= PIPI_FTF_DIST_MAX and ang <= PIPI_FTF_ANGLE_MAX:
                    result.pipi.append({
                        "p1": list(lrc), "p2": list(prc), "dist": d,
                        "angle": ang, "type": "face-to-face", "info2": info})
                elif (d <= PIPI_ETF_DIST_MAX and
                      PIPI_ETF_ANGLE_MIN <= ang <= PIPI_ETF_ANGLE_MAX):
                    result.pipi.append({
                        "p1": list(lrc), "p2": list(prc), "dist": d,
                        "angle": ang, "type": "edge-to-face", "info2": info})

    # Aromatic H-bonds: aromatic C-H ... acceptor (O/N/S)
    if do_arom_hb:
        if _lm and _lm_adj and lig_rings_info:
            lig_ring_indices = [r for r, _, _ in lig_rings_info]
            lig_ch = _get_aromatic_ch_atoms(_lm, lig_ring_indices, _lm_adj)
            for ca, cc, hc in lig_ch:
                for pa, pc in prot_atoms:
                    if _el(pa) not in AROM_HBOND_ACCEPTORS:
                        continue
                    d = _dist(cc, pc)
                    if d <= AROM_HBOND_DIST_MAX:
                        vh_c = [cc[i] - hc[i] for i in range(3)]
                        vh_a = [pc[i] - hc[i] for i in range(3)]
                        dot = _dot(_norm(vh_c), _norm(vh_a))
                        angle = math.degrees(math.acos(
                            min(1.0, max(-1.0, dot))))
                        if angle > 120.0:
                            result.arom_hbonds.append({
                                "p1": hc, "p2": pc, "dist": _dist(hc, pc),
                                "info1": f"{ca.resn} {ca.name}-H",
                                "info2": _ip(pa)})

        if _prm and _prm_adj and prot_rings_info:
            prot_ring_indices = [r for r, _, _ in prot_rings_info]
            prot_ch = _get_aromatic_ch_atoms(_prm, prot_ring_indices, _prm_adj)
            for ca, cc, hc in prot_ch:
                for la, lc in lig_atoms:
                    if _el(la) not in AROM_HBOND_ACCEPTORS:
                        continue
                    d = _dist(cc, lc)
                    if d <= AROM_HBOND_DIST_MAX:
                        vh_c = [cc[i] - hc[i] for i in range(3)]
                        vh_a = [lc[i] - hc[i] for i in range(3)]
                        dot = _dot(_norm(vh_c), _norm(vh_a))
                        angle = math.degrees(math.acos(
                            min(1.0, max(-1.0, dot))))
                        if angle > 120.0:
                            result.arom_hbonds.append({
                                "p1": hc, "p2": lc, "dist": _dist(hc, lc),
                                "info1": f"{ca.chain}/{ca.resn}{ca.resi}.{ca.name}-H",
                                "info2": _il(la)})

    # Pi-cation
    if do_pi_cation:
        pcs = (f"({prot_sel} within {PI_CATION_DIST_MAX+1} of ({lig_sel}))"
               f" and ((resn ARG and name CZ) or (resn LYS and name NZ))")
        try:
            pcm = cmd.get_model(pcs, state=state)
            for pa in pcm.atom:
                pac = tuple(pa.coord)
                for _, lrc, _ in lig_rings_info:
                    d = _dist(pac, lrc)
                    if d <= PI_CATION_DIST_MAX:
                        result.pi_cation.append({
                            "p1": list(lrc), "p2": pac, "dist": d,
                            "info1": "lig ring", "info2": _ip(pa)})
        except Exception:
            pass
        if _prm:
            for la, lc in lig_atoms:
                if la.formal_charge > 0:
                    for pri, prc, _ in prot_rings_info:
                        d = _dist(lc, prc)
                        if d <= PI_CATION_DIST_MAX:
                            a0 = _prm.atom[pri[0]]
                            result.pi_cation.append({
                                "p1": lc, "p2": list(prc), "dist": d,
                                "info1": _il(la),
                                "info2": f"{a0.chain}/{a0.resn}{a0.resi} ring"})

    return result


# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------

def _clear_contacts():
    keep = {_OBJ_SHELL, _OBJ_SURF}
    for name in list(_created_objects):
        if name not in keep:
            try: cmd.delete(name)
            except Exception: pass
            _created_objects.discard(name)


def _style(name, color, radius=DASH_RADIUS, gap=DASH_GAP, length=DASH_LENGTH):
    cmd.set("dash_radius", radius, name)
    cmd.set("dash_gap", gap, name)
    cmd.set("dash_length", length, name)
    cmd.color(color, name)
    cmd.set("label_color", "white", name)
    cmd.set("label_size", LABEL_SIZE, name)


def _add_pair(pts, pid, p1, p2, dist_name):
    r = str(pid)
    cmd.pseudoatom(pts, pos=list(p1), resi=r, name="L", chain="X")
    cmd.pseudoatom(pts, pos=list(p2), resi=r, name="P", chain="X")
    cmd.distance(dist_name,
                 f"{pts} and resi {r} and name L",
                 f"{pts} and resi {r} and name P")


def visualize(lig_sel, prot_sel, result: InteractionResult,
              show_hbonds=True, show_labels=True, state=-1):
    """Visualize all interactions.

    H-bonds are created via PyMOL's built-in polar contact detection
    (cmd.distance mode=2). All other types use pseudoatom pairs.
    """
    _register_colors()
    _clear_contacts()

    # --- H-bonds via PyMOL polar contacts (mode=2) ---
    # cmd.distance returns the number of contacts found (or -1 on error).
    # count_atoms() always returns 0 for distance/measurement objects and
    # cannot be used to check whether any H-bonds were found.
    hb_name = _INTERACTION_NAMES["hbonds"]
    if show_hbonds:
        try:
            n_hb = cmd.distance(hb_name, lig_sel, prot_sel, mode=2)
            if n_hb is not None and n_hb > 0:
                result.hbond_count = int(n_hb)
                _track(hb_name)
                _style(hb_name, "ci_hbond")
            else:
                result.hbond_count = 0
                try: cmd.delete(hb_name)
                except Exception: pass
        except Exception:
            result.hbond_count = 0
            try: cmd.delete(hb_name)
            except Exception: pass

    # --- All other types via pseudoatom pairs ---
    pts = _OBJ_PTS
    _track(pts)
    pid = 0

    def _draw(items, obj_name, color, **kw):
        nonlocal pid
        if not items:
            return
        _track(obj_name)
        for it in items:
            _add_pair(pts, pid, it["p1"], it["p2"], obj_name)
            pid += 1
        _style(obj_name, color, **kw)

    N = _INTERACTION_NAMES
    _draw(result.halogen,      N["halogen"],  "ci_halogen")
    _draw(result.salt_bridges, N["salt"],     "ci_salt")
    _draw(result.arom_hbonds,  N["arom_hb"],  "ci_arom_hb")

    _draw(result.pipi,         N["pipi"],     "ci_pipi",
          gap=0.30, length=0.25, radius=0.05)
    _draw(result.pi_cation,    N["pi_cation"],"ci_pi_cat",
          gap=0.30, length=0.25, radius=0.05)

    _draw(result.clash_good,   N["clash_good"],  "ci_clash_good",
          gap=0.15, length=0.10)
    _draw(result.clash_bad,    N["clash_bad"],   "ci_clash_bad",
          gap=0.15, length=0.10)
    _draw(result.clash_ugly,   N["clash_ugly"],  "ci_clash_ugly",
          gap=0.15, length=0.10)

    # Contacts/clashes: never show distance labels (too cluttered)
    for q in ("clash_good", "clash_bad", "clash_ugly"):
        if N[q] in _created_objects:
            cmd.hide("labels", N[q])

    if pid > 0:
        cmd.hide("everything", pts)

    if not show_labels:
        for n in N.values():
            if n in _created_objects:
                cmd.hide("labels", n)


# ---------------------------------------------------------------------------
# Scene preparation
# ---------------------------------------------------------------------------

_ELEM_COLORS = [
    ("N",  "tv_blue"), ("O",  "tv_red"),   ("S",  "tv_yellow"),
    ("P",  "orange"),  ("F",  "palegreen"),
    ("Cl", "green"),   ("Br", "firebrick"), ("I",  "purple"),
]

def _apply_elem_colors(sel):
    """Apply standard element colors to heteroatoms in sel."""
    for elem, color in _ELEM_COLORS:
        try:
            cmd.color(color, f"({sel}) and elem {elem}")
        except Exception:
            pass

def _color_rainbow_elem(sel):
    """Color carbon atoms rainbow by residue number; heteroatoms by element; H white."""
    cmd.spectrum("count", "rainbow", f"({sel}) and elem C", byres=1)
    _apply_elem_colors(sel)
    cmd.color("white", f"({sel}) and elem H")


def _prepare_scene(protein_sel, ligand_sels):
    """Color protein rainbow (C atoms) + element colors, cartoon, hide non-polar H."""
    try:
        cmd.show("cartoon", protein_sel)
        _color_rainbow_elem(protein_sel)

        # Hide non-polar H on protein (keep polar H on N/O/S visible)
        cmd.hide("everything",
                 f"({protein_sel}) and elem H and "
                 f"not (neighbor (elem N+O+S))")

        # Show ligands as sticks, white carbons + element colors for heteroatoms
        lig_u = _lig_union(ligand_sels)
        cmd.show("sticks", lig_u)
        cmd.hide("everything",
                 f"({lig_u}) and elem H and "
                 f"not (neighbor (elem N+O+S))")
        cmd.color("white", f"({lig_u}) and elem C")
        cmd.color("white", f"({lig_u}) and elem H")
        _apply_elem_colors(lig_u)
    except Exception as e:
        print(f"Contact Inspector: scene setup warning: {e}")


def _lig_union(ligand_sels):
    """Build a PyMOL union selection string from ligand(s)."""
    if isinstance(ligand_sels, str):
        return ligand_sels
    return " or ".join(f"({l})" for l in ligand_sels)


# ---------------------------------------------------------------------------
# Residue shell
# ---------------------------------------------------------------------------

def _create_shell(protein_sel, ligand_sels, dist=SHELL_DIST):
    lig_union = _lig_union(ligand_sels)

    if _OBJ_SHELL in _created_objects:
        try: cmd.delete(_OBJ_SHELL)
        except Exception: pass
        _created_objects.discard(_OBJ_SHELL)

    shell_sel = f"byres (({protein_sel}) within {dist} of ({lig_union}))"

    try:
        cmd.create(_OBJ_SHELL, shell_sel)
    except CmdException:
        return
    _track(_OBJ_SHELL)

    cmd.hide("everything", _OBJ_SHELL)
    cmd.show("lines", _OBJ_SHELL)
    cmd.remove(f"{_OBJ_SHELL} and elem H and not (neighbor (elem N+O+S))")

    # Rainbow C atoms + element colors for heteroatoms (matching protein)
    _color_rainbow_elem(_OBJ_SHELL)

    # Transparent light-grey surface on the ATOM-BASED 5 Å shell (no byres expansion)
    if _OBJ_SURF in _created_objects:
        try: cmd.delete(_OBJ_SURF)
        except Exception: pass
        _created_objects.discard(_OBJ_SURF)

    atom_surf_sel = f"({protein_sel}) within {dist} of ({lig_union})"
    try:
        cmd.create(_OBJ_SURF, atom_surf_sel)
        _track(_OBJ_SURF)
        cmd.hide("everything", _OBJ_SURF)
        cmd.show("surface", _OBJ_SURF)
        cmd.set("surface_color", "grey80", _OBJ_SURF)
        cmd.set("transparency", 0.4, _OBJ_SURF)
    except Exception:
        pass

    # Label CA atoms with residue name + number
    cmd.label(f"{_OBJ_SHELL} and name CA",
              '"%s %s" % (resn, resi)')


# ---------------------------------------------------------------------------
# Ligand stepper
# ---------------------------------------------------------------------------

class LigandStepper:
    def __init__(self):
        self.protein_sel = ""
        self.ligand_objects: List[str] = []
        self.current_index = 0
        self.mode = "objects"
        self.state_object = ""
        self.last_result: Optional[InteractionResult] = None

        self.show_hbonds = True
        self.show_halogen = True
        self.show_salt = True
        self.show_arom_hb = True
        self.show_pipi = True
        self.show_pi_cation = True
        self.show_clash_good = False
        self.show_clash_bad = False
        self.show_clash_ugly = False
        self.show_labels = True

    def setup_objects(self, prot, ligs):
        self.protein_sel = prot
        self.ligand_objects = ligs
        self.current_index = 0
        self.mode = "objects"
        _prepare_scene(prot, ligs if isinstance(ligs, str) else ligs)
        _create_shell(prot, ligs)
        if ligs:
            self._show_current()

    def setup_states(self, prot, obj):
        self.protein_sel = prot
        self.state_object = obj
        self.current_index = 0
        self.mode = "states"
        _prepare_scene(prot, obj)
        _create_shell(prot, obj)
        self._show_current()

    def _count(self):
        if self.mode == "objects":
            return len(self.ligand_objects)
        return cmd.count_states(self.state_object)

    def _label(self):
        if self.mode == "objects":
            return (self.ligand_objects[self.current_index]
                    if self.ligand_objects else "none")
        return f"{self.state_object} state {self.current_index + 1}"

    def next(self):
        c = self._count()
        if c:
            self.current_index = (self.current_index + 1) % c
            self._show_current()

    def prev(self):
        c = self._count()
        if c:
            self.current_index = (self.current_index - 1) % c
            self._show_current()

    def goto(self, i):
        if 0 <= i < self._count():
            self.current_index = i
            self._show_current()

    def _show_current(self):
        if self.mode == "objects":
            obj_names = set(cmd.get_names("objects"))
            for n in self.ligand_objects:
                if n in obj_names:
                    cmd.disable(n)
            if self.ligand_objects:
                cur = self.ligand_objects[self.current_index]
                if cur in obj_names:
                    cmd.enable(cur)
                cmd.zoom(f"({self.protein_sel}) within 8 of ({cur})",
                         buffer=3.0, animate=1)
                self._update(cur)
        else:
            st = self.current_index + 1
            cmd.set("state", st, self.state_object)
            cmd.zoom(f"({self.protein_sel}) within 8 of "
                     f"({self.state_object})", buffer=3.0, animate=1)
            self._update(self.state_object, state=st)

    def _update(self, lig, state=-1):
        import traceback
        try:
            r = detect_interactions(
                lig, self.protein_sel, state=state,
                do_halogen=self.show_halogen,
                do_salt=self.show_salt, do_arom_hb=self.show_arom_hb,
                do_pipi=self.show_pipi, do_pi_cation=self.show_pi_cation,
                do_clash_good=self.show_clash_good,
                do_clash_bad=self.show_clash_bad,
                do_clash_ugly=self.show_clash_ugly)
            self.last_result = r
            visualize(lig, self.protein_sel, r,
                      show_hbonds=self.show_hbonds,
                      show_labels=self.show_labels, state=state)
        except Exception as e:
            msg = f"Contact Inspector error: {e}\n{traceback.format_exc()}"
            print(msg)
            if _gui_window is not None:
                try:
                    from pymol.Qt import QtWidgets
                    QtWidgets.QMessageBox.warning(_gui_window, "Contact Inspector", str(e))
                except Exception:
                    pass

    def summary(self):
        r = self.last_result
        if r is None:
            return "No interactions detected."
        lines = [f"=== {self._label()} "
                 f"({self.current_index + 1}/{self._count()}) ==="]

        # H-bond count stored by visualize() from cmd.distance return value
        if r.hbond_count > 0:
            lines.append(f"\nH-bonds: {r.hbond_count} (PyMOL polar contacts)")

        def _s(title, items, extra_fn=None):
            if not items: return
            lines.append(f"\n{title} ({len(items)}):")
            for it in items:
                ex = extra_fn(it) if extra_fn else ""
                i1 = it.get("info1", "")
                i2 = it.get("info2", "")
                lines.append(f"  {i1} -- {i2}  {it['dist']:.2f} A{ex}")

        _s("Halogen bonds", r.halogen)
        _s("Salt bridges", r.salt_bridges)
        _s("Aromatic H-bonds", r.arom_hbonds)
        _s("Pi-pi", r.pipi,
           lambda x: f"  {x['angle']:.0f} deg [{x['type']}]")
        _s("Pi-cation", r.pi_cation)
        _s("Clash good", r.clash_good)
        _s("Clash bad", r.clash_bad)
        _s("Clash ugly", r.clash_ugly)

        total = sum(len(getattr(r, a)) for a in (
            "halogen","salt_bridges","arom_hbonds",
            "pipi","pi_cation","clash_good","clash_bad","clash_ugly"))
        if total == 0 and r.hbond_count == 0:
            lines.append("  No interactions found.")
        return "\n".join(lines)


_stepper = LigandStepper()


# ---------------------------------------------------------------------------
# Keybindings
# ---------------------------------------------------------------------------

def _bind_keys():
    cmd.set_key("right", ci_next)
    cmd.set_key("left", ci_prev)
    print("Contact Inspector: LEFT/RIGHT arrow keys bound.")

def _unbind_keys():
    try:
        cmd.set_key("right", lambda: None)
        cmd.set_key("left", lambda: None)
    except Exception:
        pass


# ---------------------------------------------------------------------------
# CLI commands
# ---------------------------------------------------------------------------

def ci_setup(protein="polymer.protein", ligands="organic", mode="auto"):
    """
DESCRIPTION
    Setup Contact Inspector.

USAGE
    ci_setup [protein [, ligands [, mode]]]

EXAMPLES
    ci_setup
    ci_setup protein=chain A, ligands=LIG1,LIG2,LIG3
    ci_setup protein=polymer.protein, ligands=poses, mode=states
    """
    if mode == "auto":
        names = cmd.get_names("objects")
        if ligands in names and cmd.count_states(ligands) > 1:
            mode = "states"
        else:
            mode = "objects"

    if mode == "states":
        _stepper.setup_states(protein, ligands)
        print(f"Contact Inspector: {cmd.count_states(ligands)} states.")
    else:
        if "," in ligands:
            ligs = [l.strip() for l in ligands.split(",")]
        else:
            all_n = cmd.get_names("objects")
            ligs = []
            for n in all_n:
                try:
                    if cmd.count_atoms(f"{n} and ({ligands})") > 0:
                        if cmd.count_atoms(f"{n} and ({protein})") == 0:
                            ligs.append(n)
                except CmdException:
                    pass
            if not ligs:
                ligs = [ligands]
        _stepper.setup_objects(protein, ligs)
        print(f"Contact Inspector: {len(ligs)} ligand(s).")

    _print_summary()
    _bind_keys()


def _print_summary():
    """Print summary to console only when the GUI is not open."""
    if _gui_window is None:
        print(_stepper.summary())

def ci_next():
    _stepper.next(); _print_summary()

def ci_prev():
    _stepper.prev(); _print_summary()

def ci_goto(index=0):
    _stepper.goto(int(index)); _print_summary()

def ci_update():
    _stepper._show_current(); _print_summary()

def ci_clear():
    _clear_all(); _unbind_keys(); print("Contact Inspector: cleared.")

def ci_gui():
    _open_gui()

cmd.extend("ci_setup", ci_setup)
cmd.extend("ci_next", ci_next)
cmd.extend("ci_prev", ci_prev)
cmd.extend("ci_goto", ci_goto)
cmd.extend("ci_update", ci_update)
cmd.extend("ci_clear", ci_clear)
cmd.extend("ci_gui", ci_gui)


# ---------------------------------------------------------------------------
# Qt GUI
# ---------------------------------------------------------------------------

_gui_window = None

def __init_plugin__(app=None):
    from pymol.plugins import addmenuitemqt
    addmenuitemqt("Contact Inspector", _open_gui)


def _open_gui():
    global _gui_window, _stepper

    if _gui_window is not None:
        try:
            _gui_window.raise_(); _gui_window.activateWindow(); return
        except RuntimeError:
            _gui_window = None

    from pymol.Qt import QtWidgets, QtCore, QtGui

    win = QtWidgets.QWidget()
    _gui_window = win
    win.setWindowTitle("Contact Inspector")
    win.setMinimumWidth(440)
    win.setAttribute(QtCore.Qt.WA_DeleteOnClose)
    win.destroyed.connect(lambda: _set_gui_none())

    root = QtWidgets.QVBoxLayout(win)
    root.setContentsMargins(8, 8, 8, 8)
    root.setSpacing(5)

    # Setup
    g_s = QtWidgets.QGroupBox("Setup")
    l_s = QtWidgets.QGridLayout(g_s)
    l_s.addWidget(QtWidgets.QLabel("Protein:"), 0, 0)
    e_prot = QtWidgets.QLineEdit("polymer.protein")
    l_s.addWidget(e_prot, 0, 1)
    l_s.addWidget(QtWidgets.QLabel("Ligand(s):"), 1, 0)
    e_lig = QtWidgets.QLineEdit("organic")
    l_s.addWidget(e_lig, 1, 1)

    hl_m = QtWidgets.QHBoxLayout()
    hl_m.addWidget(QtWidgets.QLabel("Mode:"))
    bg_m = QtWidgets.QButtonGroup(win)
    rb_m = {}
    for m in ("auto", "objects", "states"):
        rb = QtWidgets.QRadioButton(m)
        if m == "auto": rb.setChecked(True)
        bg_m.addButton(rb); rb_m[m] = rb; hl_m.addWidget(rb)
    hl_m.addStretch()
    l_s.addLayout(hl_m, 2, 0, 1, 2)

    hl_b = QtWidgets.QHBoxLayout()
    b_setup = QtWidgets.QPushButton("Setup")
    b_clear = QtWidgets.QPushButton("Clear")
    hl_b.addWidget(b_setup); hl_b.addWidget(b_clear)
    l_s.addLayout(hl_b, 3, 0, 1, 2)
    root.addWidget(g_s)

    # Navigate
    g_n = QtWidgets.QGroupBox("Navigate")
    l_n = QtWidgets.QVBoxLayout(g_n)
    lbl = QtWidgets.QLabel("Ready - click Setup")
    lbl.setAlignment(QtCore.Qt.AlignCenter)
    f = lbl.font(); f.setBold(True); f.setPointSize(f.pointSize()+1)
    lbl.setFont(f)
    l_n.addWidget(lbl)

    hl_pn = QtWidgets.QHBoxLayout()
    b_prev = QtWidgets.QPushButton("<  Prev")
    b_next = QtWidgets.QPushButton("Next  >")
    hl_pn.addWidget(b_prev); hl_pn.addWidget(b_next)
    l_n.addLayout(hl_pn)

    hl_j = QtWidgets.QHBoxLayout()
    hl_j.addWidget(QtWidgets.QLabel("Go to #:"))
    sp = QtWidgets.QSpinBox(); sp.setMinimum(0); sp.setMaximum(99999)
    sp.setFixedWidth(70); hl_j.addWidget(sp)
    b_go = QtWidgets.QPushButton("Go"); b_go.setFixedWidth(50)
    hl_j.addWidget(b_go); hl_j.addStretch()
    l_n.addLayout(hl_j)
    root.addWidget(g_n)

    # Swatch+checkbox helper
    def _cb(lay, text, hex_color, checked=True):
        row = QtWidgets.QHBoxLayout()
        sw = QtWidgets.QLabel(); sw.setFixedSize(14, 14)
        sw.setStyleSheet(f"background-color: {hex_color}; border: none;")
        row.addWidget(sw)
        cb = QtWidgets.QCheckBox(text); cb.setChecked(checked)
        row.addWidget(cb); row.addStretch()
        lay.addLayout(row)
        return cb

    # Non-covalent bonds
    g1 = QtWidgets.QGroupBox("Non-covalent bonds")
    l1 = QtWidgets.QVBoxLayout(g1)
    cb_hb = _cb(l1, "Hydrogen bonds",  "#ffd900")
    cb_xb = _cb(l1, "Halogen bonds",   "#9933e6")
    cb_sb = _cb(l1, "Salt bridges",    "#e633e6")
    cb_ah = _cb(l1, "Aromatic H-Bond", "#4dd97f")
    root.addWidget(g1)

    # Pi interactions
    g2 = QtWidgets.QGroupBox("Pi interactions")
    l2 = QtWidgets.QVBoxLayout(g2)
    cb_pp = _cb(l2, "Pi-pi stacking",  "#4dc0ff")
    cb_pc = _cb(l2, "Pi-cation",       "#33cc33")
    root.addWidget(g2)

    # Contacts / clashes (off by default)
    g3 = QtWidgets.QGroupBox("Contacts / Clashes")
    l3 = QtWidgets.QVBoxLayout(g3)
    cb_cg = _cb(l3, "Good",  "#33cc33", checked=False)
    cb_cb = _cb(l3, "Bad",   "#ff9900", checked=False)
    cb_cu = _cb(l3, "Ugly",  "#ff2626", checked=False)
    root.addWidget(g3)

    # Labels
    cb_lb = QtWidgets.QCheckBox("Show distance labels")
    cb_lb.setChecked(True)
    root.addWidget(cb_lb)

    # Summary
    g_sum = QtWidgets.QGroupBox("Summary")
    l_sum = QtWidgets.QVBoxLayout(g_sum)
    t_sum = QtWidgets.QTextEdit()
    t_sum.setReadOnly(True)
    t_sum.setFont(QtGui.QFont("Courier", 9))
    t_sum.setMaximumHeight(180)
    l_sum.addWidget(t_sum)
    root.addWidget(g_sum, stretch=1)

    # Callbacks
    def refresh():
        c = _stepper._count()
        if c > 0:
            lbl.setText(f"{_stepper._label()}  "
                        f"({_stepper.current_index + 1}/{c})")
            sp.setMaximum(c - 1)
        else:
            lbl.setText("Ready - click Setup")
        t_sum.setPlainText(_stepper.summary())

    def do_setup():
        mode = next(m for m, rb in rb_m.items() if rb.isChecked())
        ci_setup(protein=e_prot.text(), ligands=e_lig.text(), mode=mode)
        refresh()

    def do_clear():
        ci_clear()
        lbl.setText("Ready - click Setup"); t_sum.clear()

    def do_prev():
        ci_prev(); refresh()

    def do_next():
        ci_next(); refresh()

    def do_go():
        ci_goto(sp.value()); refresh()

    def _tog(cb, attr):
        def h(state):
            setattr(_stepper, attr, cb.isChecked())
            ci_update(); refresh()
        return h

    b_setup.clicked.connect(do_setup)
    b_clear.clicked.connect(do_clear)
    b_prev.clicked.connect(do_prev)
    b_next.clicked.connect(do_next)
    b_go.clicked.connect(do_go)

    cb_hb.stateChanged.connect(_tog(cb_hb, "show_hbonds"))
    cb_xb.stateChanged.connect(_tog(cb_xb, "show_halogen"))
    cb_sb.stateChanged.connect(_tog(cb_sb, "show_salt"))
    cb_ah.stateChanged.connect(_tog(cb_ah, "show_arom_hb"))
    cb_pp.stateChanged.connect(_tog(cb_pp, "show_pipi"))
    cb_pc.stateChanged.connect(_tog(cb_pc, "show_pi_cation"))
    cb_cg.stateChanged.connect(_tog(cb_cg, "show_clash_good"))
    cb_cb.stateChanged.connect(_tog(cb_cb, "show_clash_bad"))
    cb_cu.stateChanged.connect(_tog(cb_cu, "show_clash_ugly"))
    cb_lb.stateChanged.connect(_tog(cb_lb, "show_labels"))

    for key, fn in [(QtCore.Qt.Key_Right, do_next),
                    (QtCore.Qt.Key_Left, do_prev)]:
        sc = QtWidgets.QShortcut(QtGui.QKeySequence(key), win)
        sc.activated.connect(fn)

    win.show(); win.raise_()


def _set_gui_none():
    global _gui_window
    _gui_window = None


# ---------------------------------------------------------------------------
# Startup
# ---------------------------------------------------------------------------

print("Contact Inspector loaded.")
print("  ci_gui   - open GUI panel")
print("  ci_setup - setup from command line")
print("  LEFT/RIGHT arrow keys after setup")
