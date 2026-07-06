"""Single-file PyMOL Internal Coordinate Editor plugin.

This plugin is self-contained with GUI window. With this plugin, one can modify bond length, angle, and dihedral torsion.

In addition, this distribution file can be copied into PyMOL and loaded with::

    run /path/to/int_coord_edit.py


IntCoordEdit

Author : Yunwen Tao, Ph.D.
email: ywtao.smu@gmail.com
Date: july 2026
License: MIT License
Version 1.0
"""

from __future__ import annotations

import math
import sys
from collections import deque
from collections.abc import Callable, Hashable, Iterable, Mapping
from dataclasses import dataclass, field
from typing import Any, List, Literal, Sequence, TypeVar

try:
    from pymol import cmd as _pymol_cmd
except Exception:
    _pymol_cmd = None

cmd = _pymol_cmd

try:
    from pymol.Qt import QtCore, QtWidgets
except Exception:
    QtCore = None
    QtWidgets = None

try:
    from pymol.wizard import Wizard as _Wizard
except Exception:
    _Wizard = object

class _SelectionModelShim:
    pass


_selection_model = sys.modules.get(__name__) or _SelectionModelShim()


# ---------------------------------------------------------------------------
# Geometry utilities

Vec3 = tuple[float, float, float]
AtomKey = TypeVar("AtomKey", bound=Hashable)

EPSILON = 1.0e-12


class GeometryError(ValueError):
    """Raised when a requested coordinate operation is geometrically invalid."""


def vec3(point: Iterable[float]) -> Vec3:
    x, y, z = point
    return (float(x), float(y), float(z))


def add(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def sub(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def scale(v: Vec3, factor: float) -> Vec3:
    return (v[0] * factor, v[1] * factor, v[2] * factor)


def dot(a: Vec3, b: Vec3) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a: Vec3, b: Vec3) -> Vec3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def norm(v: Vec3) -> float:
    return math.sqrt(dot(v, v))


def normalize(v: Vec3) -> Vec3:
    length = norm(v)
    if length < EPSILON:
        raise GeometryError("Cannot normalize a near-zero vector.")
    return scale(v, 1.0 / length)


def distance(p1: Vec3, p2: Vec3) -> float:
    return norm(sub(p2, p1))


def angle(p1: Vec3, p2: Vec3, p3: Vec3) -> float:
    """Return angle p1-p2-p3 in degrees."""
    v1 = sub(p1, p2)
    v2 = sub(p3, p2)
    n1 = norm(v1)
    n2 = norm(v2)
    if n1 < EPSILON or n2 < EPSILON:
        raise GeometryError("Cannot measure angle with coincident atoms.")
    cosine = dot(v1, v2) / (n1 * n2)
    cosine = max(-1.0, min(1.0, cosine))
    return math.degrees(math.acos(cosine))


def normalize_angle_delta(delta: float) -> float:
    """Normalize a signed angular delta to (-180, 180]."""
    value = (float(delta) + 180.0) % 360.0 - 180.0
    if value <= -180.0 + 1.0e-12:
        return 180.0
    return value


def normalize_dihedral_value(value: float) -> float:
    return normalize_angle_delta(value)


def dihedral(p1: Vec3, p2: Vec3, p3: Vec3, p4: Vec3) -> float:
    """Return signed dihedral p1-p2-p3-p4 in degrees, normalized to (-180, 180]."""
    b0 = sub(p1, p2)
    b1 = sub(p3, p2)
    b2 = sub(p4, p3)
    b1_unit = normalize(b1)

    v = sub(b0, scale(b1_unit, dot(b0, b1_unit)))
    w = sub(b2, scale(b1_unit, dot(b2, b1_unit)))
    if norm(v) < EPSILON or norm(w) < EPSILON:
        raise GeometryError("Cannot measure dihedral with collinear axis atoms.")

    x = dot(v, w)
    y = dot(cross(b1_unit, v), w)
    return normalize_dihedral_value(math.degrees(math.atan2(y, x)))


def rotation_matrix(axis: Vec3, theta_radians: float) -> tuple[Vec3, Vec3, Vec3]:
    """Return a 3x3 Rodrigues rotation matrix for a right-handed axis rotation."""
    ux, uy, uz = normalize(axis)
    c = math.cos(theta_radians)
    s = math.sin(theta_radians)
    t = 1.0 - c
    return (
        (t * ux * ux + c, t * ux * uy - s * uz, t * ux * uz + s * uy),
        (t * ux * uy + s * uz, t * uy * uy + c, t * uy * uz - s * ux),
        (t * ux * uz - s * uy, t * uy * uz + s * ux, t * uz * uz + c),
    )


def mat_vec_mul(matrix: tuple[Vec3, Vec3, Vec3], vector: Vec3) -> Vec3:
    return (
        dot(matrix[0], vector),
        dot(matrix[1], vector),
        dot(matrix[2], vector),
    )


def rotate_point(point: Vec3, origin: Vec3, axis: Vec3, angle_degrees: float) -> Vec3:
    matrix = rotation_matrix(axis, math.radians(angle_degrees))
    return add(origin, mat_vec_mul(matrix, sub(point, origin)))


def translate_points(
    coords: Mapping[AtomKey, Vec3], atom_keys: Iterable[AtomKey], delta: Vec3
) -> dict[AtomKey, Vec3]:
    updated = dict(coords)
    for key in atom_keys:
        updated[key] = add(updated[key], delta)
    return updated


def rotate_points(
    coords: Mapping[AtomKey, Vec3],
    atom_keys: Iterable[AtomKey],
    origin: Vec3,
    axis: Vec3,
    angle_degrees: float,
) -> dict[AtomKey, Vec3]:
    updated = dict(coords)
    for key in atom_keys:
        updated[key] = rotate_point(updated[key], origin, axis, angle_degrees)
    return updated


def _ratio(ratio: float) -> float:
    value = float(ratio)
    if value < 0.0 or value > 1.0:
        raise GeometryError("Move ratio must be between 0.0 and 1.0.")
    return value


def _weights(first_active: bool, last_active: bool, ratio: float) -> tuple[float, float]:
    ratio = _ratio(ratio)
    if first_active and last_active:
        return ratio, 1.0 - ratio
    if first_active:
        return 1.0, 0.0
    if last_active:
        return 0.0, 1.0
    return 0.0, 0.0


def _key_set(keys: Iterable[AtomKey] | None) -> set[AtomKey]:
    return set(keys or ())


def _reject_overlapping_sets(first: set[AtomKey], last: set[AtomKey]) -> None:
    overlap = first & last
    if overlap:
        raise GeometryError(f"Movable atom groups overlap: {sorted(overlap)!r}")


def set_bond_length(
    coords: Mapping[AtomKey, Vec3],
    atom1: AtomKey,
    atom2: AtomKey,
    target_distance: float,
    movable_first: Iterable[AtomKey] | None = None,
    movable_last: Iterable[AtomKey] | None = None,
    ratio: float = 0.5,
    minimum_distance: float = 0.1,
) -> dict[AtomKey, Vec3]:
    target = float(target_distance)
    if target < minimum_distance:
        raise GeometryError(f"Bond distance must be at least {minimum_distance:g} A.")

    first = _key_set(movable_first)
    last = _key_set(movable_last)
    _reject_overlapping_sets(first, last)
    first_weight, last_weight = _weights(bool(first), bool(last), ratio)
    if first_weight == 0.0 and last_weight == 0.0:
        return dict(coords)

    p1 = coords[atom1]
    p2 = coords[atom2]
    current = distance(p1, p2)
    if current < EPSILON:
        raise GeometryError("Cannot set bond length for coincident atoms.")

    unit = normalize(sub(p2, p1))
    delta = target - current
    updated = dict(coords)
    first_delta = scale(unit, -delta * first_weight)
    last_delta = scale(unit, delta * last_weight)
    for key in first:
        updated[key] = add(updated[key], first_delta)
    for key in last:
        updated[key] = add(updated[key], last_delta)
    return updated


def set_angle(
    coords: Mapping[AtomKey, Vec3],
    atom1: AtomKey,
    atom2: AtomKey,
    atom3: AtomKey,
    target_angle: float,
    movable_first: Iterable[AtomKey] | None = None,
    movable_last: Iterable[AtomKey] | None = None,
    ratio: float = 0.5,
) -> dict[AtomKey, Vec3]:
    target = float(target_angle)
    if target <= 0.0 or target >= 180.0:
        raise GeometryError("Bond angle target must be strictly between 0 and 180 degrees.")

    first = _key_set(movable_first)
    last = _key_set(movable_last)
    _reject_overlapping_sets(first, last)
    first_weight, last_weight = _weights(bool(first), bool(last), ratio)
    if first_weight == 0.0 and last_weight == 0.0:
        return dict(coords)

    p1 = coords[atom1]
    p2 = coords[atom2]
    p3 = coords[atom3]
    axis = cross(sub(p1, p2), sub(p3, p2))
    if norm(axis) < EPSILON:
        raise GeometryError("Cannot set angle when the three atoms are near-collinear.")

    delta = target - angle(p1, p2, p3)
    updated = dict(coords)
    if first_weight:
        for key in first:
            updated[key] = rotate_point(updated[key], p2, axis, -delta * first_weight)
    if last_weight:
        for key in last:
            updated[key] = rotate_point(updated[key], p2, axis, delta * last_weight)
    return updated


def set_dihedral(
    coords: Mapping[AtomKey, Vec3],
    atom1: AtomKey,
    atom2: AtomKey,
    atom3: AtomKey,
    atom4: AtomKey,
    target_dihedral: float,
    movable_first: Iterable[AtomKey] | None = None,
    movable_last: Iterable[AtomKey] | None = None,
    ratio: float = 0.5,
) -> dict[AtomKey, Vec3]:
    first = _key_set(movable_first)
    last = _key_set(movable_last)
    _reject_overlapping_sets(first, last)
    first_weight, last_weight = _weights(bool(first), bool(last), ratio)
    if first_weight == 0.0 and last_weight == 0.0:
        return dict(coords)

    p1 = coords[atom1]
    p2 = coords[atom2]
    p3 = coords[atom3]
    p4 = coords[atom4]
    axis = sub(p3, p2)
    if norm(axis) < EPSILON:
        raise GeometryError("Cannot set dihedral with coincident axis atoms.")

    delta = normalize_angle_delta(float(target_dihedral) - dihedral(p1, p2, p3, p4))
    updated = dict(coords)
    if first_weight:
        for key in first:
            updated[key] = rotate_point(updated[key], p2, axis, -delta * first_weight)
    if last_weight:
        for key in last:
            updated[key] = rotate_point(updated[key], p2, axis, delta * last_weight)
    return updated

# ---------------------------------------------------------------------------
# Connectivity helpers

AtomKey = TypeVar("AtomKey", bound=Hashable)
Edge = frozenset[AtomKey]


@dataclass(frozen=True)
class FragmentGroups:
    first: set[AtomKey]
    last: set[AtomKey]
    ring_risk: bool = False
    message: str = ""


def edge_key(a: AtomKey, b: AtomKey) -> Edge:
    return frozenset((a, b))


def build_adjacency(
    bonds: Iterable[tuple[AtomKey, AtomKey]],
    atoms: Iterable[AtomKey] | None = None,
) -> dict[AtomKey, set[AtomKey]]:
    adjacency: dict[AtomKey, set[AtomKey]] = {}
    if atoms is not None:
        for atom in atoms:
            adjacency.setdefault(atom, set())
    for a, b in bonds:
        adjacency.setdefault(a, set()).add(b)
        adjacency.setdefault(b, set()).add(a)
    return adjacency


def bfs_component(
    adjacency: Mapping[AtomKey, set[AtomKey]],
    start: AtomKey,
    blocked_nodes: Iterable[AtomKey] | None = None,
    blocked_edges: Iterable[Edge] | None = None,
) -> set[AtomKey]:
    blocked_node_set = set(blocked_nodes or ())
    blocked_edge_set = set(blocked_edges or ())
    if start in blocked_node_set:
        return set()

    visited: set[AtomKey] = {start}
    queue: deque[AtomKey] = deque([start])
    while queue:
        node = queue.popleft()
        for neighbor in adjacency.get(node, set()):
            if neighbor in blocked_node_set:
                continue
            if edge_key(node, neighbor) in blocked_edge_set:
                continue
            if neighbor in visited:
                continue
            visited.add(neighbor)
            queue.append(neighbor)
    return visited


def connected_without_edge(
    adjacency: Mapping[AtomKey, set[AtomKey]], atom1: AtomKey, atom2: AtomKey
) -> bool:
    return atom2 in bfs_component(adjacency, atom1, blocked_edges={edge_key(atom1, atom2)})


def bond_side_groups(
    adjacency: Mapping[AtomKey, set[AtomKey]], atom1: AtomKey, atom2: AtomKey
) -> FragmentGroups:
    blocked = {edge_key(atom1, atom2)}
    first = bfs_component(adjacency, atom1, blocked_edges=blocked)
    last = bfs_component(adjacency, atom2, blocked_edges=blocked)
    overlap = first & last
    ring_risk = bool(overlap)
    message = (
        "Deleting the selected bond does not separate the two atoms; group movement is ambiguous."
        if ring_risk
        else ""
    )
    return FragmentGroups(first=first, last=last, ring_risk=ring_risk, message=message)


def angle_side_groups(
    adjacency: Mapping[AtomKey, set[AtomKey]], atom1: AtomKey, atom2: AtomKey, atom3: AtomKey
) -> FragmentGroups:
    first = bfs_component(adjacency, atom1, blocked_nodes={atom2})
    last = bfs_component(adjacency, atom3, blocked_nodes={atom2})
    overlap = first & last
    ring_risk = bool(overlap)
    message = (
        "The two angle arms remain connected without the center atom; group movement is ambiguous."
        if ring_risk
        else ""
    )
    return FragmentGroups(first=first, last=last, ring_risk=ring_risk, message=message)


def dihedral_side_groups(
    adjacency: Mapping[AtomKey, set[AtomKey]],
    atom1: AtomKey,
    atom2: AtomKey,
    atom3: AtomKey,
    atom4: AtomKey,
) -> FragmentGroups:
    blocked = {edge_key(atom2, atom3)}
    left_full = bfs_component(adjacency, atom2, blocked_edges=blocked)
    right_full = bfs_component(adjacency, atom3, blocked_edges=blocked)
    overlap = left_full & right_full
    ring_risk = bool(overlap)
    first = set(left_full) - {atom2, atom3}
    last = set(right_full) - {atom2, atom3}
    message = (
        "Deleting the central dihedral bond does not separate the two sides; group movement is ambiguous."
        if ring_risk
        else ""
    )
    return FragmentGroups(first=first, last=last, ring_risk=ring_risk, message=message)


def require_unambiguous_groups(groups: FragmentGroups, allow_ring: bool = False) -> None:
    if groups.ring_risk and not allow_ring:
        raise ValueError(groups.message or "Group movement is ambiguous for this connectivity.")

# ---------------------------------------------------------------------------
# Selection model

ICMode = Literal["bond", "angle", "dihedral"]
MoveSide = Literal["first", "last", "both", "none"]
Granularity = Literal["atom", "group"]


class SelectionError(ValueError):
    """Raised when selected atoms do not form a valid edit target."""


@dataclass(frozen=True)
class AtomRef:
    """Stable enough PyMOL atom identity for a single object/state edit."""

    object: str
    index: int
    state: int = 1
    id: int | None = None
    chain: str = ""
    resi: str = ""
    name: str = ""
    segi: str = ""

    @property
    def selection(self) -> str:
        return f"({object_selector(self.object)} and index {self.index})"

    @property
    def key(self) -> tuple[str, int, int]:
        return (self.object, int(self.index), int(self.state))

    def label(self) -> str:
        parts = [self.object]
        residue = f"{self.chain}/{self.resi}".strip("/")
        if residue:
            parts.append(residue)
        if self.name:
            parts.append(self.name)
        parts.append(f"index {self.index}")
        return " ".join(parts)


@dataclass
class ICSelection:
    mode: ICMode
    atoms: List[AtomRef] = field(default_factory=list)
    side: MoveSide = "last"
    granularity: Granularity = "group"

    @property
    def expected_count(self) -> int:
        return expected_pick_count(self.mode)

    @property
    def is_complete(self) -> bool:
        return len(self.atoms) == self.expected_count

    def append(self, atom: AtomRef) -> None:
        if any(existing.key == atom.key for existing in self.atoms):
            raise SelectionError("The same atom cannot be picked twice.")
        if self.atoms:
            first = self.atoms[0]
            if atom.object != first.object:
                raise SelectionError("Cross-object picks are not supported in the MVP.")
            if atom.state != first.state:
                raise SelectionError("All picked atoms must use the same state.")
        if len(self.atoms) >= self.expected_count:
            raise SelectionError("This selection is already complete.")
        self.atoms.append(atom)

    def clear(self) -> None:
        self.atoms.clear()

    def pop(self) -> AtomRef | None:
        if not self.atoms:
            return None
        return self.atoms.pop()


def normalize_mode(mode: str) -> ICMode:
    value = str(mode).strip().lower()
    aliases = {"torsion": "dihedral", "distance": "bond"}
    value = aliases.get(value, value)
    if value not in {"bond", "angle", "dihedral"}:
        raise SelectionError(f"Unsupported internal-coordinate mode: {mode!r}")
    return value  # type: ignore[return-value]


def normalize_side(side: str) -> MoveSide:
    value = str(side).strip().lower().replace("-", "_").replace(" ", "_")
    aliases = {
        "atom1": "first",
        "a": "first",
        "left": "first",
        "first": "first",
        "atom_1": "first",
        "atom2": "last",
        "atom3": "last",
        "atom4": "last",
        "last": "last",
        "right": "last",
        "d": "last",
        "both": "both",
        "all": "both",
        "measure": "none",
        "fixed": "none",
        "none": "none",
        "no": "none",
        "0": "none",
    }
    if value not in aliases:
        raise SelectionError(f"Unsupported move side: {side!r}")
    return aliases[value]  # type: ignore[return-value]


def normalize_granularity(mode: str) -> Granularity:
    value = str(mode).strip().lower()
    aliases = {
        "atom": "atom",
        "single": "atom",
        "endpoint": "atom",
        "group": "group",
        "fragment": "group",
        "frag": "group",
    }
    if value not in aliases:
        raise SelectionError(f"Unsupported move granularity: {mode!r}")
    return aliases[value]  # type: ignore[return-value]


def expected_pick_count(mode: str) -> int:
    normalized = normalize_mode(mode)
    return {"bond": 2, "angle": 3, "dihedral": 4}[normalized]


def prompt_for_pick(mode: str, picked_count: int) -> str:
    normalized = normalize_mode(mode)
    prompts = {
        "bond": ["Bond mode: pick atom 1", "Bond mode: pick atom 2"],
        "angle": [
            "Angle mode: pick atom 1",
            "Angle mode: pick center atom 2",
            "Angle mode: pick atom 3",
        ],
        "dihedral": [
            "Dihedral mode: pick atom 1",
            "Dihedral mode: pick atom 2",
            "Dihedral mode: pick atom 3",
            "Dihedral mode: pick atom 4",
        ],
    }
    mode_prompts = prompts[normalized]
    if picked_count >= len(mode_prompts):
        return f"{normalized.capitalize()} mode: ready to edit"
    return mode_prompts[picked_count]


def side_labels(mode: str) -> list[tuple[str, MoveSide]]:
    normalized = normalize_mode(mode)
    if normalized == "bond":
        return [
            ("Atom 1 side only", "first"),
            ("Atom 2 side only", "last"),
            ("Move both sides", "both"),
            ("Fixed / Measure only", "none"),
        ]
    if normalized == "angle":
        return [
            ("Atom 1 side only", "first"),
            ("Atom 3 side only", "last"),
            ("Move both sides", "both"),
            ("Fixed / Measure only", "none"),
        ]
    return [
        ("Atom 1/2 side only", "first"),
        ("Atom 3/4 side only", "last"),
        ("Move both sides", "both"),
        ("Fixed / Measure only", "none"),
    ]


def validate_distinct_atom_indices(indices: Sequence[int]) -> None:
    if len(set(indices)) != len(indices):
        raise SelectionError("Picked atoms must be distinct.")


def object_selector(object_name: str) -> str:
    """Return a PyMOL selection fragment for an object/model name."""
    return f"model {object_name}"

# ---------------------------------------------------------------------------
# Coordinate snapshots

def _default_object_selector(object_name: str) -> str:
    return f"model {object_name}"


object_selector = getattr(_selection_model, "object_selector", _default_object_selector)
if not hasattr(_selection_model, "object_selector"):
    _selection_model.object_selector = object_selector


@dataclass
class CoordinateSnapshot:
    object_name: str
    state: int
    coords: dict[int, Vec3]

    @classmethod
    def capture(cls, cmd, object_name: str, state: int) -> "CoordinateSnapshot":
        model = cmd.get_model(f"({object_selector(object_name)})", state=state)
        coords = {int(atom.index): vec3(atom.coord) for atom in model.atom}
        return cls(object_name=object_name, state=int(state), coords=coords)

    def restore(self, cmd) -> None:
        for index, (x, y, z) in self.coords.items():
            cmd.alter_state(
                self.state,
                f"({object_selector(self.object_name)} and index {index})",
                f"(x, y, z)=({x:.17g}, {y:.17g}, {z:.17g})",
            )
        cmd.rebuild()


@dataclass
class MemoryCoordinateState:
    coords: dict[Hashable, Vec3]

    @classmethod
    def capture(cls, coords: Mapping[Hashable, Vec3]) -> "MemoryCoordinateState":
        return cls(coords=dict(coords))

    def restored(self) -> dict[Hashable, Vec3]:
        return dict(self.coords)

# ---------------------------------------------------------------------------
# PyMOL command API

class CommandError(RuntimeError):
    """Raised when a PyMOL command cannot be completed."""


def _default_object_selector(object_name: str) -> str:
    return f"model {object_name}"


object_selector = getattr(_selection_model, "object_selector", _default_object_selector)
if not hasattr(_selection_model, "object_selector"):
    _selection_model.object_selector = object_selector


def _get_cmd(cmd_module=None):
    command = cmd_module or _pymol_cmd
    if command is None:
        raise CommandError("PyMOL is not available. Run these commands inside PyMOL.")
    return command


def _bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def _state(cmd, state: str | int) -> int:
    if state is None or str(state).strip().lower() in {"", "current"}:
        return int(cmd.get_state())
    if str(state).strip().lower() == "all":
        raise CommandError("state=all is not implemented in the MVP; use state=current or an integer.")
    return int(state)


def _resolve_atom(cmd, selection: str, state: int) -> AtomRef:
    model = cmd.get_model(f"({selection})", state=state)
    if len(model.atom) != 1:
        raise SelectionError(f"Selection must resolve to exactly one atom: {selection!r}")
    atom = model.atom[0]
    object_name = getattr(atom, "model", "") or _object_for_selection(cmd, selection)
    if not object_name:
        raise SelectionError(f"Could not determine object for selection: {selection!r}")
    return AtomRef(
        object=object_name,
        index=int(atom.index),
        id=getattr(atom, "id", None),
        chain=getattr(atom, "chain", "") or "",
        resi=getattr(atom, "resi", "") or "",
        name=getattr(atom, "name", "") or "",
        segi=getattr(atom, "segi", "") or "",
        state=state,
    )


def _object_for_selection(cmd, selection: str) -> str:
    names = cmd.get_object_list(f"({selection})")
    return names[0] if names else ""


def _validate_refs(refs: list[AtomRef]) -> None:
    if not refs:
        raise SelectionError("No atoms were selected.")
    validate_distinct_atom_indices([ref.index for ref in refs])
    object_name = refs[0].object
    state = refs[0].state
    for ref in refs[1:]:
        if ref.object != object_name:
            raise SelectionError("Cross-object edits are not supported in the MVP.")
        if ref.state != state:
            raise SelectionError("All selected atoms must use the same state.")


def _unique_object_name(cmd, base: str) -> str:
    existing = set(cmd.get_names("objects"))
    candidate = base
    counter = 1
    while candidate in existing:
        counter += 1
        candidate = f"{base}_{counter}"
    return candidate


def _maybe_copy_object(cmd, object_name: str, state: int, object_mode: str, preview: bool) -> tuple[str, int]:
    mode = str(object_mode).strip().lower()
    if preview:
        new_name = _unique_object_name(cmd, f"{object_name}_ic_preview")
        cmd.create(new_name, f"({object_selector(object_name)})", state, 1)
        return new_name, 1
    if mode == "inplace":
        return object_name, state
    if mode == "copy":
        new_name = _unique_object_name(cmd, f"{object_name}_ic")
        cmd.create(new_name, f"({object_selector(object_name)})", state, 1)
        return new_name, 1
    raise CommandError("object_mode must be 'inplace' or 'copy'.")


def _model_coords_and_bonds(cmd, object_name: str, state: int) -> tuple[dict[int, Vec3], dict[int, set[int]]]:
    model = cmd.get_model(f"({object_selector(object_name)})", state=state)
    coords: dict[int, Vec3] = {}
    position_to_index: list[int] = []
    for atom in model.atom:
        index = int(atom.index)
        position_to_index.append(index)
        coords[index] = vec3(atom.coord)

    bonds: list[tuple[int, int]] = []
    for bond in getattr(model, "bond", []):
        raw = list(getattr(bond, "index", []))
        if len(raw) != 2:
            continue
        a, b = int(raw[0]), int(raw[1])
        if 0 <= a < len(position_to_index) and 0 <= b < len(position_to_index):
            bonds.append((position_to_index[a], position_to_index[b]))
        elif a in coords and b in coords:
            bonds.append((a, b))
    return coords, build_adjacency(bonds, atoms=coords.keys())


def _write_coords(cmd, object_name: str, state: int, coords: Mapping[int, Vec3], changed: Iterable[int]) -> None:
    for index in sorted(set(changed)):
        x, y, z = coords[index]
        cmd.alter_state(
            state,
            f"({object_selector(object_name)} and index {index})",
            f"(x, y, z)=({x:.17g}, {y:.17g}, {z:.17g})",
        )
    cmd.rebuild()


def _select_movable(
    side: str,
    granularity: str,
    adjacency: Mapping[int, set[int]],
    indices: list[int],
    allow_ring: bool,
) -> tuple[set[int], set[int]]:
    normalized_side = normalize_side(side)
    if normalized_side == "none":
        return set(), set()

    normalized_granularity = normalize_granularity(granularity)
    if normalized_granularity == "atom":
        first = {indices[0]}
        last = {indices[-1]}
    elif len(indices) == 2:
        groups = bond_side_groups(adjacency, indices[0], indices[1])
        require_unambiguous_groups(groups, allow_ring=allow_ring)
        first, last = groups.first, groups.last
    elif len(indices) == 3:
        groups = angle_side_groups(adjacency, indices[0], indices[1], indices[2])
        require_unambiguous_groups(groups, allow_ring=allow_ring)
        first, last = groups.first, groups.last
    elif len(indices) == 4:
        groups = dihedral_side_groups(adjacency, indices[0], indices[1], indices[2], indices[3])
        require_unambiguous_groups(groups, allow_ring=allow_ring)
        first, last = groups.first, groups.last
    else:
        raise CommandError("Internal-coordinate edits require 2, 3, or 4 atoms.")

    if normalized_side == "first":
        return first, set()
    if normalized_side == "last":
        return set(), last
    if normalized_side == "both":
        return first, last
    raise CommandError(f"Unsupported move side: {side!r}")


def _remap_refs_to_copy(refs: list[AtomRef], object_name: str, state: int) -> list[AtomRef]:
    return [
        AtomRef(
            object=object_name,
            index=ref.index,
            state=state,
            id=ref.id,
            chain=ref.chain,
            resi=ref.resi,
            name=ref.name,
            segi=ref.segi,
        )
        for ref in refs
    ]


def _report(cmd, message: str) -> None:
    print(message)


def ic_bond(
    atom1: str,
    atom2: str,
    target_distance: float,
    side: str = "last",
    mode: str = "group",
    state: str | int = "current",
    object_mode: str = "inplace",
    preview: int | bool = 0,
    ratio: float = 0.5,
    force: int | bool = 0,
    cmd_module=None,
) -> float:
    cmd = _get_cmd(cmd_module)
    edit_state = _state(cmd, state)
    refs = [_resolve_atom(cmd, atom1, edit_state), _resolve_atom(cmd, atom2, edit_state)]
    _validate_refs(refs)
    target_object, target_state = _maybe_copy_object(cmd, refs[0].object, edit_state, object_mode, _bool(preview))
    refs = _remap_refs_to_copy(refs, target_object, target_state)
    coords, adjacency = _model_coords_and_bonds(cmd, target_object, target_state)
    indices = [ref.index for ref in refs]
    movable_first, movable_last = _select_movable(side, mode, adjacency, indices, allow_ring=_bool(force))
    updated = set_bond_length(
        coords, indices[0], indices[1], float(target_distance), movable_first, movable_last, float(ratio)
    )
    _write_coords(cmd, target_object, target_state, updated, movable_first | movable_last)
    value = distance(updated[indices[0]], updated[indices[1]])
    _report(cmd, f"ic_bond: {target_object} distance = {value:.6f} A")
    return value


def ic_angle(
    atom1: str,
    atom2: str,
    atom3: str,
    target_angle: float,
    side: str = "last",
    mode: str = "group",
    state: str | int = "current",
    object_mode: str = "inplace",
    preview: int | bool = 0,
    ratio: float = 0.5,
    force: int | bool = 0,
    cmd_module=None,
) -> float:
    cmd = _get_cmd(cmd_module)
    edit_state = _state(cmd, state)
    refs = [
        _resolve_atom(cmd, atom1, edit_state),
        _resolve_atom(cmd, atom2, edit_state),
        _resolve_atom(cmd, atom3, edit_state),
    ]
    _validate_refs(refs)
    target_object, target_state = _maybe_copy_object(cmd, refs[0].object, edit_state, object_mode, _bool(preview))
    refs = _remap_refs_to_copy(refs, target_object, target_state)
    coords, adjacency = _model_coords_and_bonds(cmd, target_object, target_state)
    indices = [ref.index for ref in refs]
    movable_first, movable_last = _select_movable(side, mode, adjacency, indices, allow_ring=_bool(force))
    updated = set_angle(coords, indices[0], indices[1], indices[2], float(target_angle), movable_first, movable_last, float(ratio))
    _write_coords(cmd, target_object, target_state, updated, movable_first | movable_last)
    value = measure_angle(updated[indices[0]], updated[indices[1]], updated[indices[2]])
    _report(cmd, f"ic_angle: {target_object} angle = {value:.6f} deg")
    return value


def ic_dihedral(
    atom1: str,
    atom2: str,
    atom3: str,
    atom4: str,
    target_dihedral: float,
    side: str = "last",
    mode: str = "group",
    state: str | int = "current",
    object_mode: str = "inplace",
    preview: int | bool = 0,
    ratio: float = 0.5,
    force: int | bool = 0,
    cmd_module=None,
) -> float:
    cmd = _get_cmd(cmd_module)
    edit_state = _state(cmd, state)
    refs = [
        _resolve_atom(cmd, atom1, edit_state),
        _resolve_atom(cmd, atom2, edit_state),
        _resolve_atom(cmd, atom3, edit_state),
        _resolve_atom(cmd, atom4, edit_state),
    ]
    _validate_refs(refs)
    target_object, target_state = _maybe_copy_object(cmd, refs[0].object, edit_state, object_mode, _bool(preview))
    refs = _remap_refs_to_copy(refs, target_object, target_state)
    coords, adjacency = _model_coords_and_bonds(cmd, target_object, target_state)
    indices = [ref.index for ref in refs]
    movable_first, movable_last = _select_movable(side, mode, adjacency, indices, allow_ring=_bool(force))
    updated = set_dihedral(
        coords,
        indices[0],
        indices[1],
        indices[2],
        indices[3],
        float(target_dihedral),
        movable_first,
        movable_last,
        float(ratio),
    )
    _write_coords(cmd, target_object, target_state, updated, movable_first | movable_last)
    value = measure_dihedral(updated[indices[0]], updated[indices[1]], updated[indices[2]], updated[indices[3]])
    _report(cmd, f"ic_dihedral: {target_object} dihedral = {value:.6f} deg")
    return value


def ic_pick(mode: str = "dihedral", cmd_module=None) -> None:
    cmd = _get_cmd(cmd_module)
    normalized = normalize_mode(mode)
    wizard = ICPickWizard(mode=normalized, cmd_module=cmd)
    set_active_wizard(cmd, wizard)


def register_commands(cmd_module=None) -> None:
    cmd = _get_cmd(cmd_module)
    cmd.extend("ic_bond", ic_bond)
    cmd.extend("ic_angle", ic_angle)
    cmd.extend("ic_dihedral", ic_dihedral)
    cmd.extend("ic_pick", ic_pick)

# ---------------------------------------------------------------------------
# Pick wizard

PICK_HIGHLIGHT_NAMES = tuple(f"_ic_pick_{number}" for number in range(1, 5))
TRANSIENT_PICK_SELECTIONS = ("sele", "pk1", "pk2", "pkbond")


def clear_pick_artifacts(cmd_module) -> None:
    """Remove temporary pick labels/selections left by the picker."""
    if cmd_module is None:
        return
    for name in PICK_HIGHLIGHT_NAMES:
        _clear_label(cmd_module, name)
        _delete_selection(cmd_module, name)


class ICPickWizard(_Wizard):
    """Wizard state machine for bond, angle, and dihedral atom picks."""

    def __init__(
        self,
        mode: str = "dihedral",
        on_complete: Callable[[list[AtomRef]], None] | None = None,
        cmd_module=None,
        state: int | None = None,
        auto_exit: bool = False,
    ) -> None:
        if _Wizard is not object:
            super().__init__(_pymol_cmd or cmd_module)
        self.cmd = cmd_module or _pymol_cmd
        self.selection = ICSelection(mode=normalize_mode(mode))
        self.on_complete = on_complete
        self.state = state
        self.auto_exit = auto_exit
        self._highlight_names: list[str] = []
        self._previous_mouse_selection_mode = None
        self._set_atom_selection_mode()
        self._message(prompt_for_pick(self.selection.mode, 0))

    def get_prompt(self):
        return [prompt_for_pick(self.selection.mode, len(self.selection.atoms))]

    def get_panel(self):
        return [
            [1, "Internal Coordinate Picker", ""],
            [2, "Undo Pick", "cmd.get_wizard().undo_pick()"],
            [2, "Clear Picks", "cmd.get_wizard().clear_picks()"],
            [2, "Cancel", "cmd.get_wizard().cancel()"],
        ]

    def do_select(self, name):
        return self._handle_selection(name)

    def do_pick(self, bondFlag=0):  # noqa: N803 - PyMOL Wizard API name
        return self._handle_selection("pk1")

    def clear_picks(self) -> None:
        self.selection.clear()
        self._clear_highlights()
        self._message(prompt_for_pick(self.selection.mode, 0))
        self._refresh()

    def undo_pick(self) -> None:
        removed = self.selection.pop()
        if removed is not None and self._highlight_names:
            name = self._highlight_names.pop()
            _clear_label(self.cmd, name)
            _delete_selection(self.cmd, name)
        self._message(prompt_for_pick(self.selection.mode, len(self.selection.atoms)))
        self._refresh()

    def cancel(self) -> None:
        self.cleanup()
        if self.cmd is not None:
            try:
                self.cmd.set_wizard()
            except Exception:
                self.cmd.wizard(None)

    def cleanup(self) -> None:
        self._clear_highlights()
        self._restore_mouse_selection_mode()

    def _handle_selection(self, name: str) -> None:
        if self.cmd is None:
            raise RuntimeError("PyMOL command module is not available.")
        try:
            atom = self._atom_from_selection(name)
            self.selection.append(atom)
            self._highlight(atom, len(self.selection.atoms))
        except Exception as exc:
            self._message(str(exc))
            self._refresh()
            return
        finally:
            self._clear_click_state(name)

        if self.selection.is_complete:
            picks = list(self.selection.atoms)
            self._message(f"{self.selection.mode.capitalize()} mode: ready to edit")
            if self.on_complete is not None:
                self.on_complete(picks)
            if self.auto_exit:
                self.cancel()
            else:
                self._restore_mouse_selection_mode()
        else:
            self._message(prompt_for_pick(self.selection.mode, len(self.selection.atoms)))
        self._refresh()

    def _set_atom_selection_mode(self) -> None:
        if self.cmd is None:
            return
        try:
            self._previous_mouse_selection_mode = self.cmd.get("mouse_selection_mode")
            self.cmd.set("mouse_selection_mode", 0)
        except Exception:
            self._previous_mouse_selection_mode = None

    def _restore_mouse_selection_mode(self) -> None:
        if self.cmd is None or self._previous_mouse_selection_mode is None:
            return
        try:
            self.cmd.set("mouse_selection_mode", self._previous_mouse_selection_mode)
        except Exception:
            pass
        self._previous_mouse_selection_mode = None

    def _atom_from_selection(self, selection: str) -> AtomRef:
        state = int(self.state or self.cmd.get_state())
        model = self.cmd.get_model(f"({selection})", state=state)
        if len(model.atom) != 1:
            raise SelectionError("Pick exactly one atom.")
        atom = model.atom[0]
        object_name = getattr(atom, "model", "")
        if not object_name:
            names = self.cmd.get_object_list(f"({selection})")
            object_name = names[0] if names else ""
        if not object_name:
            raise SelectionError("Could not resolve the picked atom's object.")
        return AtomRef(
            object=object_name,
            index=int(atom.index),
            id=getattr(atom, "id", None),
            chain=getattr(atom, "chain", "") or "",
            resi=getattr(atom, "resi", "") or "",
            name=getattr(atom, "name", "") or "",
            segi=getattr(atom, "segi", "") or "",
            state=state,
        )

    def _highlight(self, atom: AtomRef, number: int) -> None:
        name = f"_ic_pick_{number}"
        self._highlight_names.append(name)
        try:
            _clear_label(self.cmd, name)
            _delete_selection(self.cmd, name)
            self.cmd.select(name, atom.selection)
            self.cmd.show("spheres", name)
            self.cmd.label(name, f'"{number}"')
        except Exception:
            pass

    def _clear_highlights(self) -> None:
        if self.cmd is None:
            return
        for name in self._highlight_names:
            _clear_label(self.cmd, name)
            _delete_selection(self.cmd, name)
        self._highlight_names.clear()

    def _clear_click_state(self, source_selection: str) -> None:
        if self.cmd is None:
            return
        try:
            self.cmd.unpick()
        except Exception:
            pass

        try:
            self.cmd.deselect()
        except Exception:
            pass

        transient_names = set(TRANSIENT_PICK_SELECTIONS)
        if source_selection in transient_names:
            transient_names.add(source_selection)
        for name in transient_names:
            _delete_selection(self.cmd, name)

    def _message(self, message: str) -> None:
        print(message)

    def _refresh(self) -> None:
        if self.cmd is not None:
            try:
                self.cmd.refresh_wizard()
            except Exception:
                pass


def set_active_wizard(cmd_module, wizard: ICPickWizard) -> None:
    """Install an already-created Wizard instance in PyMOL.

    PyMOL's ``cmd.wizard("name")`` path imports a wizard by name. For plugin
    instances with callbacks, PyMOL 3.x needs ``cmd.set_wizard(instance)``.
    """
    if hasattr(cmd_module, "set_wizard"):
        cmd_module.set_wizard(wizard)
    else:
        cmd_module.wizard(wizard)
    cmd_module.refresh_wizard()


def _clear_label(cmd_module, selection: str) -> None:
    try:
        cmd_module.label(selection, '""')
    except Exception:
        pass


def _delete_selection(cmd_module, selection: str) -> None:
    try:
        cmd_module.delete(selection)
    except Exception:
        pass

# ---------------------------------------------------------------------------
# Qt panel

_dialog = None
_BaseDialog = QtWidgets.QDialog if QtWidgets is not None else object


def open_dialog() -> None:
    if QtWidgets is None or cmd is None:
        raise RuntimeError("The Internal Coordinate Editor UI must be opened inside PyMOL with PyQt.")

    global _dialog
    if _dialog is None:
        _dialog = InternalCoordinateEditorPanel()
    _dialog.show()
    _dialog.raise_()
    _dialog.activateWindow()


class InternalCoordinateEditorPanel(_BaseDialog):  # type: ignore[misc]
    def __init__(self, parent=None) -> None:
        if QtWidgets is None or QtCore is None or cmd is None:
            raise RuntimeError("The Internal Coordinate Editor UI must be opened inside PyMOL with PyQt.")
        super().__init__(parent)
        self.setWindowTitle("Internal Coordinate Editor")
        self.resize(420, 360)
        self.picks: list[AtomRef] = []
        self.snapshot: CoordinateSnapshot | None = None
        self._wizard = None
        self._building = False
        self._build_ui()
        self._update_mode_controls()

    def _build_ui(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)

        form = QtWidgets.QFormLayout()
        self.mode_combo = QtWidgets.QComboBox()
        self.mode_combo.addItems(["Bond", "Angle", "Dihedral"])
        self.mode_combo.currentTextChanged.connect(self._update_mode_controls)
        form.addRow("Mode", self.mode_combo)

        self.side_combo = QtWidgets.QComboBox()
        form.addRow("Move side", self.side_combo)

        self.granularity_combo = QtWidgets.QComboBox()
        self.granularity_combo.addItem("Atom", "atom")
        self.granularity_combo.addItem("Group / Fragment", "group")
        self.granularity_combo.setCurrentIndex(1)
        form.addRow("Move granularity", self.granularity_combo)

        layout.addLayout(form)

        self.start_button = QtWidgets.QPushButton("Start Picking")
        self.start_button.clicked.connect(self._start_picking)
        self.start_button.setAutoDefault(False)
        self.start_button.setDefault(False)
        layout.addWidget(self.start_button)

        self.picked_list = QtWidgets.QListWidget()
        self.picked_list.setMinimumHeight(80)
        layout.addWidget(self.picked_list)

        value_form = QtWidgets.QFormLayout()
        self.current_label = QtWidgets.QLabel("-")
        value_form.addRow("Current value", self.current_label)
        self.target_edit = QtWidgets.QLineEdit()
        self.target_edit.editingFinished.connect(self._target_edit_finished)
        value_form.addRow("Target value", self.target_edit)
        layout.addLayout(value_form)

        self.slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.slider.valueChanged.connect(self._slider_changed)
        layout.addWidget(self.slider)

        self.auto_preview = QtWidgets.QCheckBox("Auto Preview")
        self.auto_preview.setChecked(True)
        layout.addWidget(self.auto_preview)

        button_row = QtWidgets.QHBoxLayout()
        self.apply_button = QtWidgets.QPushButton("Apply")
        self.apply_button.clicked.connect(self._apply)
        self.reset_button = QtWidgets.QPushButton("Reset")
        self.reset_button.clicked.connect(self._reset)
        self.cancel_button = QtWidgets.QPushButton("Cancel")
        self.cancel_button.clicked.connect(self._cancel)
        for button in (self.apply_button, self.reset_button, self.cancel_button):
            button.setAutoDefault(False)
            button.setDefault(False)
        button_row.addWidget(self.apply_button)
        button_row.addWidget(self.reset_button)
        button_row.addWidget(self.cancel_button)
        layout.addLayout(button_row)

        self.status_label = QtWidgets.QLabel("")
        self.status_label.setWordWrap(True)
        layout.addWidget(self.status_label)
        self._set_edit_enabled(False)

    def _mode(self) -> str:
        return normalize_mode(self.mode_combo.currentText())

    def _update_mode_controls(self, *_args) -> None:
        if not self._building and (self.picks or self._wizard is not None):
            self._clear_pick_state()
        self._building = True
        mode = self._mode()
        self.side_combo.clear()
        for label, value in side_labels(mode):
            self.side_combo.addItem(label, value)
        self.side_combo.setCurrentIndex(1)
        if mode == "bond":
            self.slider.setRange(10, 500)
        elif mode == "angle":
            self.slider.setRange(1, 1790)
        else:
            self.slider.setRange(-1800, 1800)
        self._building = False

    def _clear_pick_state(self) -> None:
        self._clear_active_wizard()
        self.picks.clear()
        self.snapshot = None
        self.picked_list.clear()
        self.current_label.setText("-")
        self.target_edit.clear()
        self._set_edit_enabled(False)

    def _set_edit_enabled(self, enabled: bool) -> None:
        for widget in (
            self.target_edit,
            self.slider,
            self.apply_button,
            self.reset_button,
            self.cancel_button,
            self.side_combo,
            self.granularity_combo,
        ):
            widget.setEnabled(enabled)

    def _start_picking(self) -> None:
        self._clear_pick_state()
        wizard = ICPickWizard(mode=self._mode(), on_complete=self.on_atoms_picked, cmd_module=cmd)
        self._wizard = wizard
        set_active_wizard(cmd, wizard)
        self.status_label.setText(prompt_for_pick(self._mode(), 0))

    def on_atoms_picked(self, picks: list[AtomRef]) -> None:
        self.picks = picks
        self.picked_list.clear()
        for idx, pick in enumerate(picks, 1):
            self.picked_list.addItem(f"{idx}. {pick.label()}")
        self.snapshot = CoordinateSnapshot.capture(cmd, picks[0].object, picks[0].state)
        value = self._current_value()
        self.current_label.setText(self._format_value(value))
        self._set_target_value(value)
        self._set_edit_enabled(True)
        self.status_label.setText("Ready")

    def _coords(self) -> list[tuple[float, float, float]]:
        coords = []
        for pick in self.picks:
            coords.append(tuple(float(v) for v in cmd.get_atom_coords(pick.selection, state=pick.state)))
        return coords

    def _current_value(self) -> float:
        coords = self._coords()
        mode = self._mode()
        if mode == "bond":
            return distance(coords[0], coords[1])
        if mode == "angle":
            return angle(coords[0], coords[1], coords[2])
        return dihedral(coords[0], coords[1], coords[2], coords[3])

    def _format_value(self, value: float) -> str:
        unit = "A" if self._mode() == "bond" else "deg"
        return f"{value:.6f} {unit}"

    def _set_target_value(self, value: float) -> None:
        self._building = True
        self.target_edit.setText(f"{value:.6f}")
        if self._mode() == "bond":
            self.slider.setValue(int(round(value * 100.0)))
        else:
            self.slider.setValue(int(round(value * 10.0)))
        self._building = False

    def _target_value(self) -> float:
        return float(self.target_edit.text())

    def _slider_changed(self, value: int) -> None:
        if self._building:
            return
        target = value / 100.0 if self._mode() == "bond" else value / 10.0
        self.target_edit.setText(f"{target:.6f}")
        if self.auto_preview.isChecked():
            self._preview()

    def _target_edit_finished(self) -> None:
        if self._building:
            return
        try:
            self._set_target_value(self._target_value())
            if self.auto_preview.isChecked():
                self._preview()
        except Exception as exc:
            self.status_label.setText(self._error_text(exc))

    def _command(self) -> Callable:
        return {
            "bond": ic_bond,
            "angle": ic_angle,
            "dihedral": ic_dihedral,
        }[self._mode()]

    def _selection_args(self) -> list[str]:
        return [pick.selection for pick in self.picks]

    def _run_edit(self) -> None:
        if not self.picks:
            return
        if self.snapshot is not None:
            self.snapshot.restore(cmd)
        command = self._command()
        command(
            *self._selection_args(),
            self._target_value(),
            side=self.side_combo.currentData(),
            mode=self.granularity_combo.currentData(),
            state=self.picks[0].state,
            object_mode="inplace",
            preview=0,
            cmd_module=cmd,
        )
        value = self._current_value()
        self.current_label.setText(self._format_value(value))
        self._set_target_value(value)

    def _preview(self) -> None:
        if not self.picks:
            self.status_label.setText("Pick atoms before previewing.")
            return
        try:
            self._run_edit()
            self.status_label.setText("Preview updated")
        except Exception as exc:
            self.status_label.setText(self._error_text(exc))

    def _apply(self) -> None:
        try:
            self._run_edit()
            if self.picks:
                self.snapshot = CoordinateSnapshot.capture(cmd, self.picks[0].object, self.picks[0].state)
            self.status_label.setText("Applied")
        except Exception as exc:
            self.status_label.setText(self._error_text(exc))

    def _reset(self) -> None:
        if self.snapshot is not None:
            self.snapshot.restore(cmd)
            value = self._current_value()
            self.current_label.setText(self._format_value(value))
            self._set_target_value(value)
            self.status_label.setText("Reset")

    def _cancel(self) -> None:
        self._reset()
        self._clear_active_wizard()
        self.close()

    def closeEvent(self, event) -> None:  # noqa: N802 - Qt API name
        self._clear_active_wizard()
        event.accept()

    def _clear_active_wizard(self) -> None:
        if self._wizard is None:
            return
        wizard = self._wizard
        try:
            wizard.cleanup()
        except Exception:
            pass
        try:
            if cmd.get_wizard() is wizard:
                if hasattr(cmd, "set_wizard"):
                    cmd.set_wizard()
                else:
                    cmd.wizard(None)
        except Exception:
            pass
        self._wizard = None

    def _error_text(self, exc: Exception) -> str:
        message = str(exc).strip()
        if message:
            return message
        return exc.__class__.__name__


# ---------------------------------------------------------------------------
# Plugin entry points and safe auto-registration

__all__ = [
    "__init_plugin__",
    "launch",
    "open_dialog",
    "register_commands",
    "ic_bond",
    "ic_angle",
    "ic_dihedral",
    "ic_pick",
]


def launch() -> None:
    """Open the Qt editor panel inside PyMOL."""
    open_dialog()


def __init_plugin__(app=None) -> None:
    """PyMOL plugin entry point."""
    register_commands()

    try:
        from pymol.plugins import addmenuitemqt
    except Exception:
        return

    addmenuitemqt("Internal Coordinate Editor", launch)


def _try_auto_register_commands() -> None:
    """Register commands when this file is executed with PyMOL's run command."""
    if _pymol_cmd is None:
        return
    try:
        register_commands(_pymol_cmd)
    except Exception as exc:
        print(f"Internal Coordinate Editor: command registration failed: {exc}")


_try_auto_register_commands()
