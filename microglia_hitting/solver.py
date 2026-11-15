"""3D signal-weighted hitting-probability solvers.

The solver checks all 26 neighbors of each voxel. Directional behavior is
controlled by ``mode``:

``standard``
    Allow every in-bounds neighbor.
``inner``
    Allow a neighbor when its squared XY distance to the soma is no greater
    than the current squared distance plus ``inner_squared_tolerance``.
``hook_friendly``
    Allow a neighbor when its XY distance to the soma is no greater than the
    current distance plus ``outward_tolerance_xy``.

Allowed inactive neighbors remain in the transition denominator but contribute
zero hitting probability. Thus, the outside of the active domain acts as a
zero-probability boundary.
"""

from __future__ import annotations

from time import perf_counter
from typing import Any

import numba as nb
import numpy as np


_MODE_STANDARD = 0
_MODE_INNER = 1
_MODE_HOOK_FRIENDLY = 2

_MODE_CODES = {
    "standard": _MODE_STANDARD,
    "inner": _MODE_INNER,
    "hook_friendly": _MODE_HOOK_FRIENDLY,
}


@nb.njit(inline="always")
def _transition_is_allowed(
    x: int,
    y: int,
    neighbor_x: int,
    neighbor_y: int,
    soma_x: float,
    soma_y: float,
    mode_code: int,
    outward_tolerance_xy: float,
    inner_squared_tolerance: float,
) -> bool:
    """Return whether one XY neighbor satisfies the selected direction rule."""

    if mode_code == _MODE_STANDARD:
        return True

    current_distance_squared = (
        (x - soma_x) ** 2
        + (y - soma_y) ** 2
    )
    neighbor_distance_squared = (
        (neighbor_x - soma_x) ** 2
        + (neighbor_y - soma_y) ** 2
    )

    if mode_code == _MODE_INNER:
        return (
            neighbor_distance_squared
            <= current_distance_squared + inner_squared_tolerance
        )

    current_distance = np.sqrt(current_distance_squared)
    neighbor_distance = np.sqrt(neighbor_distance_squared)
    return neighbor_distance <= current_distance + outward_tolerance_xy


@nb.njit(parallel=True, cache=True)
def _compute_denominators(
    weights: np.ndarray,
    coordinates: np.ndarray,
    soma_x: float,
    soma_y: float,
    mode_code: int,
    outward_tolerance_xy: float,
    inner_squared_tolerance: float,
) -> np.ndarray:
    """Calculate outgoing transition-weight denominators."""

    z_size, y_size, x_size = weights.shape
    denominators = np.zeros(coordinates.shape[0], dtype=np.float64)

    for index in nb.prange(coordinates.shape[0]):
        z, y, x = coordinates[index]
        total_weight = 0.0

        for dz in range(-1, 2):
            for dy in range(-1, 2):
                for dx in range(-1, 2):
                    if dz == 0 and dy == 0 and dx == 0:
                        continue

                    neighbor_z = z + dz
                    neighbor_y = y + dy
                    neighbor_x = x + dx

                    if not (
                        0 <= neighbor_z < z_size
                        and 0 <= neighbor_y < y_size
                        and 0 <= neighbor_x < x_size
                    ):
                        continue

                    if not _transition_is_allowed(
                        x,
                        y,
                        neighbor_x,
                        neighbor_y,
                        soma_x,
                        soma_y,
                        mode_code,
                        outward_tolerance_xy,
                        inner_squared_tolerance,
                    ):
                        continue

                    total_weight += weights[
                        neighbor_z,
                        neighbor_y,
                        neighbor_x,
                    ]

        denominators[index] = total_weight

    return denominators


@nb.njit(parallel=True, cache=True)
def _update_once(
    weights: np.ndarray,
    active: np.ndarray,
    coordinates: np.ndarray,
    denominators: np.ndarray,
    soma_x: float,
    soma_y: float,
    mode_code: int,
    outward_tolerance_xy: float,
    inner_squared_tolerance: float,
    previous: np.ndarray,
    updated: np.ndarray,
) -> None:
    """Perform one Jacobi update of all unknown active voxels."""

    z_size, y_size, x_size = weights.shape

    for index in nb.prange(coordinates.shape[0]):
        z, y, x = coordinates[index]
        denominator = denominators[index]
        new_probability = 0.0

        if denominator > 0.0:
            for dz in range(-1, 2):
                for dy in range(-1, 2):
                    for dx in range(-1, 2):
                        if dz == 0 and dy == 0 and dx == 0:
                            continue

                        neighbor_z = z + dz
                        neighbor_y = y + dy
                        neighbor_x = x + dx

                        if not (
                            0 <= neighbor_z < z_size
                            and 0 <= neighbor_y < y_size
                            and 0 <= neighbor_x < x_size
                        ):
                            continue

                        if not _transition_is_allowed(
                            x,
                            y,
                            neighbor_x,
                            neighbor_y,
                            soma_x,
                            soma_y,
                            mode_code,
                            outward_tolerance_xy,
                            inner_squared_tolerance,
                        ):
                            continue

                        if active[neighbor_z, neighbor_y, neighbor_x]:
                            new_probability += (
                                weights[neighbor_z, neighbor_y, neighbor_x]
                                * previous[neighbor_z, neighbor_y, neighbor_x]
                                / denominator
                            )

        updated[z, y, x] = new_probability


def solve_hitting_probability(
    weights: np.ndarray,
    active: np.ndarray,
    target: np.ndarray,
    soma_x: float,
    soma_y: float,
    *,
    zero_target: np.ndarray | None = None,
    mode: str = "hook_friendly",
    outward_tolerance_xy: float = 1.0,
    inner_squared_tolerance: float = 3.0,
    maximum_iterations: int = 40_000,
    convergence_tolerance: float = 1e-5,
    convergence_check_interval: int = 10,
    progress_interval: int = 500,
    verbose: bool = True,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Solve a 3D signal-weighted soma hitting probability.

    Parameters
    ----------
    weights
        Nonnegative transition weights in Z-Y-X order. For the current
        workflow, pass ``MG_probability ** 5 + epsilon``.
    active
        Boolean active-domain mask in Z-Y-X order.
    target
        Boolean absorbing soma mask. Target voxels are assigned probability 1.
    zero_target
        Optional competing absorbing mask assigned probability 0. This is used
        for a connected two-cell problem: when solving Cell A, pass Cell B's
        soma here, and vice versa.
    soma_x, soma_y
        Soma center in local X-Y coordinates.
    mode
        One of ``"standard"``, ``"inner"``, or ``"hook_friendly"``.
    outward_tolerance_xy
        Actual XY-distance tolerance used only by ``hook_friendly``.
    inner_squared_tolerance
        Squared XY-distance tolerance used only by ``inner``.

    Returns
    -------
    probability, information
        A float32 Z-Y-X probability array and a convergence summary.
    """

    weights = np.asarray(weights, dtype=np.float64)
    active = np.asarray(active, dtype=bool)
    target = np.asarray(target, dtype=bool)
    if zero_target is None:
        zero_target = np.zeros(target.shape, dtype=bool)
    else:
        zero_target = np.asarray(zero_target, dtype=bool)

    if mode not in _MODE_CODES:
        choices = ", ".join(repr(name) for name in _MODE_CODES)
        raise ValueError(f"Unknown mode {mode!r}; expected one of {choices}.")
    if weights.ndim != 3:
        raise ValueError("weights must be a three-dimensional Z-Y-X array.")
    if (
        active.shape != weights.shape
        or target.shape != weights.shape
        or zero_target.shape != weights.shape
    ):
        raise ValueError(
            "weights, active, target, and zero_target must have identical shapes."
        )
    if not target.any():
        raise ValueError("target must contain at least one soma voxel.")
    if np.any(target & zero_target):
        raise ValueError("target and zero_target must not overlap.")
    if not np.isfinite(weights).all() or np.any(weights < 0):
        raise ValueError("weights must be finite and nonnegative.")
    if outward_tolerance_xy < 0 or inner_squared_tolerance < 0:
        raise ValueError("Directional tolerances must be nonnegative.")
    if maximum_iterations < 1:
        raise ValueError("maximum_iterations must be at least 1.")
    if convergence_tolerance <= 0:
        raise ValueError("convergence_tolerance must be positive.")
    if convergence_check_interval < 1 or progress_interval < 1:
        raise ValueError("Iteration intervals must be positive integers.")

    mode_code = _MODE_CODES[mode]
    active_with_target = active | target | zero_target
    unknown = active_with_target & ~target & ~zero_target
    coordinates = np.argwhere(unknown).astype(np.int32)
    start_time = perf_counter()

    if verbose:
        print(f"Hitting-probability mode: {mode}")
        if mode == "hook_friendly":
            print(f"Outward XY tolerance: {outward_tolerance_xy:g} voxel(s)")
        elif mode == "inner":
            print(f"Squared-distance tolerance: {inner_squared_tolerance:g}")
        print(f"Active voxels: {int(active_with_target.sum()):,}")
        print(f"Unknown voxels: {int(unknown.sum()):,}")
        print(f"Target voxels: {int(target.sum()):,}")
        print(f"Zero-target voxels: {int(zero_target.sum()):,}")

    previous = np.zeros(weights.shape, dtype=np.float64)
    previous[target] = 1.0

    if coordinates.shape[0] == 0:
        information = {
            "mode": mode,
            "converged": True,
            "iterations": 0,
            "final_maximum_change": 0.0,
            "runtime_minutes": (perf_counter() - start_time) / 60.0,
            "active_voxels": int(active_with_target.sum()),
            "unknown_voxels": 0,
            "target_voxels": int(target.sum()),
            "zero_target_voxels": int(zero_target.sum()),
            "outward_tolerance_xy": float(outward_tolerance_xy),
            "inner_squared_tolerance": float(inner_squared_tolerance),
        }
        return previous.astype(np.float32), information

    denominators = _compute_denominators(
        weights,
        coordinates,
        float(soma_x),
        float(soma_y),
        mode_code,
        float(outward_tolerance_xy),
        float(inner_squared_tolerance),
    )

    updated = previous.copy()
    maximum_change = np.inf
    converged = False
    completed_iterations = maximum_iterations

    for iteration in range(1, maximum_iterations + 1):
        _update_once(
            weights,
            active_with_target,
            coordinates,
            denominators,
            float(soma_x),
            float(soma_y),
            mode_code,
            float(outward_tolerance_xy),
            float(inner_squared_tolerance),
            previous,
            updated,
        )
        updated[target] = 1.0

        if iteration % convergence_check_interval == 0:
            maximum_change = float(
                np.max(np.abs(updated[unknown] - previous[unknown]))
            )

            if verbose and (
                iteration % progress_interval == 0
                or maximum_change < convergence_tolerance
            ):
                elapsed_minutes = (perf_counter() - start_time) / 60.0
                print(
                    f"iteration={iteration:,}; "
                    f"max_change={maximum_change:.3e}; "
                    f"elapsed={elapsed_minutes:.2f} min"
                )

            if maximum_change < convergence_tolerance:
                previous, updated = updated, previous
                completed_iterations = iteration
                converged = True
                break

        previous, updated = updated, previous

    runtime_minutes = (perf_counter() - start_time) / 60.0
    probability = np.clip(previous, 0.0, 1.0).astype(np.float32)

    information = {
        "mode": mode,
        "converged": converged,
        "iterations": int(completed_iterations),
        "final_maximum_change": float(maximum_change),
        "runtime_minutes": float(runtime_minutes),
        "active_voxels": int(active_with_target.sum()),
        "unknown_voxels": int(unknown.sum()),
        "target_voxels": int(target.sum()),
        "zero_target_voxels": int(zero_target.sum()),
        "outward_tolerance_xy": float(outward_tolerance_xy),
        "inner_squared_tolerance": float(inner_squared_tolerance),
        "convergence_tolerance": float(convergence_tolerance),
        "maximum_iterations": int(maximum_iterations),
    }

    if verbose:
        print(f"Converged: {converged}")
        print(f"Completed iterations: {completed_iterations:,}")
        print(f"Runtime: {runtime_minutes:.2f} minutes")

    return probability, information
