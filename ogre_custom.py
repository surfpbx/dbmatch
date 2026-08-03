"""Two kinds of customization over OgreInterface's own matching machinery:

- `MillerSearch` subclasses `OgreInterface.miller.MillerSearch`, forcing the
  scan over the fixed `HEX_INDS` list (not just the miller indices implied by
  the substrate's own symmetry) and computing the `n_sub_reps`/`n_flm_reps`/
  transform strings that match.py writes to each matches-db row.
- `interface_from_row` (and its helpers) rebuild the exact `Interface` a
  matches-db row already describes directly from that row's own
  sub_transform/flm_transform, instead of re-deriving it via
  `InterfaceGenerator`'s own `ZurMcGill` search (`generate_interface(index=0)`),
  which depends on the search's tolerances happening to include the known
  match and on `ZurMcGill.run()`'s sort order happening to put it first.

  `ZurMcGill.run()` always sorts its results by (area, strain) ascending, and
  `MillerSearch.run_scan` above took exactly that first (smallest) match when
  it wrote each matches-db row -- so index 0 already reproduces it in
  practice. `match_from_transforms` instead reconstructs the exact
  `OgreMatch` in closed form from the stored transform, with no search and no
  tolerance dependence at all: every field OgreInterface actually reads off
  `match` to build atomic positions (*_sl_transform, *_sl_basis,
  *_sl_scale_factors, *_align_transform) is a deterministic function of the
  substrate/film `Surface`'s own inplane_vectors/crystallographic_basis and
  the known transform. area/strain are display-only (`Interface.area` is its
  own computed property; `match.strain` only feeds the printed summary
  string) so they're just taken from the matches-db row, cross-checked
  against a recomputed strain as a sanity assertion.
"""

import math

import numpy as np
from OgreInterface.lattice_match import OgreMatch, ZurMcGill
from OgreInterface.miller import MillerSearch as _MillerSearch
from OgreInterface.surfaces import OrientedBulk

HEX_INDS = [
    [0, 1, 0],
    [1, 1, 0],
    [0, 0, 1],
    [1, 0, 1],
    [1, 1, 1],
]


def transform_to_str(array):
    coefs_2D = array[:2, :2].flatten()
    return ','.join(str(i) for i in coefs_2D.flatten())


def transform_from_str(s):
    """Inverse of transform_to_str: 'a,b,c,d' -> [[a,b],[c,d]]."""
    a, b, c, d = (int(x) for x in s.split(','))
    return np.array([[a, b], [c, d]])


def hkl_from_str(hkl):
    """Inverse of the ','.join(...) used for sub_hkl/flm_hkl: 'h,k,l' -> [h,k,l]."""
    return [int(i) for i in hkl.split(',')]


class MillerSearch(_MillerSearch):
    def run_scan(self):
        substrates = []
        films = []

        gamma = self.film.lattice.gamma

        if np.isclose(gamma, 120.0, atol=0.1, rtol=0.0):
            film_inds = HEX_INDS
        else:
            film_inds = np.unique(np.abs(self.film_inds), axis=0)

        for inds in HEX_INDS:
            sub_obs = OrientedBulk(
                bulk=self.substrate,
                miller_index=inds,
                make_planar=True,
            )
            sub_inplane_vectors = sub_obs.inplane_vectors
            sub_basis = sub_obs.crystallographic_basis
            sub_miller = inds

            substrates.append([sub_inplane_vectors, sub_miller, sub_basis])

        for inds in film_inds:
            film_obs = OrientedBulk(
                bulk=self.film,
                miller_index=inds,
                make_planar=True,
            )
            film_inplane_vectors = film_obs.inplane_vectors
            film_basis = film_obs.crystallographic_basis
            film_miller = inds

            films.append([film_inplane_vectors, film_miller, film_basis])

        results = []

        for substrate in substrates:
            for film in films:
                zm = ZurMcGill(
                    film_vectors=film[0],
                    substrate_vectors=substrate[0],
                    film_basis=film[2],
                    substrate_basis=substrate[2],
                    max_area=self.max_area,
                    max_strain=self.max_strain,
                    max_area_mismatch=self.max_area_mismatch,
                    max_area_scale_factor=self.max_area_scale_factor,
                )
                matches = zm.run()

                if len(matches) > 0:
                    min_area_match = matches[0]
                    sub_trasf = min_area_match.substrate_sl_transform
                    flm_trasf = min_area_match.film_sl_transform
                    n_sub_reps = np.abs(np.linalg.det(sub_trasf[:2, :2]))
                    n_flm_reps = np.abs(np.linalg.det(flm_trasf[:2, :2]))

                    results.append({
                        'area'         : min_area_match.area,
                        'strain'       : min_area_match.strain,
                        'sub_hkl'      : ','.join(str(i) for i in substrate[1]),
                        'flm_hkl'      : ','.join(str(i) for i in film[1]),
                        'sub_transform': transform_to_str(sub_trasf),
                        'flm_transform': transform_to_str(flm_trasf),
                        'n_sub_reps'   : int(n_sub_reps.round(0)),
                        'n_flm_reps'   : int(n_flm_reps.round(0)),
                    })

        return results


def _embed_3x3(transform_2d):
    full = np.eye(3, dtype=int)
    full[:2, :2] = transform_2d
    return full


def _align_to_x_axis(vectors):
    """3x3 rotation putting vectors[0] on the +x axis -- the single-pair
    version of ZurMcGill._build_a_to_i (which operates on a batch)."""
    ax, ay = vectors[0, :2] / np.linalg.norm(vectors[0, :2])
    return np.array([
        [ax, -ay, 0.0],
        [ay, ax, 0.0],
        [0.0, 0.0, 1.0],
    ])


def _inv2x2(m):
    a, b = m[0]
    c, d = m[1]
    det = a * d - b * c
    return (1.0 / det) * np.array([[d, -b], [-c, a]])


def match_from_transforms(sub_surface, film_surface, sub_transform, film_transform, area, strain):
    """
    Build a real OgreMatch straight from a matches-db row's own
    sub_transform/film_transform (already-parsed 2x2 int arrays, e.g. via
    transform_from_str), applied to sub_surface/film_surface -- the Surface
    objects (e.g. subs[0]/films[0]) that will be passed to InterfaceGenerator.

    area/strain are the matches-db row's own values; used only for display
    (see module docstring), not for building atomic positions.
    """
    film_vectors = film_surface.inplane_vectors
    sub_vectors = sub_surface.inplane_vectors
    film_basis = film_surface.crystallographic_basis
    sub_basis = sub_surface.crystallographic_basis

    film_sl_vectors = film_transform.astype(float) @ film_vectors
    sub_sl_vectors = sub_transform.astype(float) @ sub_vectors

    # _get_sl_basis is a pure function of the transform + film/substrate
    # basis -- no tolerances involved, so it's safe to call directly with
    # placeholder max_strain/max_area/max_area_mismatch (unused: we never
    # call zm.run()/_is_same).
    zm = ZurMcGill(
        film_vectors=film_vectors,
        substrate_vectors=sub_vectors,
        film_basis=film_basis,
        substrate_basis=sub_basis,
        max_area=1.0,
        max_strain=1.0,
        max_area_mismatch=1.0,
    )
    film_sl_basis, film_sl_scale_factors, sub_sl_basis, sub_sl_scale_factors = (
        zm._get_sl_basis(film_transform[None], sub_transform[None])
    )

    film_a_to_i = _align_to_x_axis(film_sl_vectors)
    sub_a_to_i = _align_to_x_axis(sub_sl_vectors)
    aligned_film = film_sl_vectors @ film_a_to_i
    aligned_sub = sub_sl_vectors @ sub_a_to_i

    film_area = np.linalg.norm(np.cross(film_sl_vectors[0], film_sl_vectors[1]))
    sub_area = np.linalg.norm(np.cross(sub_sl_vectors[0], sub_sl_vectors[1]))

    raw_strain_transform = _inv2x2(aligned_film[:, :2]) @ aligned_sub[:, :2]
    if film_area > sub_area:
        strain_transform_2d = _inv2x2(raw_strain_transform)
    else:
        strain_transform_2d = raw_strain_transform

    computed_strain = (1 / math.sqrt(2)) * np.linalg.norm(np.eye(2) - strain_transform_2d)
    assert math.isclose(computed_strain, strain, rel_tol=1e-3, abs_tol=1e-6), (
        f'recomputed strain {computed_strain} does not match stored strain {strain} '
        f'-- sub_transform/film_transform likely do not correspond to this '
        f'sub_surface/film_surface'
    )

    strain_transform = np.eye(3)
    strain_transform[:2, :2] = strain_transform_2d

    return OgreMatch(
        area=area,
        strain=strain,
        film_vectors=film_vectors,
        film_sl_vectors=film_sl_vectors,
        film_zur_mcgill_transform=film_transform,
        film_sl_transform=_embed_3x3(film_transform),
        substrate_vectors=sub_vectors,
        substrate_sl_vectors=sub_sl_vectors,
        substrate_zur_mcgill_transform=sub_transform,
        substrate_sl_transform=_embed_3x3(sub_transform),
        substrate_basis=sub_basis,
        substrate_sl_basis=sub_sl_basis[0],
        substrate_sl_scale_factors=sub_sl_scale_factors[0],
        film_basis=film_basis,
        film_sl_basis=film_sl_basis[0],
        film_sl_scale_factors=film_sl_scale_factors[0],
        substrate_align_transform=sub_a_to_i,
        film_align_transform=film_a_to_i,
        film_to_substrate_strain_transform=strain_transform,
    )


def interface_from_row(interface_generator, row, sub_surface, film_surface):
    """
    Build the single Interface corresponding to row's own sub_transform/
    flm_transform, using interface_generator (already constructed with
    substrate=sub_surface, film=film_surface) -- bypassing its match_list
    (and so its ZurMcGill search's tolerances/sort order) entirely.
    """
    match = match_from_transforms(
        sub_surface=sub_surface,
        film_surface=film_surface,
        sub_transform=transform_from_str(row.sub_transform),
        film_transform=transform_from_str(row.flm_transform),
        area=row.area,
        strain=row.strain,
    )
    return interface_generator._build_interface(match=match)
