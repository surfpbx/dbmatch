import numpy as np
from OgreInterface.miller import MillerSearch as _MillerSearch
from OgreInterface.surfaces import OrientedBulk
from OgreInterface.lattice_match import ZurMcGill

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
                    n_flm_reps = np.abs(np.linalg.det(sub_trasf[:2, :2]))

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
