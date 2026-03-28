#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
r"""
Leaf Area Index (LAI)
*********************

Derive the LAI for each of the corresponding functional types of the land cover
type fraction ancillary, that is the vegetative classes.

The following steps are taken:

* The total LAI source is regridded onto the grid of the land cover type
  fraction grid.
* The total LAI where missing across some but not all months is set to 1e-10.
* Any remaining unmasked LAI which is less than 1e-10 is set to 1e-10.
* At each grid point, calculate LAI corresponding to each of our functional
  types through substitution, where the total vegetation across functional
  types is non-zero.
* Where the total vegetation fraction (for the functional types) equals less
  than 20% of a grid box but non-zero, each LAI PFT is masked.  This avoids
  unreasonably hight LAI PFTs.
* Make the resulting LAI fields consistent with the land cover type fraction
  by using the spiral search (nearest neighbour).

The Leaf Area Index (LAI) of each plant functional type (PFT) is related to
the total LAI through:

.. math::

    L_{tot} & = & \sum_{i} L_i f_i \\

Where :math:`L_{tot}` is the total LAI, :math:`L_i` is the LAI contribution
from PFT :math:`i`, and :math:`f_i` is the fraction of land covered by PFT
:math:`i`.

The LAI of each PFT is found using the assumption that the relative LAIs of
each plant functional type are fixed.  If the LAI of a PFT is expressed
relative to the first PFT then this is equivalent to:

.. math::

  \begin{eqnarray}
    L_2 & = & \gamma_2 L_1 \\
    L_3 & = & \gamma_3 L_1 \\
        & \vdots & \\
  \end{eqnarray}

Where :math:`\gamma_i` is the LAI of type :math:`i` relative to PFT 1.

Substitute this second set of equations into the first gives:

.. math::

  L_{tot} & = & \sum_i \gamma_i L_1 f_i \\

(where :math:`\gamma_1 = 1`)

So

.. math::

  \begin{eqnarray}
    L_{tot} & = & L_1 \sum_i \gamma_i  f_i \\
    L_1 & = & L_{tot} / \sum_i \gamma_i  f_i \\
  \end{eqnarray}

Then substitute :math:`L_1` back into the expressions for the individual
:math:`L_i`.

.. math::

   L_i = \gamma_i L_{tot} / \sum_i \gamma_i f_i \\

This is calculated at each grid point.

"""
from functools import partial

import ants
import ants.fileformats.json
import ants.decomposition as decomp
import ants.io.save as save
import iris
import numpy as np

import lai


def load_data(sources, relative_weights_path):
    lai_source = ants.load_cube(sources, 'leaf_area_index')
    lct_stash_con = iris.AttributeConstraint(STASH='m01s00i216')
    lct = ants.load_cube(sources, lct_stash_con)
    relative_weights = ants.fileformats.json.load(
        relative_weights_path, keys=['relative_weights', 'jules_classes'],
        dtypes=['float', 'int'])
    if (relative_weights['relative_weights'].size !=
            relative_weights['jules_classes'].size):
        raise ValueError('Expecting lai weights for each corresponding '
                         'jules class specified.')
    return lai_source, lct, relative_weights


def main(sources, relative_weights_path, output_filepath, use_new_saver):
    lai_source, lct, relative_weights = load_data(sources,
                                                       relative_weights_path)
    operation = partial(lai.calc_lai, relative_weights)
    lai_tar = decomp.decompose(operation, lai_source, lct)

    # Final Nearest neighbour fill to ensure the lsm is consistent with the
    # lct.
    lbm = lct[0].copy(
        np.ma.getmaskarray(lct[0].data))
    nfiller = ants.analysis.FillMissingPoints(lai_tar[0], target_mask=lbm)
    for time_step in range(lai_tar.shape[0]):
        # Iterate over each time step - until such a time that horizontal grid
        # re-order has 4D+ support.
        # We copy the cube and slice the array to make it a view of the
        # original cubes array.

        # move the nfiller inside so it can deal with the missing points difference for differnt time
        
        #nfiller = ants.analysis.FillMissingPoints(lai_tar[time_step], target_mask=lbm)
        lai_slice = lai_tar[time_step].copy(lai_tar.data[time_step])
        nfiller(lai_slice)

    
    save.ancil(lai_tar, output_filepath)
    if use_new_saver:
        save.netcdf(lai_tar, output_filepath)

    return lai_tar


if __name__ == '__main__':
    parser = ants.AntsArgParser()
    parser.add_argument('--relative-weights', type=str,
                        help='Relative weights for LAI funtional types.',
                        required=True)
    args = parser.parse_args()
    main(args.sources, args.relative_weights, args.output, args.use_new_saver)
