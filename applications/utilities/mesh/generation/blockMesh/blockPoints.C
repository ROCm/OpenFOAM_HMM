/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    private member of block. Creates vertices for cells filling the block.

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "block.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void block::blockPoints()
{
    // set local variables for mesh specification
    label ni = blockDef_.n().x();
    label nj = blockDef_.n().y();
    label nk = blockDef_.n().z();

    point start = blockDef_.points()[blockDef_.blockShape()[0]];
    point xEnd = blockDef_.points()[blockDef_.blockShape()[1]];
    point xyEnd = blockDef_.points()[blockDef_.blockShape()[2]];
    point yEnd = blockDef_.points()[blockDef_.blockShape()[3]];

    point zEnd = blockDef_.points()[blockDef_.blockShape()[4]];
    point xzEnd = blockDef_.points()[blockDef_.blockShape()[5]];
    point xyzEnd = blockDef_.points()[blockDef_.blockShape()[6]];
    point yzEnd = blockDef_.points()[blockDef_.blockShape()[7]];

    // set reference to the list of edge point and weighting factors
    const List<List<point> >& edgePoints = blockDef_.blockEdgePoints();
    const scalarListList& edgeWeights = blockDef_.blockEdgeWeights();

    // generate vertices

    for (label k = 0; k <= nk; k++)
    {
        for (label j = 0; j <= nj; j++)
        {
            for (label i = 0; i <= ni; i++)
            {
                label vertexNo = vtxLabel(i, j, k);

                vector edgex1 = start*(1.0 - edgeWeights[0][i])
                      + xEnd*edgeWeights[0][i];

                vector edgex2 = yEnd*(1.0 - edgeWeights[1][i])
                      + xyEnd*edgeWeights[1][i];

                vector edgex3 = yzEnd*(1.0 - edgeWeights[2][i])
                      + xyzEnd*edgeWeights[2][i];

                vector edgex4 = zEnd*(1.0 - edgeWeights[3][i])
                      + xzEnd*edgeWeights[3][i];



                vector edgey1 = start*(1.0 - edgeWeights[4][j])
                      + yEnd*edgeWeights[4][j];

                vector edgey2 = xEnd*(1.0 - edgeWeights[5][j])
                      + xyEnd*edgeWeights[5][j];

                vector edgey3 = xzEnd*(1.0 - edgeWeights[6][j])
                      + xyzEnd*edgeWeights[6][j];

                vector edgey4 = zEnd*(1.0 - edgeWeights[7][j])
                      + yzEnd*edgeWeights[7][j];



                vector edgez1 = start*(1.0 - edgeWeights[8][k])
                      + zEnd*edgeWeights[8][k];

                vector edgez2 = xEnd*(1.0 - edgeWeights[9][k])
                      + xzEnd*edgeWeights[9][k];

                vector edgez3 = xyEnd*(1.0 - edgeWeights[10][k])
                      + xyzEnd*edgeWeights[10][k];

                vector edgez4 = yEnd*(1.0 - edgeWeights[11][k])
                      + yzEnd*edgeWeights[11][k];

                // calculate the importance factors for all edges
                // x - direction
                scalar impx1 =
                (
                    (1.0 - edgeWeights[0][i])
                        *(1.0 - edgeWeights[4][j])
                        *(1.0 - edgeWeights[8][k])
                  + edgeWeights[0][i]
                        *(1.0 - edgeWeights[5][j])
                        *(1.0 - edgeWeights[9][k])
                );

                scalar impx2 =
                (
                    (1.0 - edgeWeights[1][i])
                        *edgeWeights[4][j]
                        *(1.0 - edgeWeights[11][k])
                  + edgeWeights[1][i]
                        *edgeWeights[5][j]
                        *(1.0 - edgeWeights[10][k])
                );

               scalar impx3 =
               (
                   (1.0 - edgeWeights[2][i])
                       *edgeWeights[7][j]
                       *edgeWeights[11][k]
                 + edgeWeights[2][i]
                       *edgeWeights[6][j]
                       *edgeWeights[10][k]
               );


               scalar impx4 =
               (
                   (1.0 - edgeWeights[3][i])
                       *(1.0 - edgeWeights[7][j])
                       *edgeWeights[8][k]
                 + edgeWeights[3][i]
                       *(1.0 - edgeWeights[6][j])
                       *edgeWeights[9][k]
               );



                // y - direction
                scalar impy1 =
                (
                    (1.0 - edgeWeights[4][j])
                        *(1.0 - edgeWeights[0][i])
                        *(1.0 - edgeWeights[8][k])
                  + edgeWeights[4][j]
                        *(1.0 - edgeWeights[1][i])
                        *(1.0 - edgeWeights[11][k])
                );

                scalar impy2 =
                (
                    (1.0 - edgeWeights[5][j])
                        *edgeWeights[0][i]
                        *(1.0 - edgeWeights[9][k])
                  + edgeWeights[5][j]
                        *edgeWeights[1][i]
                        *(1.0 - edgeWeights[10][k])
                );

                scalar impy3 =
                (
                    (1.0 - edgeWeights[6][j])
                        *edgeWeights[3][i]
                        *edgeWeights[9][k]
                  + edgeWeights[6][j]
                        *edgeWeights[2][i]
                        *edgeWeights[10][k]
                );

                scalar impy4 =
                (
                    (1.0 - edgeWeights[7][j])
                        *(1.0 - edgeWeights[3][i])
                        *edgeWeights[8][k]
                  + edgeWeights[7][j]
                        *(1.0 - edgeWeights[2][i])
                        *edgeWeights[11][k]
                );



                // z - direction
                scalar impz1 =
                (
                    (1.0 - edgeWeights[8][k])
                        *(1.0 - edgeWeights[0][i])
                        *(1.0 - edgeWeights[4][j])
                  + edgeWeights[8][k]
                        *(1.0 - edgeWeights[3][i])
                        *(1.0 - edgeWeights[7][j])
                );

                scalar impz2 =
                (
                    (1.0 - edgeWeights[9][k])
                        *edgeWeights[0][i]
                        *(1.0 - edgeWeights[5][j])
                  + edgeWeights[9][k]
                        *edgeWeights[3][i]
                        *(1.0 - edgeWeights[6][j])
                );

                scalar impz3 =
                (
                    (1.0 - edgeWeights[10][k])
                        *edgeWeights[1][i]
                        *edgeWeights[5][j]
                  + edgeWeights[10][k]
                        *edgeWeights[2][i]
                        *edgeWeights[6][j]
                );

                scalar impz4 =
                (
                    (1.0 - edgeWeights[11][k])
                        *(1.0 - edgeWeights[1][i])
                        *edgeWeights[4][j]
                  + edgeWeights[11][k]
                        *(1.0 - edgeWeights[2][i])
                        *edgeWeights[7][j]
                );

                // calculate the correction vectors
                vector corx1 = impx1*(edgePoints[0][i] - edgex1);
                vector corx2 = impx2*(edgePoints[1][i] - edgex2);
                vector corx3 = impx3*(edgePoints[2][i] - edgex3);
                vector corx4 = impx4*(edgePoints[3][i] - edgex4);

                vector cory1 = impy1*(edgePoints[4][j] - edgey1);
                vector cory2 = impy2*(edgePoints[5][j] - edgey2);
                vector cory3 = impy3*(edgePoints[6][j] - edgey3);
                vector cory4 = impy4*(edgePoints[7][j] - edgey4);

                vector corz1 = impz1*(edgePoints[8][k] - edgez1);
                vector corz2 = impz2*(edgePoints[9][k] - edgez2);
                vector corz3 = impz3*(edgePoints[10][k] - edgez3);
                vector corz4 = impz4*(edgePoints[11][k] - edgez4);


                // multiply by the importance factor
                // x - direction
                edgex1 *= impx1;
                edgex2 *= impx2;
                edgex3 *= impx3;
                edgex4 *= impx4;

                // y - direction
                edgey1 *= impy1;
                edgey2 *= impy2;
                edgey3 *= impy3;
                edgey4 *= impy4;

                // z - direction
                edgez1 *= impz1;
                edgez2 *= impz2;
                edgez3 *= impz3;
                edgez4 *= impz4;


                // add the contributions
                vertices_[vertexNo] = edgex1 + edgex2 + edgex3 + edgex4;
                vertices_[vertexNo] += edgey1 + edgey2 + edgey3 + edgey4;
                vertices_[vertexNo] += edgez1 + edgez2 + edgez3 + edgez4;

                vertices_[vertexNo] /= 3.0;

                vertices_[vertexNo] += corx1 + corx2 + corx3 + corx4;
                vertices_[vertexNo] += cory1 + cory2 + cory3 + cory4;
                vertices_[vertexNo] += corz1 + corz2 + corz3 + corz4;

            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
