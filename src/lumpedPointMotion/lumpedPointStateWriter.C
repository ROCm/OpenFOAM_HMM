/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "lumpedPointState.H"
#include "OFstream.H"
#include "axesRotation.H"
#include "coordinateSystem.H"
#include "foamVtkOutput.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// file-local
const static Foam::FixedList<Foam::point, 4> standardCorners
{
    {-0.1, -0.1, 0},
    {+0.1, -0.1, 0},
    {+0.1, +0.1, 0},
    {-0.1, +0.1, 0}
};


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::lumpedPointState::writeVTP
(
    const fileName& outputFile,
    const vector& axis
) const
{
    // local-to-global transformation tensors
    const tensorField& localToGlobal = rotations();

    OFstream fos(outputFile);
    std::ostream& os = fos.stdStream();

    autoPtr<vtk::formatter> format = vtk::newFormatter
    (
        os,
        vtk::formatType::INLINE_ASCII
    );

    format().xmlHeader()
        .beginVTKFile<vtk::fileTag::POLY_DATA>();

    //
    // The 'spine' of lumped mass points
    //
    {
        format()
            .tag
            (
                vtk::fileTag::PIECE,
                vtk::fileAttr::NUMBER_OF_POINTS, points_.size(),
                vtk::fileAttr::NUMBER_OF_VERTS,  points_.size(),
                vtk::fileAttr::NUMBER_OF_LINES,  1
            );

        // 'points'
        {
            const uint64_t payLoad = vtk::sizeofData<float, 3>(points_.size());

            format()
                .tag(vtk::fileTag::POINTS)
                .beginDataArray<float, 3>(vtk::dataArrayAttr::POINTS);

            format().writeSize(payLoad);
            vtk::writeList(format(), points_);
            format().flush();

            format()
                .endDataArray()
                .endTag(vtk::fileTag::POINTS);
        }

        // <Verts>
        format().tag(vtk::fileTag::VERTS);

        //
        // 'connectivity'
        //
        {
            const uint64_t payLoad = vtk::sizeofData<label>(points_.size());

            format().beginDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY);
            format().writeSize(payLoad);

            forAll(points_, i)
            {
                format().write(i);
            }
            format().flush();

            format().endDataArray();
        }

        //
        // 'offsets'  (connectivity offsets)
        // = linear mapping onto points (with 1 offset)
        //
        {
            const uint64_t payLoad = vtk::sizeofData<label>(points_.size());

            format().beginDataArray<label>(vtk::dataArrayAttr::OFFSETS);
            format().writeSize(payLoad);

            forAll(points_, i)
            {
                format().write(i+1);
            }
            format().flush();

            format().endDataArray();
        }

        format().endTag(vtk::fileTag::VERTS);
        // </Verts>


        // <Lines>
        format().tag(vtk::fileTag::LINES);

        //
        // 'connectivity'
        //
        {
            const uint64_t payLoad = vtk::sizeofData<label>(points_.size());

            format().beginDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY);
            format().writeSize(payLoad);

            forAll(points_, i)
            {
                format().write(i);
            }
            format().flush();

            format().endDataArray();
        }

        //
        // 'offsets'  (connectivity offsets)
        // = single line
        //
        {
            const uint64_t payLoad = vtk::sizeofData<label>(1);

            format().beginDataArray<label>(vtk::dataArrayAttr::OFFSETS);
            format().writeSize(payLoad);

            format().write(points_.size());
            format().flush();

            format().endDataArray();
        }

        format().endTag(vtk::fileTag::LINES);
        format().endPiece();
    }

    // Standard corners in local axis
    FixedList<point, 4> corners;

    {
        coordinateRotations::axes orient(axis);
        coordinateSystem cornerTransform(orient);

        forAll(standardCorners, corni)
        {
            corners[corni] = cornerTransform.transform(standardCorners[corni]);
        }
    }

    //
    // Planes to visualize location/rotation
    //
    {
        const label nPoints = 4*points_.size();  // 4 points per quad
        const label nPolys  = points_.size();

        format()
            .tag
            (
                vtk::fileTag::PIECE,
                vtk::fileAttr::NUMBER_OF_POINTS, nPoints,
                vtk::fileAttr::NUMBER_OF_POLYS,  nPolys
            );

        // 'points'
        {
            const uint64_t payLoad = vtk::sizeofData<float, 3>(nPoints);

            format()
                .tag(vtk::fileTag::POINTS)
                .beginDataArray<float, 3>(vtk::dataArrayAttr::POINTS);

            format().writeSize(payLoad);

            forAll(points_, posI)
            {
                const point& origin = points_[posI];
                const tensor& rotTensor =
                (
                    posI < localToGlobal.size()
                  ? localToGlobal[posI]
                  : pTraits<tensor>::I
                );


                for (const point& cornerPt : corners)
                {
                    // Local-to-global rotation and translation
                    const point pt = (rotTensor & cornerPt) + origin;

                    vtk::write(format(), pt);
                }
            }

            format().flush();

            format()
                .endDataArray()
                .endTag(vtk::fileTag::POINTS);
        }

        // <Polys>
        format().tag(vtk::fileTag::POLYS);

        //
        // 'connectivity' - 4 points (ie, quad)
        //
        {
            const uint64_t payLoad = vtk::sizeofData<label>(4*nPolys);

            format().beginDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY);
            format().writeSize(payLoad);

            for (label i=0; i < 4*nPolys; ++i)
            {
                format().write(i);
            }
            format().flush();

            format().endDataArray();
        }

        //
        // 'offsets'  (connectivity offsets)
        // = single quad
        //
        {
            const uint64_t payLoad = vtk::sizeofData<label>(nPolys);

            format().beginDataArray<label>(vtk::dataArrayAttr::OFFSETS);
            format().writeSize(payLoad);

            for (label i=0; i < nPolys; ++i)
            {
                const label off = 4 * (i+1);
                format().write(off);
            }
            format().flush();

            format().endDataArray();
        }

        format().endTag(vtk::fileTag::POLYS);

#if 0
        format().beginCellData();

        // zone Id
        {
            const uint64_t payLoad = vtk::sizeofData<label>(nPolys);

            format().beginDataArray<label>("zoneId");
            format().writeSize(payLoad);

            for (label i=0; i < nPolys; ++i)
            {
                format().write(i);
            }
            format().flush();

            format().endDataArray();
        }

        format().endCellData();
#endif

        format().endPiece();
    }

    // Finally
    // could add a 'ghost' level above to visualize extrapolated values
    // draw as two triangles to distingush from real levels ...

    format().endTag(vtk::fileTag::POLY_DATA)
        .endVTKFile();
}


// ************************************************************************* //
