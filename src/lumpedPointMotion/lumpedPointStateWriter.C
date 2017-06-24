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
        .beginVTKFile(vtk::fileTag::POLY_DATA, "0.1");

    //
    // The 'spine' of lumped mass points
    //
    {
        format()
            .openTag(vtk::fileTag::PIECE)
            .xmlAttr(vtk::fileAttr::NUMBER_OF_POINTS, points_.size())
            .xmlAttr(vtk::fileAttr::NUMBER_OF_VERTS,  points_.size())
            .xmlAttr(vtk::fileAttr::NUMBER_OF_LINES,  1)
            .closeTag();

        // 'points'
        {
            const uint64_t payLoad = (points_.size()*3* sizeof(float));

            format()
                .tag(vtk::fileTag::POINTS)
                .openDataArray<float, 3>(vtk::dataArrayAttr::POINTS)
                .closeTag();

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
            const uint64_t payLoad = (points_.size()*sizeof(label));

            format().openDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY)
                .closeTag();

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
            const uint64_t payLoad = (points_.size()*sizeof(label));

            format().openDataArray<label>(vtk::dataArrayAttr::OFFSETS)
                .closeTag();

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
            const uint64_t payLoad = (points_.size()*sizeof(label));

            format().openDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY)
                .closeTag();

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
            const uint64_t payLoad = (1*sizeof(label));

            format().openDataArray<label>(vtk::dataArrayAttr::OFFSETS)
                .closeTag();

            format().writeSize(payLoad);
            format().write(points_.size());
            format().flush();

            format().endDataArray();
        }

        format().endTag(vtk::fileTag::LINES);
        format().endTag(vtk::fileTag::PIECE);
    }

    // Standard corners in local axis
    const axesRotation cornerTransform(axis);

    FixedList<point, 4> corners;
    forAll(standardCorners, corni)
    {
        corners[corni] = cornerTransform.transform(standardCorners[corni]);
    }


    //
    // Planes to visualize location/rotation
    //
    {
        const label nPoints = 4*points_.size();  // 4 points per quad
        const label nPolys  = points_.size();

        format()
            .openTag(vtk::fileTag::PIECE)
            .xmlAttr(vtk::fileAttr::NUMBER_OF_POINTS, nPoints)
            .xmlAttr(vtk::fileAttr::NUMBER_OF_POLYS,  nPolys)
            .closeTag();

        // 'points'
        {
            const uint64_t payLoad = (nPoints*3*sizeof(float));

            format()
                .tag(vtk::fileTag::POINTS)
                .openDataArray<float, 3>(vtk::dataArrayAttr::POINTS)
                .closeTag();

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
            const uint64_t payLoad = (4*nPolys*sizeof(label));

            format().openDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY)
                .closeTag();

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
            const uint64_t payLoad = (nPolys*sizeof(label));

            format().openDataArray<label>(vtk::dataArrayAttr::OFFSETS)
                .closeTag();

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
        format().tag(vtk::fileTag::CELL_DATA);

        // zone Id
        {
            const uint64_t payLoad = (points_.size()*sizeof(label));

            format().openDataArray<label>("zoneId")
                .closeTag();

            format().writeSize(payLoad);
            for (label i=0; i < nPolys; ++i)
            {
                format().write(i);
            }
            format().flush();

            format().endDataArray();
        }

        format().endTag(vtk::fileTag::CELL_DATA);
#endif

        format().endTag(vtk::fileTag::PIECE);
    }

    // Finally
    // could add a 'ghost' level above to visualize extrapolated values
    // draw as two triangles to distingush from real levels ...

    format().endTag(vtk::fileTag::POLY_DATA);
    format().endVTKFile();
}


// ************************************************************************* //
