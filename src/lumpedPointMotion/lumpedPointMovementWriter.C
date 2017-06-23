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

#include "lumpedPointMovement.H"
#include "polyMesh.H"
#include "OFstream.H"
#include "uindirectPrimitivePatch.H"

#include "foamVtkOutput.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::lumpedPointMovement::writeStateVTP(const fileName& file) const
{
    state().writeVTP(file, axis());
}


void Foam::lumpedPointMovement::writeForcesAndMomentsVTP
(
    const fileName& file,
    const UList<vector>& forces,
    const UList<vector>& moments
) const
{
    OFstream fos(file);
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
    const label nPoints = state().points().size();

    {
        format()
            .openTag(vtk::fileTag::PIECE)
            .xmlAttr(vtk::fileAttr::NUMBER_OF_POINTS, nPoints)
            .xmlAttr(vtk::fileAttr::NUMBER_OF_VERTS,  nPoints)
            .closeTag();

        // 'points'
        {
            const uint64_t payLoad = (nPoints*3* sizeof(float));

            format()
                .tag(vtk::fileTag::POINTS)
                .openDataArray<float, 3>(vtk::dataArrayAttr::POINTS)
                .closeTag();

            format().writeSize(payLoad);
            vtk::writeList(format(), state().points());
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
            const uint64_t payLoad = (nPoints*sizeof(label));

            format().openDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY)
                .closeTag();

            format().writeSize(payLoad);
            for (label i=0; i<nPoints; ++i)
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
            const uint64_t payLoad = (nPoints*sizeof(label));

            format().openDataArray<label>(vtk::dataArrayAttr::OFFSETS)
                .closeTag();

            format().writeSize(payLoad);
            for (label i=0; i<nPoints; ++i)
            {
                format().write(i+1);
            }
            format().flush();

            format().endDataArray();
        }

        format().endTag(vtk::fileTag::VERTS);
        // </Verts>
    }

    format().tag(vtk::fileTag::POINT_DATA);

    // forces
    if (forces.size() == nPoints)
    {
        const uint64_t payLoad = (nPoints * 3 * sizeof(float));

        format().openDataArray<float, 3>("forces")
            .closeTag();

        format().writeSize(payLoad);
        vtk::writeList(format(), forces);
        format().flush();

        format().endDataArray();
    }

    // moments
    if (moments.size() == nPoints)
    {
        const uint64_t payLoad = (nPoints * 3 * sizeof(float));

        format().openDataArray<float, 3>("moments")
            .closeTag();

        format().writeSize(payLoad);
        vtk::writeList(format(), moments);
        format().flush();

        format().endDataArray();
    }

    format().endTag(vtk::fileTag::POINT_DATA);

    format().endTag(vtk::fileTag::PIECE);
    format().endTag(vtk::fileTag::POLY_DATA);
    format().endVTKFile();
}


void Foam::lumpedPointMovement::writeZonesVTP
(
    const fileName& file,
    const polyMesh& mesh,
    const pointField& points0
) const
{
    OFstream fos(file);
    std::ostream& os = fos.stdStream();

    autoPtr<vtk::formatter> format = vtk::newFormatter
    (
        os,
        vtk::formatType::INLINE_ASCII
    );

    format().xmlHeader()
        .beginVTKFile(vtk::fileTag::POLY_DATA, "0.1");

    forAll(faceZones_, zoneI)
    {
        uindirectPrimitivePatch pp
        (
            UIndirectList<face>(mesh.faces(), faceZones_[zoneI]),
            points0
        );

        format()
            .openTag(vtk::fileTag::PIECE)
            .xmlAttr(vtk::fileAttr::NUMBER_OF_POINTS, pp.nPoints())
            .xmlAttr(vtk::fileAttr::NUMBER_OF_POLYS,  pp.size())
            .closeTag();

        // 'points'
        {
            const uint64_t payLoad = (pp.nPoints()*3* sizeof(float));

            format()
                .tag(vtk::fileTag::POINTS)
                .openDataArray<float, 3>(vtk::dataArrayAttr::POINTS)
                .closeTag();

            format().writeSize(payLoad);
            vtk::writeList(format(), pp.localPoints());
            format().flush();

            format()
                .endDataArray()
                .endTag(vtk::fileTag::POINTS);
        }

        // <Polys>
        format().tag(vtk::fileTag::POLYS);

        //
        // 'connectivity'
        //
        {
            uint64_t payLoad = 0;
            forAll(pp, facei)
            {
                const face& f = pp[facei];
                payLoad += f.size();
            }

            format().openDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY)
                .closeTag();

            format().writeSize(payLoad * sizeof(label));

            for (const face& f : pp.localFaces())
            {
                vtk::writeList(format(), f);
            }
            format().flush();

            format().endDataArray();
        }

        //
        // 'offsets'  (connectivity offsets)
        //
        {
            const uint64_t payLoad = (pp.size() * sizeof(label));

            format().openDataArray<label>(vtk::dataArrayAttr::OFFSETS)
                .closeTag();

            format().writeSize(payLoad);

            label off = 0;
            forAll(pp, facei)
            {
                const face& f = pp[facei];
                off += f.size();

                format().write(off);
            }
            format().flush();

            format().endDataArray();
        }

        format().endTag(vtk::fileTag::POLYS);


        format().tag(vtk::fileTag::CELL_DATA);

        // zone Id
        {
            const uint64_t payLoad = (pp.size() * sizeof(label));

            format().openDataArray<label>("zoneId")
                .closeTag();

            format().writeSize(payLoad);

            forAll(pp, facei)
            {
                format().write(zoneI);
            }
            format().flush();

            format().endDataArray();
        }

        format().endTag(vtk::fileTag::CELL_DATA);

        format().endTag(vtk::fileTag::PIECE);
    }

    format().endTag(vtk::fileTag::POLY_DATA);
    format().endVTKFile();
}


void Foam::lumpedPointMovement::writeVTP
(
    const fileName& file,
    const polyMesh& mesh,
    const labelUList& patchIds,
    const pointField& points0
) const
{
    writeVTP(file, state(), mesh, patchIds, points0);
}


void Foam::lumpedPointMovement::writeVTP
(
    const fileName& file,
    const lumpedPointState& state,
    const polyMesh& mesh,
    const labelUList& patchIds,
    const pointField& points0
) const
{
    const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();

    OFstream fos(file);
    std::ostream& os = fos.stdStream();

    autoPtr<vtk::formatter> format = vtk::newFormatter
    (
        os,
        vtk::formatType::INLINE_ASCII
    );

    format().xmlHeader()
        .beginVTKFile(vtk::fileTag::POLY_DATA, "0.1");

    for (const label patchId : patchIds)
    {
        const polyPatch& pp = boundaryMesh[patchId];

        format()
            .openTag(vtk::fileTag::PIECE)
            .xmlAttr(vtk::fileAttr::NUMBER_OF_POINTS, pp.nPoints())
            .xmlAttr(vtk::fileAttr::NUMBER_OF_POLYS,  pp.size())
            .closeTag();

        // 'points'
        {
            const uint64_t payLoad = (pp.nPoints()*3* sizeof(float));

            format()
                .tag(vtk::fileTag::POINTS)
                .openDataArray<float, 3>(vtk::dataArrayAttr::POINTS)
                .closeTag();

            // Could be more efficient, but not often needed
            tmp<pointField> tpts = displacePoints
            (
                state,
                points0,
                pp.meshPoints()
            ) + pointField(points0, pp.meshPoints());

            const pointField& pts = tpts();

            format().writeSize(payLoad);
            vtk::writeList(format(), pts);
            format().flush();

            format()
                .endDataArray()
                .endTag(vtk::fileTag::POINTS);
        }

        // <Polys>
        format().tag(vtk::fileTag::POLYS);


        //
        // 'connectivity'
        //
        {
            uint64_t payLoad = 0;
            for (const face& f : pp)
            {
                payLoad += f.size();
            }

            format().openDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY)
                .closeTag();

            format().writeSize(payLoad * sizeof(label));

            for (const face& f : pp.localFaces())
            {
                vtk::writeList(format(), f);
            }
            format().flush();

            format().endDataArray();
        }

        //
        // 'offsets'  (connectivity offsets)
        //
        {
            const uint64_t payLoad = (pp.size() * sizeof(label));

            format().openDataArray<label>(vtk::dataArrayAttr::OFFSETS)
                .closeTag();

            format().writeSize(payLoad);

            label off = 0;
            forAll(pp, facei)
            {
                const face& f = pp[facei];
                off += f.size();

                format().write(off);
            }
            format().flush();

            format().endDataArray();
        }

        format().endTag(vtk::fileTag::POLYS);

        format().endTag(vtk::fileTag::PIECE);
    }

    format().endTag(vtk::fileTag::POLY_DATA);
    format().endVTKFile();
}


// ************************************************************************* //
