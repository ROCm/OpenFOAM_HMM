/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

#include "VTPsurfaceFormat.H"
#include "OFstream.H"
#include "foamVtkOutput.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// File-scope constant.
//
// TODO: make this run-time selectable
// - No append mode supported
// - Legacy mode is dispatched via 'VTKsurfaceFormat' instead

static const Foam::foamVtkOutput::formatType fmtType =
    Foam::foamVtkOutput::formatType::INLINE_ASCII;
    // Foam::foamVtkOutput::formatType::INLINE_BASE64;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Face>
void Foam::fileFormats::VTPsurfaceFormat<Face>::writePolys
(
    foamVtkOutput::formatter& format,
    const UList<Face>& faces
)
{
    format.tag(vtkFileTag::POLYS);

    //
    // 'connectivity'
    //
    {
        uint64_t payLoad = 0;
        for (const auto& f : faces)
        {
            payLoad += f.size();
        }

        format.openDataArray<label>("connectivity")
            .closeTag();

        format.writeSize(payLoad * sizeof(label));

        for (const Face& f : faces)
        {
            foamVtkOutput::writeList(format, f);
        }

        format.flush();
        format.endDataArray();
    }


    //
    // 'offsets'  (connectivity offsets)
    //
    {
        const uint64_t payLoad(faces.size() * sizeof(label));

        format
            .openDataArray<label>("offsets")
            .closeTag();

        format.writeSize(payLoad);

        label off = 0;
        for (const auto& f : faces)
        {
            off += f.size();

            format.write(off);
        }

        format.flush();
        format.endDataArray();
    }

    format.endTag(vtkFileTag::POLYS);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::VTPsurfaceFormat<Face>::VTPsurfaceFormat()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
void Foam::fileFormats::VTPsurfaceFormat<Face>::write
(
    const fileName& filename,
    const MeshedSurfaceProxy<Face>& surf
)
{
    const pointField& pointLst = surf.points();
    const List<Face>&  faceLst = surf.surfFaces();
    const List<label>& faceMap = surf.faceMap();

    const List<surfZone>& zones =
    (
        surf.surfZones().empty()
      ? surfaceFormatsCore::oneZone(faceLst)
      : surf.surfZones()
    );

    const bool useFaceMap = (surf.useFaceMap() && zones.size() > 1);

    std::ofstream os(filename.c_str(), std::ios::binary);

    autoPtr<foamVtkOutput::formatter> format =
        foamVtkOutput::newFormatter(os, fmtType);

    writeHeader(format(), pointLst, faceLst.size());

    if (useFaceMap)
    {
        format().tag(vtkFileTag::POLYS);

        //
        // 'connectivity'
        //
        {
            uint64_t payLoad = 0;
            for (const auto& f : faceLst)
            {
                payLoad += f.size();
            }

            format().openDataArray<label>("connectivity")
                .closeTag();

            format().writeSize(payLoad * sizeof(label));

            label faceIndex = 0;
            for (const surfZone& zone : zones)
            {
                forAll(zone, i)
                {
                    const Face& f = faceLst[faceMap[faceIndex++]];

                    foamVtkOutput::writeList(format(), f);
                }
            }

            format().flush();
            format().endDataArray();
        }


        //
        // 'offsets'  (connectivity offsets)
        //
        {
            const uint64_t payLoad(faceLst.size() * sizeof(label));

            format()
                .openDataArray<label>("offsets")
                    .closeTag();

            format().writeSize(payLoad);

            label off = 0, faceIndex = 0;
            for (const surfZone& zone : zones)
            {
                forAll(zone, i)
                {
                    const Face& f = faceLst[faceMap[faceIndex++]];

                    off += f.size();

                    format().write(off);
                }
            }

            format().flush();
            format().endDataArray();
        }

        format().endTag(vtkFileTag::POLYS);
    }
    else
    {
        // Easy to write polys without a faceMap
        writePolys(format(), faceLst);
    }

    // Write regions (zones) as CellData
    if (zones.size() > 1)
    {
        writeCellData(format(), zones);
    }

    writeFooter(format());
}


template<class Face>
void Foam::fileFormats::VTPsurfaceFormat<Face>::write
(
    const fileName& filename,
    const UnsortedMeshedSurface<Face>& surf
)
{
    std::ofstream os(filename.c_str(), std::ios::binary);

    autoPtr<foamVtkOutput::formatter> format =
        foamVtkOutput::newFormatter(os, fmtType);

    const List<Face>& faceLst = surf.surfFaces();

    writeHeader(format(), surf.points(), faceLst.size());

    // Easy to write polys without a faceMap
    writePolys(format(), faceLst);

    // Write regions (zones) as CellData
    writeCellData(format(), surf.zoneIds());

    writeFooter(format());
}


// ************************************************************************* //
