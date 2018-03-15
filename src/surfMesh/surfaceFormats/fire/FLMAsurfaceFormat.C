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

#include "FLMAsurfaceFormat.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//! \cond fileScope
//- Output newline in ascii mode, no-op in binary mode
inline static void newline(Foam::OSstream& os)
{
    if (os.format() == Foam::IOstream::ASCII)
    {
        os  << '\n';
    }
}


template<class Face>
inline static int countFaces(const Face& f)
{
    int n = (f.size() - 2); // number triangles can be determined directly
    return n == 2 ? 1 : n;  // quads don't need triangulation
}

//! \endcond


template<class Face>
inline void Foam::fileFormats::FLMAsurfaceFormat<Face>::writeShell
(
    OSstream& os,
    const Face& f
)
{
    if (os.format() == IOstream::BINARY)
    {
        if (f.size() == 3 || f.size() == 4)
        {
            putFireLabel(os, f.size());
            for (const label verti : f)
            {
                putFireLabel(os, verti);
            }
        }
        else
        {
            // simple triangulation about f[0].
            // better triangulation should have been done before
            for (label fp1 = 1; fp1 < f.size() - 1; ++fp1)
            {
                const label fp2 = f.fcIndex(fp1);

                putFireLabel(os, 3);
                putFireLabel(os, f[0]);
                putFireLabel(os, f[fp1]);
                putFireLabel(os, f[fp2]);
            }
        }
    }
    else
    {
        // ASCII
        if (f.size() == 3 || f.size() == 4)
        {
            os  << ' ' << f.size();
            for (const label verti : f)
            {
                os  << ' ' << verti;
            }
            os  << nl;
        }
        else
        {
            for (label fp1 = 1; fp1 < f.size() - 1; ++fp1)
            {
                const label fp2 = f.fcIndex(fp1);
                os  << ' ' << 3 << ' '
                    << f[0] << ' ' << f[fp1] << ' ' << f[fp2]
                    << nl;
            }
        }
    }
}


template<class Face>
inline void Foam::fileFormats::FLMAsurfaceFormat<Face>::writeType
(
    OSstream& os,
    const Face& f
)
{
    if (os.format() == IOstream::BINARY)
    {
        if (f.size() == 4)
        {
            putFireLabel(os, fireQuad);
        }
        else
        {
            const label n = countFaces(f);
            for (label i=0; i < n; ++i)
            {
                putFireLabel(os, fireTri);
            }
        }
    }
    else
    {
        // ASCII
        if (f.size() == 4)
        {
            os  << ' ' << fireQuad;
        }
        else
        {
            const label n = countFaces(f);
            for (label i=0; i < n; ++i)
            {
                os  << ' ' << fireTri;
            }
        }
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Face>
void Foam::fileFormats::FLMAsurfaceFormat<Face>::write
(
    OSstream& os,
    const MeshedSurfaceProxy<Face>& surf
)
{
    if (!os.good())
    {
        FatalErrorInFunction
            << "bad output state "
            << exit(FatalError);
    }

    const UList<point>& pointLst = surf.points();
    const UList<Face>&  faceLst = surf.surfFaces();
    const UList<label>& faceMap = surf.faceMap();

    // for no zones, suppress the group name
    const List<surfZone>& zones =
    (
        surf.surfZones().empty()
      ? surfaceFormatsCore::oneZone(faceLst, word::null)
      : surf.surfZones()
    );

    const bool useFaceMap = (surf.useFaceMap() && zones.size() > 1);


    // determine the number of faces by counting the
    // tri/quads/triangulated) faces in each zone
    label nFaces = 0;
    List<label> zoneCount(zones.size());

    {
        label faceIndex = 0;
        forAll(zones, zoneI)
        {
            const surfZone& zone = zones[zoneI];

            label selCount = 0;
            if (useFaceMap)
            {
                forAll(zone, localFaceI)
                {
                    selCount += countFaces(faceLst[faceMap[faceIndex++]]);
                }
            }
            else
            {
                forAll(zone, localFaceI)
                {
                    selCount += countFaces(faceLst[faceIndex++]);
                }
            }

            zoneCount[zoneI] = selCount;
            nFaces += selCount;
        }
    }


    // Points
    // ~~~~~~

    // Set the precision of the points data to 10
    os.precision(10);

    Info<< nl << "points: " << pointLst.size() << endl;
    putFireLabel(os, pointLst.size());
    newline(os);

    for (const point& pt : pointLst)
    {
        // scaling is normally 1
        putFirePoint(os, pt);
    }
    newline(os); // readability

    // Faces indices
    {
        Info<< "faces:  " << nFaces << endl;
        putFireLabel(os, nFaces);
        newline(os);

        label faceIndex = 0;
        for (const surfZone& zone : zones)
        {
            const label nLocalFaces = zone.size();

            if (useFaceMap)
            {
                for (label i=0; i<nLocalFaces; ++i)
                {
                    writeShell(os, faceLst[faceMap[faceIndex++]]);
                }
            }
            else
            {
                for (label i=0; i<nLocalFaces; ++i)
                {
                    writeShell(os, faceLst[faceIndex++]);
                }
            }
        }
        newline(os);
        newline(os); // readability
    }


    // Face types
    {
        putFireLabel(os, nFaces);
        newline(os);

        label faceIndex = 0;
        for (const surfZone& zone : zones)
        {
            const label nLocalFaces = zone.size();

            if (useFaceMap)
            {
                for (label i=0; i<nLocalFaces; ++i)
                {
                    writeType(os, faceLst[faceMap[faceIndex++]]);
                }
            }
            else
            {
                for (label i=0; i<nLocalFaces; ++i)
                {
                    writeType(os, faceLst[faceIndex++]);
                }
            }
        }
        newline(os);
        newline(os); // readability
    }

    // Selections (cell)
    {
        putFireLabel(os, zones.size());
        newline(os);

        label faceIndex = 0;
        forAll(zones, zoneI)
        {
            const surfZone& zone = zones[zoneI];
            const label selCount = zoneCount[zoneI];

            putFireString(os, zone.name());
            putFireLabel(os, static_cast<int>(FIRECore::cellSelection));
            newline(os);

            putFireLabels(os, selCount, faceIndex);
            faceIndex += selCount;

            newline(os); // readability
        }
    }
}


template<class Face>
void Foam::fileFormats::FLMAsurfaceFormat<Face>::write
(
    bool compress,
    const fileName& filename,
    const MeshedSurfaceProxy<Face>& surf
)
{
    autoPtr<OFstream> osPtr
    (
        compress
      ? new OFstream
        (
            filename,
            IOstream::ASCII,
            IOstream::currentVersion,
            IOstream::COMPRESSED
        )
      : new OFstream(filename)
    );

    if (osPtr->good())
    {
        FLMAsurfaceFormat<Face>::write(osPtr(), surf);
        osPtr.clear();    // implicitly close the file

        if (compress)
        {
            // rename .flmaz.gz -> .flmaz
            // The '.gz' is automatically added by OFstream in compression mode
            Foam::mv(filename + ".gz", filename);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Cannot open file for writing " << filename
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
void Foam::fileFormats::FLMAsurfaceFormat<Face>::write
(
    const fileName& filename,
    const MeshedSurfaceProxy<Face>& surf,
    const dictionary&
)
{
    write(false, filename, surf);
}


template<class Face>
void Foam::fileFormats::FLMAZsurfaceFormat<Face>::write
(
    const fileName& filename,
    const MeshedSurfaceProxy<Face>& surf,
    const dictionary&
)
{
    FLMAsurfaceFormat<Face>::write(true, filename, surf);
}


// ************************************************************************* //
