/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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

Class
    Foam::MeshedSurfaceIOAllocator

Description
    A helper class for storing points, faces and zones with IO capabilities.

SourceFiles
    MeshedSurfaceIOAllocator.C

\*---------------------------------------------------------------------------*/

#ifndef MeshedSurfaceIOAllocator_H
#define MeshedSurfaceIOAllocator_H

#include "pointIOField.H"
#include "faceIOList.H"
#include "surfZoneIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class MeshedSurfaceIOAllocator Declaration
\*---------------------------------------------------------------------------*/

class MeshedSurfaceIOAllocator
{
    // Private data

        //- Points
        pointIOField points_;

        //- Faces
        faceCompactIOList faces_;

        //- Surface zones
        surfZoneIOList zones_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        MeshedSurfaceIOAllocator(const MeshedSurfaceIOAllocator&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const MeshedSurfaceIOAllocator&) = delete;


public:

    // Constructors

        //- Read construct from IOobjects
        MeshedSurfaceIOAllocator
        (
            const IOobject& ioPoints,
            const IOobject& ioFaces,
            const IOobject& ioZones
        );

        //- Construct from IOobjects, copying components
        MeshedSurfaceIOAllocator
        (
            const IOobject& ioPoints,
            const pointField& points,
            const IOobject& ioFaces,
            const faceList& faces,
            const IOobject& ioZones,
            const surfZoneList& zones
        );

        //- Construct from IOobjects, possibly transferring components
        MeshedSurfaceIOAllocator
        (
            const IOobject& ioPoints,
            const Xfer<pointField>& points,
            const IOobject& ioFaces,
            const Xfer<faceList>& faces,
            const IOobject& ioZones,
            const Xfer<surfZoneList>& zones
        );


    //- Destructor
    virtual ~MeshedSurfaceIOAllocator();


    // Member Functions

        // Database

            //- Set the instance for mesh files
            void setInstance(const fileName&);

            //- Adjust the write option for all components
            void setWriteOption(IOobject::writeOption);


        // Access

            //- Non-const access to the points
            pointIOField& storedIOPoints()
            {
                return points_;
            }

            //- Non-const access to the faces
            faceCompactIOList& storedIOFaces()
            {
                return faces_;
            }

            //- Non-const access to the zones
            surfZoneIOList& storedIOZones()
            {
                return zones_;
            }

            //- Const access to the points
            const pointIOField& storedIOPoints() const
            {
                return points_;
            }

            //- Const access to the faces
            const faceCompactIOList& storedIOFaces() const
            {
                return faces_;
            }

            //- Const access to the zones
            const surfZoneIOList& storedIOZones() const
            {
                return zones_;
            }


        //  Storage management

            //- Clear primitive data (points, faces and zones)
            void clear();

            //- Reset primitive data (faces and zones)
            //  Note, optimized to avoid overwriting data (with Xfer::null)
            void resetFaces
            (
                const Xfer<faceList>& faces,
                const Xfer<surfZoneList>& zones
            );

            //- Reset primitive data (points, faces and zones)
            //  Note, optimized to avoid overwriting data (with Xfer::null)
            void reset
            (
                const Xfer<pointField>& points,
                const Xfer<faceList>& faces,
                const Xfer<surfZoneList>& zones
            );

            //- Reset primitive data (points, faces and zones)
            //  Note, optimized to avoid overwriting data (with Xfer::null)
            void reset
            (
                const Xfer<List<point>>& points,
                const Xfer<faceList>& faces,
                const Xfer<surfZoneList>& zones
            );


        // Writing

            //- Write using given format, version and compression
            virtual bool writeObject
            (
                IOstream::streamFormat,
                IOstream::versionNumber,
                IOstream::compressionType,
                const bool valid
            ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
