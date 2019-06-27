/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::functionObjects::fluxSummary

Group
    grpFieldFunctionObjects

Description
    This function object calculates the flux across selections of faces.

    Output comprises, per set of faces, the fluxes:
    - positive
    - negative
    - net
    - absolute

Usage
    \verbatim
    fluxSummary1
    {
        type        fluxSummary;
        libs        ("libfieldFunctionObjects.so");
        ...
        write       yes;
        log         yes;
        mode        cellZoneAndDirection;
        cellZoneAndDirection
        (
            (porosity (1 0 0))
        );
        scaleFactor 1.2;
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property    | Description                           | Required | Default
        type        | Type name: fluxSummary                | yes |
        write       | Write flux data to file               | no  | yes
        log         | Write flux data to standard output    | no  | yes
        mode        | Mode to generate faces to test        | yes |
        scaleFactor | Optional factor to scale result       | no  | 1
        tolerance   | Tolerance for the reference direction | no  | 0.8
        phi         | Surface flux field                    | no  | phi
    \endtable

    The mode is one of:
    - faceZone
    - faceZoneAndDirection
    - cellZoneAndDirection
    - surface
    - surfaceAndDirection

    Output data is written to files of the form \<timeDir\>/\<faceZoneName\>.dat

Note
    For surface and direction, can use phi = 'U' for determining the fluxes.

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile
    Foam::functionObjects::timeControl

SourceFiles
    fluxSummary.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fluxSummary_H
#define functionObjects_fluxSummary_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "vector.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dimensionSet;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class fluxSummary Declaration
\*---------------------------------------------------------------------------*/

class fluxSummary
:
    public fvMeshFunctionObject,
    public writeFile
{
public:

    // Public enumerations

        //- Face mode type
        enum modeType
        {
            mdFaceZone,             //!< Face zone
            mdFaceZoneAndDirection, //!< Face zone with prescribed direction
            mdCellZoneAndDirection, //!< Cell zone with prescribed direction
            mdSurface,              //!< A functionObject surface
            mdSurfaceAndDirection   //!< A surface with prescribed direction
        };

        //- Face mode names
        static const Enum<modeType> modeTypeNames_;


protected:

    // Protected Data

        //- Track if the surface needs an update
        bool needsUpdate_;

        //- Mode for face determination
        modeType mode_;

        //- Scale factor
        scalar scaleFactor_;

        //- Name of flux field, default = phi
        word phiName_;


        // Per-faceZone/surface information

            //- Region (zone/surface) names
            List<word> zoneNames_;

            //- Region (zone/surface) directions
            List<vector> zoneDirections_;

            //- Face IDs
            List<labelList> faceID_;

            //- Face patch IDs
            List<labelList> facePatchID_;

            //- Face flip map signs
            List<boolList> faceFlip_;

            //- Output file per face zone
            PtrList<OFstream> filePtrs_;


        //- Tolerance applied when matching face normals
        scalar tolerance_;


    // Private Member Functions

        //- Check if surface mode instead of zone mode
        bool isSurfaceMode() const;

        //- Check flowType (mass or volume).
        //  Return name on success, fatal error on failure.
        word checkFlowType
        (
            const dimensionSet& fieldDims,
            const word& fieldName
        ) const;

        //- Initialise for given surface name
        void initialiseSurface
        (
            const word& surfName,
            DynamicList<word>& names,
            DynamicList<vector>& dir,
            DynamicList<boolList>& faceFlip
        ) const;

        //- Initialise for given surface name and direction
        void initialiseSurfaceAndDirection
        (
            const word& surfName,
            const vector& refDir,
            DynamicList<word>& names,
            DynamicList<vector>& dir,
            DynamicList<boolList>& faceFlip
        ) const;

        //- Initialise face set from face zone
        void initialiseFaceZone
        (
            const word& faceZoneName,
            DynamicList<word>& names,
            DynamicList<vector>& dir,
            DynamicList<labelList>& faceID,
            DynamicList<labelList>& facePatchID,
            DynamicList<boolList>& faceFlip
        ) const;

        //- Initialise face set from face zone and direction
        void initialiseFaceZoneAndDirection
        (
            const word& faceZoneName,
            const vector& refDir,
            DynamicList<word>& names,
            DynamicList<vector>& dir,
            DynamicList<labelList>& faceID,
            DynamicList<labelList>& facePatchID,
            DynamicList<boolList>& faceFlip
        ) const;

        //- Initialise face set from cell zone and direction
        void initialiseCellZoneAndDirection
        (
            const word& cellZoneName,
            const vector& refDir,
            DynamicList<word>& names,
            DynamicList<vector>& dir,
            DynamicList<labelList>& faceID,
            DynamicList<labelList>& facePatchID,
            DynamicList<boolList>& faceFlip
        ) const;

        //- Calculate the total area for the surface or derived faceZone
        scalar totalArea(const label idx) const;

        //- Initialise - after read(), before write()
        bool update();


        //- Output file header information
        virtual void writeFileHeader
        (
            const word& zoneName,
            const scalar area,
            const vector& refDir,
            Ostream& os
        ) const;

        //- Specialized write for surfaces
        bool surfaceModeWrite();

        //- No copy construct
        fluxSummary(const fluxSummary&) = delete;

        //- No copy assignment
        void operator=(const fluxSummary&) = delete;


public:

    //- Runtime type information
    TypeName("fluxSummary");


    // Constructors

        //- Construct from Time and dictionary
        fluxSummary
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~fluxSummary() = default;


    // Member Functions

        //- Read the field fluxSummary data
        virtual bool read(const dictionary& dict);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the fluxSummary
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
