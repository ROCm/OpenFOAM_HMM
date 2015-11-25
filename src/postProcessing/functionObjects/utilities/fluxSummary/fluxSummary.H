/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

Class
    Foam::fluxSummary

Group
    grpUtilityFunctionObjects

Description
    This function object calculates the flux across selections of faces.

    Output comprises, per set of faces, the:
    - positive
    - negative
    - net
    - absolute
    fluxes

    Example of function object specification:
    \verbatim
    fluxSummary1
    {
        type        fluxSummary;
        functionObjectLibs ("libutilityFunctionObjects.so");
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

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: fluxSummary  | yes         |
        write        | write flux data to file |  no  | yes
        log          | write flux data to standard output | no | yes
        mode         | mode to generate faces to test | yes |
        scaleFactor  | optional factor to scale result | no | 1
    \endtable

    The mode is one of:
    - faceZone
    - faceZoneAndDirection
    - cellZoneAndDirection

    Output data is written to files of the form \<timeDir\>/<faceZoneName>.dat

SeeAlso
    Foam::functionObject
    Foam::OutputFilterFunctionObject

SourceFiles
    fluxSummary.C
    IOfluxSummary.H

\*---------------------------------------------------------------------------*/

#ifndef fluxSummary_H
#define fluxSummary_H

#include "functionObjectFile.H"
#include "Switch.H"
#include "vector.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                    Class fluxSummary Declaration
\*---------------------------------------------------------------------------*/

class fluxSummary
:
    public functionObjectFile
{
public:

    // Public enumerations

        //- Face mode type
        enum modeType
        {
            mdFaceZone,
            mdFaceZoneAndDirection,
            mdCellZoneAndDirection
        };

        //- Mode type names
        static const NamedEnum<modeType, 3> modeTypeNames_;


protected:

    // Protected data

        //- Name of function object
        //  Also used as the name of the output directory
        word name_;

        //- Reference to the database
        const objectRegistry& obr_;

        //- on/off switch
        bool active_;

        //- Switch to send output to Info as well
        Switch log_;

        //- Mode for face determination
        modeType mode_;

        //- Scale factor
        scalar scaleFactor_;

        //- Name of flux field, default = phi
        word phiName_;


        // Per-faceZone information

            //- Region names
            List<word> faceZoneName_;

            //- Reference direction
            List<vector> refDir_;

            //- Face IDs
            List<List<label> > faceID_;

            //- Face patch IDs
            List<List<label> > facePatchID_;

            //- Face flip map signs
            List<List<scalar> > faceSign_;

            //- Sum of face areas
            List<scalar> faceArea_;

            //- Output file per face zone
            PtrList<OFstream> filePtrs_;


        //- Tolerance applied when matching face normals
        scalar tolerance_;


    // Private Member Functions

        //- Initialise face set from face zone
        void initialiseFaceZone
        (
            const word& faceZoneName,
            DynamicList<word>& faceZoneNames,
            DynamicList<List<label> >& faceID,
            DynamicList<List<label> >& facePatchID,
            DynamicList<List<scalar> >& faceSign
        ) const;

        //- Initialise face set from face zone and direction
        void initialiseFaceZoneAndDirection
        (
            const word& faceZoneName,
            const vector& refDir,
            DynamicList<vector>& dir,
            DynamicList<word>& faceZoneNames,
            DynamicList<List<label> >& faceID,
            DynamicList<List<label> >& facePatchID,
            DynamicList<List<scalar> >& faceSign
        ) const;

        //- Initialise face set from cell zone and direction
        void initialiseCellZoneAndDirection
        (
            const word& cellZoneName,
            const vector& refDir,
            DynamicList<vector>& dir,
            DynamicList<word>& faceZoneNames,
            DynamicList<List<label> >& faceID,
            DynamicList<List<label> >& facePatchID,
            DynamicList<List<scalar> >& faceSign
        ) const;

        //- Initialise the total area per derived faceZone
        void initialiseFaceArea();

        //- Output file header information
        virtual void writeFileHeader
        (
            const word& fzName,
            const scalar area,
            const vector& refDir,
            Ostream& os
        ) const;

        //- Disallow default bitwise copy construct
        fluxSummary(const fluxSummary&);

        //- Disallow default bitwise assignment
        void operator=(const fluxSummary&);


public:

    //- Runtime type information
    TypeName("fluxSummary");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        fluxSummary
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~fluxSummary();


    // Member Functions

        //- Return name of the set of field min/max
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the field min/max data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Write the fluxSummary
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
