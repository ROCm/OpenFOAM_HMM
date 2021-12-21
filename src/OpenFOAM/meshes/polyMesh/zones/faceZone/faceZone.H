/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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
    Foam::faceZone

Description
    A subset of mesh faces organised as a primitive patch.

    For quick check whether a face belongs to the zone use the lookup
    mechanism in faceZoneMesh, where all the zoned faces are registered
    with their zone number.

SourceFiles
    faceZone.C
    faceZoneNew.C

\*---------------------------------------------------------------------------*/

#ifndef faceZone_H
#define faceZone_H

#include "zone.H"
#include "faceZoneMeshFwd.H"
#include "primitiveFacePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class mapPolyMesh;
class faceZone;
Ostream& operator<<(Ostream& os, const faceZone& zn);


/*---------------------------------------------------------------------------*\
                          Class faceZone Declaration
\*---------------------------------------------------------------------------*/

class faceZone
:
    public zone
{
    // Private Data

        //- Flip map for all faces in the zone.
        //  True if the face needs flipping for the correct orientation.
        boolList flipMap_;

        //- Reference to zone list
        const faceZoneMesh& zoneMesh_;

        //- Demand-driven: Primitive patch of correctly flipped faces
        mutable primitiveFacePatch* patchPtr_;

        //- Demand-driven: Master cell layer
        mutable labelList* masterCellsPtr_;

        //- Demand-driven: Slave cell layer
        mutable labelList* slaveCellsPtr_;

        //- Demand-driven: Global edge addressing
        mutable labelList* mePtr_;


    // Private Member Functions

        //- Set flip-map to uniform value
        void setFlipMap(const bool val);

        //- Build primitive patch
        void calcFaceZonePatch() const;

        //- Calculate master and slave face layer
        void calcCellLayers() const;

        //- Check addressing
        void checkAddressing() const;


        //- No copy construct
        faceZone(const faceZone&) = delete;

        //- No copy assignment
        void operator=(const faceZone&) = delete;


public:

    // Static Data Members

        //- The name associated with the zone-labels dictionary entry
        //- ("faceLabels")
        static const char * const labelsName;


    //- Runtime type information
    TypeName("faceZone");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            faceZone,
            dictionary,
            (
                const word& name,
                const dictionary& dict,
                const label index,
                const faceZoneMesh& zm
            ),
            (name, dict, index, zm)
        );


    // Constructors

        //- Construct an empty zone
        faceZone
        (
            const word& name,
            const label index,
            const faceZoneMesh& zm
        );

        //- Construct from components with uniform flip map value
        faceZone
        (
            const word& name,
            const labelUList& addr,
            const bool flipMapValue,
            const label index,
            const faceZoneMesh& zm
        );

        //- Construct from components with uniform flip map value,
        //- transferring addressing.
        faceZone
        (
            const word& name,
            labelList&& addr,
            const bool flipMapValue,
            const label index,
            const faceZoneMesh& zm
        );

        //- Construct from components
        faceZone
        (
            const word& name,
            const labelUList& addr,
            const boolUList& fm,
            const label index,
            const faceZoneMesh& zm
        );

        //- Construct from components, transferring addressing
        faceZone
        (
            const word& name,
            labelList&& addr,
            boolList&& fm,
            const label index,
            const faceZoneMesh& zm
        );

        //- Construct from dictionary
        faceZone
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const faceZoneMesh& zm
        );

        //- Construct with a new index and zone mesh information, the name
        //- of the original zone, resetting the face addressing
        //- and flip-map.
        faceZone
        (
            const faceZone& origZone,
            const labelUList& addr,
            const boolUList& fm,
            const label index,
            const faceZoneMesh& zm
        );

        //- Construct with a new index and zone mesh information, the name
        //- of the original zone, (move) resetting the face addressing
        //- and flip-map.
        faceZone
        (
            const faceZone& origZone,
            labelList&& addr,
            boolList&& fm,
            const label index,
            const faceZoneMesh& zm
        );

        //- Construct and return a clone, resetting the zone mesh
        virtual autoPtr<faceZone> clone(const faceZoneMesh& zm) const
        {
            return autoPtr<faceZone>::New(*this, *this, flipMap(), index(), zm);
        }

        //- Construct and return a clone,
        //- resetting the face list and zone mesh
        virtual autoPtr<faceZone> clone
        (
            const labelUList& addr,
            const boolUList& fm,
            const label index,
            const faceZoneMesh& zm
        ) const
        {
            return autoPtr<faceZone>::New(*this, addr, fm, index, zm);
        }


    // Selectors

        //- Return a pointer to a new face zone
        //- created on freestore from dictionary
        static autoPtr<faceZone> New
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const faceZoneMesh& zm
        );


    //- Destructor
    virtual ~faceZone();


    // Member Functions

        //- Return reference to the zone mesh
        const faceZoneMesh& zoneMesh() const noexcept
        {
            return zoneMesh_;
        }

        //- Return face flip map
        const boolList& flipMap() const noexcept
        {
            return flipMap_;
        }

        //- Helper function to re-direct to zone::localID(...)
        label whichFace(const label globalCellID) const;

        //- Return reference to primitive patch
        const primitiveFacePatch& operator()() const;


    // Addressing into mesh

        //- Return labels of master cells (cells next to the master face
        //- zone in the prescribed direction)
        const labelList& masterCells() const;

        //- Return labels of slave cells
        const labelList& slaveCells() const;

        //- Return global edge index for local edges
        const labelList& meshEdges() const;


        //- Clear addressing
        virtual void clearAddressing();

        //- Reset addressing - use uniform flip map value
        //  Clears demand-driven data.
        virtual void resetAddressing
        (
            const labelUList& addr,
            const bool flipMapValue
        );

        //- Reset addressing and flip map.
        //  Clears demand-driven data.
        virtual void resetAddressing
        (
            const labelUList& addr,
            const boolUList& flipMap
        );

        //- Move reset addressing - use uniform flip map value
        //  Clears demand-driven data.
        virtual void resetAddressing
        (
            labelList&& addr,
            const bool flipMapValue
        );


        //- Check zone definition. Return true if in error.
        virtual bool checkDefinition(const bool report = false) const;

        //- Check whether all procs have faces synchronised.
        //  \return True if any errors.
        virtual bool checkParallelSync(const bool report = false) const;

        //- Correct patch after moving points
        virtual void movePoints(const pointField& pts);

        //- Update for changes in topology
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Write
        virtual void write(Ostream& os) const;

        //- Write dictionary
        virtual void writeDict(Ostream& os) const;


    // I-O

        //- Ostream Operator
        friend Ostream& operator<<(Ostream& os, const faceZone& zn);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
