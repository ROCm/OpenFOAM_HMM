/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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
    Foam::faEdgeMapper

Description
    FA edge mapper.  Currently, edge-based finite area data is not mapped,
    but only resized, since edge-based mapping data is not available

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    faEdgeMapper.C

\*---------------------------------------------------------------------------*/

#ifndef faEdgeMapper_H
#define faEdgeMapper_H

#include "morphFieldMapper.H"
#include "faMesh.H"
#include "faceMapper.H"
#include "HashSet.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class faEdgeMapper Declaration
\*---------------------------------------------------------------------------*/

class faEdgeMapper
:
    public morphFieldMapper
{
    // Private data

        //- Reference to mesh
        const faMesh& mesh_;

        //- Reference to mapPolyMesh
        //const mapPolyMesh& mpm_;

        //- Old mesh size
        label sizeBeforeMapping_;


    // Demand-driven private data

        mutable bool hasUnmapped_;

        //- Direct addressing
        mutable labelList* directAddrPtr_;


    // Private Member Functions

        //- No copy construct
        faEdgeMapper(const faEdgeMapper&) = delete;

        //- No copy assignment
        void operator=(const faEdgeMapper&) = delete;

        //- Calculate addressing
        void calcAddressing() const;

        //- Clear out local storage
        void clearOut();


public:

    //- Construct from components
    faEdgeMapper
    (
        const faMesh& mesh,
        const mapPolyMesh& mpm
    );


    //- Destructor
    virtual ~faEdgeMapper();


    // Member Functions

        //- Return size
        virtual label size() const
        {
            return mesh_.nInternalEdges();
        }

        //- Return size of field before mapping
        virtual label sizeBeforeMapping() const
        {
            return sizeBeforeMapping_;
        }

        //- Is the mapping direct
        virtual bool direct() const
        {
            return true;
        }

        virtual bool hasUnmapped() const
        {
            return hasUnmapped_;
        }

        //- Return direct addressing
        virtual const labelUList& directAddressing() const;

        //- Return interpolated addressing
        virtual const labelListList& addressing() const;

        //- Return interpolation weights
        virtual const scalarListList& weights() const;

        //- Are there any inserted faces
        virtual bool insertedObjects() const
        {
            return false;
        }

        //- Return list of inserted faces
        virtual const labelList& insertedObjectLabels() const
        {
            return labelList::null();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
