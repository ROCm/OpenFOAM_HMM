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
    Foam::faPatchMapper

Description
    Mapping class for a faPatchField.  Edge mapping is calculated based on
    faceCells comparison of old and new patch

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    faPatchMapper.C

\*---------------------------------------------------------------------------*/

#ifndef faPatchMapper_H
#define faPatchMapper_H

#include "faPatchFieldMapper.H"
#include "faceMapper.H"
#include "faPatch.H"
#include "primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class faPatch;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                        Class faPatchMapper Declaration
\*---------------------------------------------------------------------------*/

class faPatchMapper
:
    public faPatchFieldMapper
{
    // Private data

        //- Reference to patch
        const faPatch& patch_;

        //- Reference to mapPolyMesh
        const mapPolyMesh& mpm_;

        //- Size before mapping
        const label sizeBeforeMapping_;

        //- faceCells before mapping
        const labelList oldEdgeFaces_;


    // Demand-driven private data

        mutable bool hasUnmapped_;

        //- Direct addressing
        mutable labelList* directAddrPtr_;


    // Private Member Functions

        //- No copy construct
        faPatchMapper(const faPatchMapper&) = delete;

        //- No copy assignment
        void operator=(const faPatchMapper&) = delete;


        //- Calculate addressing for mapping with inserted cells
        void calcAddressing() const;

        //- Clear out local storage
        void clearOut();


public:

    //- Construct from mappers
    faPatchMapper
    (
        const faPatch& patch,
        const mapPolyMesh& mpm
    );


    //- Destructor
    virtual ~faPatchMapper();


    // Member Functions

        //- Return size
        virtual label size() const
        {
            return patch_.size();
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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
