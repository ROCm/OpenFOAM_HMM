/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::distributedUnallocatedDirectFvPatchFieldMapper

Description
    FieldMapper with direct mapping from remote quantities.

\*---------------------------------------------------------------------------*/

#ifndef distributedUnallocatedDirectFvPatchFieldMapper_H
#define distributedUnallocatedDirectFvPatchFieldMapper_H

#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class distributedUnallocatedDirectFvPatchFieldMapper Declaration
\*---------------------------------------------------------------------------*/

class distributedUnallocatedDirectFvPatchFieldMapper
:
    public fvPatchFieldMapper
{
    const labelUList& directAddressing_;

    const mapDistributeBase& distMap_;

    bool hasUnmapped_;

public:

    // Constructors

        //- Construct given addressing
        distributedUnallocatedDirectFvPatchFieldMapper
        (
            const labelUList& directAddressing,
            const mapDistributeBase& distMap
        )
        :
            directAddressing_(directAddressing),
            distMap_(distMap),
            hasUnmapped_(directAddressing_.size() && min(directAddressing_) < 0)
        {}


    //- Destructor
    virtual ~distributedUnallocatedDirectFvPatchFieldMapper() = default;


    // Member Functions

        virtual label size() const
        {
            return
            (
                notNull(directAddressing_)
              ? directAddressing_.size()
              : distMap_.constructSize()
            );
        }

        virtual bool direct() const
        {
            return true;
        }

        virtual bool distributed() const
        {
            return true;
        }

        virtual const mapDistributeBase& distributeMap() const
        {
            return distMap_;
        }

        virtual bool hasUnmapped() const
        {
            return hasUnmapped_;
        }

        virtual const labelUList& directAddressing() const
        {
            return directAddressing_;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
