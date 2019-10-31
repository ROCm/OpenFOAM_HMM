/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::lagrangianReconstructor

Description
    Reconstructor for lagrangian positions and fields

SourceFiles
    lagrangianReconstructor.C
    lagrangianReconstructorFields.C

\*---------------------------------------------------------------------------*/

#ifndef lagrangianReconstructor_H
#define lagrangianReconstructor_H

#include "cloud.H"
#include "polyMesh.H"
#include "IOobjectList.H"
#include "CompactIOField.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class lagrangianReconstructor Declaration
\*---------------------------------------------------------------------------*/

class lagrangianReconstructor
{
    // Private data

        //- Mesh reference
        const fvMesh& mesh_;

        //- List of processor meshes
        const PtrList<fvMesh>& procMeshes_;

        //- List of processor face addressing lists
        const PtrList<labelIOList>& faceProcAddressing_;

        //- List of processor cell addressing lists
        const PtrList<labelIOList>& cellProcAddressing_;


    // Private Member Functions

        //- No copy construct
        lagrangianReconstructor(const lagrangianReconstructor&) = delete;

        //- No copy assignment
        void operator=(const lagrangianReconstructor&) = delete;


public:

    // Constructors

        //- Construct from components
        lagrangianReconstructor
        (
            const fvMesh& mesh,
            const PtrList<fvMesh>& procMeshes,
            const PtrList<labelIOList>& faceProcAddressing,
            const PtrList<labelIOList>& cellProcAddressing
        );


    // Member Functions

        //- Reconstruct positions for given cloud
        label reconstructPositions(const word& cloudName) const;

        //- Reconstruct a single field for given cloud
        template<class Type>
        tmp<IOField<Type>> reconstructField
        (
            const word& cloudName,
            const word& fieldName
        );

        //- Reconstruct a single field-field for given cloud
        template<class Type>
        tmp<CompactIOField<Field<Type>, Type>> reconstructFieldField
        (
            const word& cloudName,
            const word& fieldName
        );

        //- Reconstruct multiple fields for given cloud
        template<class Type>
        label reconstructFields
        (
            const word& cloudName,
            const IOobjectList& objects,
            const UList<word>& fieldNames
        );

        //- Reconstruct multiple fields for given cloud
        template<class Type>
        label reconstructFields
        (
            const word& cloudName,
            const IOobjectList& objects,
            const wordRes& selectedFields = wordRes()
        );

        //- Reconstruct multiple field-field for given cloud
        template<class Type>
        label reconstructFieldFields
        (
            const word& cloudName,
            const IOobjectList& objects,
            const wordRes& selectedFields = wordRes()
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "lagrangianReconstructorFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
