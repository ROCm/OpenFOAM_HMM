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
    Foam::faFieldReconstructor

Description
    Finite Area area and edge fields reconstruction.

Author
    Zeljko Tukovic, FSB Zagreb
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    faFieldReconstructor.C
    faFieldReconstructorFields.C

\*---------------------------------------------------------------------------*/

#ifndef faFieldReconstructor_H
#define faFieldReconstructor_H

#include "PtrList.H"
#include "faMesh.H"
#include "IOobjectList.H"
#include "faPatchFieldMapper.H"
#include "labelIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class faFieldReconstructor Declaration
\*---------------------------------------------------------------------------*/

class faFieldReconstructor
{
    // Private data

        //- Reconstructed mesh reference
        faMesh& mesh_;

        //- List of processor meshes
        const PtrList<faMesh>& procMeshes_;

        //- List of processor edge addressing lists
        const PtrList<labelIOList>& edgeProcAddressing_;

        //- List of processor face addressing lists
        const PtrList<labelIOList>& faceProcAddressing_;

        //- List of processor boundary addressing lists
        const PtrList<labelIOList>& boundaryProcAddressing_;


    // Private Member Functions

        //- No copy construct
        faFieldReconstructor(const faFieldReconstructor&) = delete;

        //- No copy assignment
        void operator=(const faFieldReconstructor&) = delete;


public:

        class faPatchFieldReconstructor
        :
            public faPatchFieldMapper
        {
            label size_;
            label sizeBeforeMapping_;

        public:

            // Constructors

                //- Construct given size
                faPatchFieldReconstructor
                (
                    const label size,
                    const label sizeBeforeMapping
                )
                :
                    size_(size),
                    sizeBeforeMapping_(sizeBeforeMapping)
                {}


            // Member functions

                virtual label size() const
                {
                    return size_;
                }

                virtual label sizeBeforeMapping() const
                {
                    return sizeBeforeMapping_;
                }

                virtual bool direct() const
                {
                    return true;
                }

                virtual bool hasUnmapped() const
                {
                    return false;
                }

                virtual const labelUList& directAddressing() const
                {
                    return labelUList::null();
                }
        };


    // Constructors

        //- Construct from components
        faFieldReconstructor
        (
            faMesh& mesh,
            const PtrList<faMesh>& procMeshes,
            const PtrList<labelIOList>& edgeProcAddressing,
            const PtrList<labelIOList>& faceProcAddressing,
            const PtrList<labelIOList>& boundaryProcAddressing
        );


    // Member Functions

        //- Reconstruct area field
        template<class Type>
        tmp<GeometricField<Type, faPatchField, areaMesh>>
        reconstructFaAreaField
        (
            const IOobject& fieldIoObject
        );

        //- Reconstruct edge field
        template<class Type>
        tmp<GeometricField<Type, faePatchField, edgeMesh>>
        reconstructFaEdgeField
        (
            const IOobject& fieldIoObject
        );

        //- Reconstruct and write all area fields
        template<class Type>
        void reconstructFaAreaFields
        (
            const IOobjectList& objects
        );

        //- Reconstruct and write all area fields
        template<class Type>
        void reconstructFaEdgeFields
        (
            const IOobjectList& objects
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "faFieldReconstructorFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
