/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::lduPrimitiveMesh

Description
    Simplest contrete lduMesh which stores the addressing needed by lduMatrix.

SourceFiles
    lduPrimitiveMesh.C

\*---------------------------------------------------------------------------*/

#ifndef lduPrimitiveMesh_H
#define lduPrimitiveMesh_H

#include "lduMesh.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//// Forward declaration of friend functions and operators
//
//class lduPrimitiveMesh;
//
//Ostream& operator<<(Ostream&, const lduPrimitiveMesh&);

/*---------------------------------------------------------------------------*\
                       Class lduPrimitiveMesh Declaration
\*---------------------------------------------------------------------------*/

class lduPrimitiveMesh
:
    public lduMesh,
    public lduAddressing
{
    // Private data

        //- Lower addressing
        labelList lowerAddr_;

        //- Upper addressing
        labelList upperAddr_;

        //- Patch to internal addressing
        labelListList patchAddr_;

        //- List of pointers for each patch
        //  with only those pointing to interfaces being set
        const lduInterfacePtrsList interfaces_;

        //- Patch field evaluation schedule
        const lduSchedule patchSchedule_;

        //- Communicator to use for any parallel communication
        const label comm_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        lduPrimitiveMesh(const lduPrimitiveMesh&);

        //- Disallow default bitwise assignment
        void operator=(const lduPrimitiveMesh&);


public:

    // Constructors

        //- Construct from components as copies
        lduPrimitiveMesh
        (
            const label nCells,
            const labelUList& l,
            const labelUList& u,
            const labelListList& pa,
            const lduInterfacePtrsList& interfaces,
            const lduSchedule& ps,
            const label comm
        );

        //- Construct from components and re-use storage as specified.
        lduPrimitiveMesh
        (
            const label nCells,
            labelList& l,
            labelList& u,
            labelListList& pa,
            lduInterfacePtrsList& interfaces,
            const lduSchedule& ps,
            const label comm,
            bool reUse
        );

        //- Construct from Istream
        lduPrimitiveMesh(Istream& is);


    //- Destructor
    virtual ~lduPrimitiveMesh()
    {}


    // Member Functions

        // Access

            //- Return ldu addressing
            virtual const lduAddressing& lduAddr() const
            {
                return *this;
            }

            //- Return a list of pointers for each patch
            //  with only those pointing to interfaces being set
            virtual lduInterfacePtrsList interfaces() const
            {
                return interfaces_;
            }

            //- Return communicator used for parallel communication
            virtual label comm() const
            {
                return comm_;
            }

            //- Return Lower addressing
            virtual const labelUList& lowerAddr() const
            {
                return lowerAddr_;
            }

            //- Return Upper addressing
            virtual const labelUList& upperAddr() const
            {
                return upperAddr_;
            }

            //- Return patch addressing
            virtual const labelUList& patchAddr(const label i) const
            {
                return patchAddr_[i];
            }

            //- Return patch evaluation schedule
            virtual const lduSchedule& patchSchedule() const
            {
                return patchSchedule_;
            }
//
//    // Ostream operator
//
//        friend Ostream& operator<<(Ostream&, const procLduMatrix&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
