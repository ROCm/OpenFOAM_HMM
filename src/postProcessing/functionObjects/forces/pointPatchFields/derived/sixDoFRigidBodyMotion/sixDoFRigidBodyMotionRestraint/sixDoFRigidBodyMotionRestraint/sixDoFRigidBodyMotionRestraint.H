/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Namespace
    Foam::sixDoFRigidBodyMotionRestraints

Description
    Namespace for six DoF motion restraints


Class
    Foam::sixDoFRigidBodyMotionRestraint

Description
    Base class for defining restraints for sixDoF motions

SourceFiles
    sixDoFRigidBodyMotionRestraint.C
    newDynamicFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef sixDoFRigidBodyMotionRestraint_H
#define sixDoFRigidBodyMotionRestraint_H

#include "Time.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "vector.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class sixDoFRigidBodyMotion;


/*---------------------------------------------------------------------------*\
                Class sixDoFRigidBodyMotionRestraint Declaration
\*---------------------------------------------------------------------------*/

class sixDoFRigidBodyMotionRestraint
{

protected:

    // Protected data

        dictionary sDoFRBMRCoeffs_;


public:

    //- Runtime type information
    TypeName("sixDoFRigidBodyMotionRestraint");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            sixDoFRigidBodyMotionRestraint,
            dictionary,
            (const dictionary& sDoFRBMRCoeffs),
            (sDoFRBMRCoeffs)
        );


    // Constructors

        //- Construct from the sDoFRBMRCoeffs dictionary and Time
        sixDoFRigidBodyMotionRestraint
        (
            const dictionary& sDoFRBMRCoeffs
        );

        //- Construct and return a clone
        virtual autoPtr<sixDoFRigidBodyMotionRestraint> clone() const = 0;


    // Selectors

        //- Select constructed from the sDoFRBMRCoeffs dictionary and Time
        static autoPtr<sixDoFRigidBodyMotionRestraint> New
        (
            const dictionary& sDoFRBMRCoeffs
        );


    // Destructor

        virtual ~sixDoFRigidBodyMotionRestraint();


    // Member Functions

        //- Calculate the restraint position, force and moment.
        //  Global reference frame vectors.
        virtual void restrain
        (
            const sixDoFRigidBodyMotion& motion,
            vector& restraintPosition,
            vector& restraintForce,
            vector& restraintMoment
        ) const = 0;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& sDoFRBMRCoeffs) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
