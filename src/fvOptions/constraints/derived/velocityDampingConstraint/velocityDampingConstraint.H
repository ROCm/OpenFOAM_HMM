/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2015 OpenCFD Ltd.
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
    Foam::fv::velocityDampingConstraint

Group
    grpFvOptionsConstraints

Description
    Constraint on velocity field to dampen velocity fluctuations.

    This constraint is primarily used to dampen velocity fluctuations in
    the start-up phase of simulations. When the local velocity magnitude
    exceeds the user-supplied maximum value a sink term is activated in
    the affected region to lower the velocity to the limiting value.

Usage
    Example of functionality specification in \c fvOptions:
    \verbatim
    velocityDamper
    {
        // Mandatory entries
        type            velocityDampingConstraint;
        UMax            200;

        // Optional entries
        selectionMode   all;
        UNames          (U);    // names of given velocity fields
        U               U;      // name of given velocity field if
                                //`UNames` is not present (default: U)

        // Optional fvOptions entries
        active          true;
    }
    \endverbatim

Note
    When active, this constraint manipulates the system of equations.
    Users should ensure that it is not active when the case is converged
    (steady-state) or during the period of interest (transient) to ensure
    that its presence does not pollute the results.

SourceFiles
    velocityDampingConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef velocityDampingConstraint_H
#define velocityDampingConstraint_H

#include "cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                  Class velocityDampingConstraint Declaration
\*---------------------------------------------------------------------------*/

class velocityDampingConstraint
:
    public cellSetOption
{

protected:

    // Protected Data

        //- Maximum velocity magnitude
        scalar UMax_;


    // Protected Member Functions

        void addDamping(fvMatrix<vector>& eqn);


private:

    // Private Member Functions

        //- No copy construct
        velocityDampingConstraint(const velocityDampingConstraint&) = delete;

        //- No copy assignment
        void operator=(const velocityDampingConstraint&) = delete;


public:

    //- Runtime type information
    TypeName("velocityDampingConstraint");


    // Constructors

        //- Construct from components
        velocityDampingConstraint
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~velocityDampingConstraint() = default;


    // Member Functions

        // Set values directly

            //- Constrain vector matrix
            virtual void constrain(fvMatrix<vector>& eqn, const label fieldi);


        // I-O

            //- Write data
            virtual void writeData(Ostream& os) const;

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
