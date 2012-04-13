/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::targetForceTrim

Description
    Target force trim coefficients

SourceFiles
    targetForceTrim.C

\*---------------------------------------------------------------------------*/

#ifndef targetForceTrim_H
#define targetForceTrim_H

#include "trimModel.H"
#include "tensor.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class targetForceTrim Declaration
\*---------------------------------------------------------------------------*/

class targetForceTrim
:
    public trimModel
{

protected:

    // Protected data

        //- Number of iterations between calls to 'correct'
        label calcFrequency_;

        //- Target force [N]
        vector target_;

        //- Pitch angles (collective, roll, pitch) [rad]
        vector theta_;

        //- Maximum number of iterations in trim routine
        label nIter_;

        //- Convergence tolerance
        scalar tol_;

        //- Under-relaxation coefficient
        scalar relax_;

        //- Perturbation angle used to determine jacobian
        scalar dTheta_;


    // Protected member functions

        //- Calculate the rotor forces
        vector calcForce
        (
            const vectorField& U,
            const scalarField& alphag,
            vectorField& force    
        ) const;


public:

    //- Run-time type information
    TypeName("targetForceTrim");

    //- Constructor
    targetForceTrim(const rotorDiskSource& rotor, const dictionary& dict);

    //- Destructor
    virtual ~targetForceTrim();


    // Member functions

        //- Read
        void read(const dictionary& dict);

        //- Return the geometric angle of attack [rad]
        virtual tmp<scalarField> thetag() const;

        //- Correct the model
        virtual void correct(const vectorField& U, vectorField& force);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
