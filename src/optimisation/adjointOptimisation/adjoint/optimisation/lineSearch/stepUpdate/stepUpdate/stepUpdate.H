/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2007-2019 PCOpt/NTUA
                            | Copyright (C) 2013-2019 FOSS GP
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
    Foam::stepUpdate

Description
    Abstract base class for optimisation methods

SourceFiles
    stepUpdate.C

\*---------------------------------------------------------------------------*/

#ifndef stepUpdate_H
#define stepUpdate_H

#include "runTimeSelectionTables.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class stepUpdate Declaration
\*---------------------------------------------------------------------------*/

class stepUpdate
{
protected:

    // Protected data

        const dictionary dict_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        stepUpdate(const stepUpdate&);

        //- Disallow default bitwise assignment
        void operator=(const stepUpdate&);


public:

    //- Runtime type information
    TypeName("stepUpdate");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            stepUpdate,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );


    // Constructors

        //- Construct from components
        stepUpdate(const dictionary& dict);


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<stepUpdate> New(const dictionary& dict);


    //- Destructor
    virtual ~stepUpdate() = default;


    // Member Functions

       //- Update the line search step
       virtual void updateStep(scalar& step) = 0;

       //- Set objective derivative
       virtual void setDeriv(const scalar deriv);

       //- Set new objective value
       virtual void setNewMeritValue(const scalar value);

       //- Set old objective value
       virtual void setOldMeritValue(const scalar value);

       //- Set old objective value
       virtual void setInitialStep(const scalar value);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
