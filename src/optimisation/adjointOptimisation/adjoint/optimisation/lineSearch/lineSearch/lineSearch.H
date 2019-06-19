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
    Foam::lineSearch

Description
    Abstract base class for optimisation methods

SourceFiles
    lineSearch.C

\*---------------------------------------------------------------------------*/

#ifndef lineSearch_H
#define lineSearch_H

#include "runTimeSelectionTables.H"
#include "IOdictionary.H"
#include "scalarField.H"
#include "stepUpdate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class lineSearch Declaration
\*---------------------------------------------------------------------------*/

class lineSearch
{
protected:

    // Protected data

        const dictionary dict_;
        IOdictionary lineSearchDict_;
        scalar directionalDeriv_;
        scalarField direction_;

        //- Old merit value from this opt cycle
        scalar oldMeritValue_;

        //- New merit value from this opt cycle
        scalar newMeritValue_;

        //- Merit directional deriv from the previous opt cycle
        scalar prevMeritDeriv_;
        scalar initialStep_;
        scalar minStep_;
        scalar step_;
        label iter_;
        label maxIters_;
        bool extrapolateInitialStep_;
        autoPtr<stepUpdate> stepUpdate_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        lineSearch(const lineSearch&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const lineSearch&) = delete;


public:

    //- Runtime type information
    TypeName("lineSearch");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            lineSearch,
            dictionary,
            (
                const dictionary& dict,
                const Time& time
            ),
            (dict, time)
        );


    // Constructors

        //- Construct from components
        lineSearch(const dictionary& dict, const Time& time);

    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<lineSearch> New
        (
            const dictionary& dict,
            const Time& time
        );


    //- Destructor
    virtual ~lineSearch() = default;


    // Member Functions

       //- Set objective derivative
       virtual void setDeriv(const scalar deriv);

       //- Set direction
       void setDirection(const scalarField& direction);

       //- Set new objective value
       void setNewMeritValue(const scalar value);

       //- Set old objective value
       void setOldMeritValue(const scalar value);

       //- Reset step to initial value
       virtual void reset();

       //- Get max number of iterations
       label maxIters() const;

       //- Get current step
       scalar step() const;

       //- Return the correction of the design variables
       virtual bool converged() = 0;

       //- Update the line search step based on the specific line search
       //- strategy, e.g. bisection, quadratic fit, etc.
       virtual void updateStep() = 0;

       //- Update the step using the supplied value
       virtual void updateStep(const scalar newStep);


    // Member operators

       //- Increment iteration number and store merit value corresponding to
       //- the previous optimisation cycle
       virtual lineSearch& operator++();

       //- Postfix increment. Necessary for compilation
       virtual lineSearch& operator++(int);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
