/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenCFD Ltd.
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
    Foam::multiphaseSystem

Description
    Class which solves the volume fraction equations for multiple phases

SourceFiles
    multiphaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseSystem_H
#define multiphaseSystem_H

#include "phaseSystem.H"
#include "UPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class multiphaseSystem Declaration
\*---------------------------------------------------------------------------*/

class multiphaseSystem
:
    public phaseSystem
{
protected:

    // Protected data

        //- Unallocated phase list
        UPtrList<phaseModel> phases_;

        //- dmdt continuity error
        volScalarField dmdtContinuityError_;


public:

    //- Runtime type information
    TypeName("multiphaseSystem");

    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            multiphaseSystem,
            dictionary,
            (
                const fvMesh& mesh
            ),
            (mesh)
        );


    // Constructors

        //- Construct from fvMesh
        multiphaseSystem(const fvMesh&);


    //- Destructor
    virtual ~multiphaseSystem();


    // Selectors

        static autoPtr<multiphaseSystem> New(const fvMesh& mesh);


     // Member Functions

        //- Solve for the phase fractions
        virtual void solve();


        // Access


            //- Constant access the dmdt continuity error
            const volScalarField& dmdtContinuityError() const;

            //- Access the dmdt continuity error
            volScalarField& dmdtContinuityError();

            //- Return phases
            const UPtrList<phaseModel>& phases() const;

            //- Return phases
            UPtrList<phaseModel>& phases();

            //- Constant access phase model i
            virtual const phaseModel& phase(const label i) const;

            //- Access phase model i
            virtual phaseModel& phase(const label i);


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
