/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2011 OpenCFD Ltd.
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

    controllers

Description

    Stores all the information for the controllers models defined within
    the controllersDict, and selects & builds the models automatically.

\*---------------------------------------------------------------------------*/

#ifndef controllers_H
#define controllers_H

#include "List.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "polyMesh.H"
#include "timeData.H"
#include "stateController.H"
#include "fluxController.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class controllers Declaration
\*---------------------------------------------------------------------------*/

class controllers
{
    // Private data

        Time& time_;

        //- The entire dictionary (containing multiple subDictionaries)
        IOdictionary controllersDict_;

        //- state controllers
        PtrList<entry> stateControllersList_;
        List<word> sCNames_;
        List<label> sCIds_;
        List<fileName> sCFixedPathNames_;
        List< autoPtr<stateController> > stateControllers_;

        //- flux controllers
        PtrList<entry> fluxControllersList_;
        List<word> fCNames_;
        List<label> fCIds_;
        List<fileName> fCFixedPathNames_;
        List< autoPtr<fluxController> > fluxControllers_;


public:

    // Constructors

        //- Null Constructor
        controllers
        (
            const polyMesh& mesh
        );

        //- Constructor for with cloud
        controllers
        (
            const polyMesh& mesh,
            polyatomicCloud& cloud
        );

    //- Destructor
    ~controllers();

    // Member Functions

        //- Initial configuration call this function after the polyatomicCloud
        //  is completely initialised
        void initialConfig();

        //- this function is to be called at the beginning of the MD time-step.
        //  since we have placed a non-referenced time-data class in the
        //  state-controller class.
        void updateTimeInfo();

        //- control molecular state -- call this after the intermolecular force
        //  calulation
        void controlState();

        //-
        void controlVelocitiesI();

        //-
        void controlVelocitiesII();

        //-
        void controlPriorToForces();

        //- calculate properties -- call this at the end of the MD time-step.
        void calculateStateProps();

        //- output -- call this function at the end of the MD time-step
        void outputStateResults();

        // Access

            //-
            inline List< autoPtr<stateController> >& stateControllers();

            //-
            inline const List< autoPtr<stateController> >&
            stateControllers() const;

            //-
            inline List< autoPtr<fluxController> >& fluxControllers();

            //-
            inline const List< autoPtr<fluxController> >&
            fluxControllers() const;

            //-
            inline const List<word>& stateControllersNames() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "controllersI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
