/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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
    stateController

Description

    Basic/abstract class of a state controller

SourceFiles
    stateControllerI.H
    stateController.C
    stateControllerIO.C

\*---------------------------------------------------------------------------*/

#ifndef stateController_H
#define stateController_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "Random.H"
#include "polyatomic.H"
#include "timeData.H"
#include "writeTimeData.H"
#include "selectIds.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class stateController Declaration
\*---------------------------------------------------------------------------*/

class stateController
{

protected:

    // Protected data

        //-
        const fvMesh& mesh_;

        //-
        polyatomicCloud& cloud_;

        //-
        Random rndGen_;

        //- subDictionary containing the properties
        dictionary controllerDict_;

        //-
        dictionary timeDict_;

        //-
        timeData time_;

        //-
        scalar timePeriod_;

        //-
        scalar initialTime_;

        //- name of control zone
        word regionName_;

        //-
        label regionId_;

        //-
        bool control_;

        //-
        bool readStateFromFile_;

        //- set all the properties below from model if required

            //-
            bool singleValueController_;

            //- target values
            scalar density_;
            vector velocity_;
            scalar temperature_;
            scalar pressure_;
            tensor strainRate_;
            vector tempGradient_;

            //- set this in model
            bool fieldController_;

            //- targeted fields
            scalarField densities_;
            vectorField velocities_;
            scalarField temperatures_;
            scalarField pressures_;

            // set these in model
            bool writeInTimeDir_;
            bool writeInCase_;


public:

    //- Runtime type information
    TypeName("stateController");

    //- Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            stateController,
            dictionary,
            (
                polyatomicCloud& cloud,
                const dictionary& dict
            ),
            (t, cloud, dict)
        );

    // Constructors

        //- Construct from components
        stateController
        (
            polyatomicCloud& cloud,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<stateController> New
        (
            polyatomicCloud& cloud,
            const dictionary& dict
        );


    // Destructor

        virtual ~stateController();


    // Member Functions


        void updateTime();

        //- create an initial configuration
        virtual void initialConfiguration() = 0;

        //- calculate any required properties
        virtual void calculateProperties() = 0;

        //- control molecules at different stages of the integration time-step
        virtual void controlMolsBeg() = 0;

        virtual void controlBeforeForces() = 0;

        virtual void controlMols() = 0;

        virtual void controlMolsEnd() = 0;


        //- output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        ) = 0;

        //- update properties from a modified dictionary
        virtual void updateProperties(const dictionary&) = 0;

        void updateStateControllerProperties(const dictionary&);

        // Access

            //- return the control zone cells
            const labelList& controlZone() const;

            //- return the control zone name
            const word& regionName() const;

            //- return the targeted fields
            scalar density() const;
            scalar& density();

            const vector& velocity() const;
            vector& velocity();

            scalar temperature() const;
            scalar& temperature();

            scalar pressure() const;
            scalar& pressure();

            const tensor& strainRate() const;
            tensor& strainRate();

            const vector& tempGradient() const;
            vector& tempGradient();

            //- return the targeted fields
            const scalarField& densityField() const;
            scalarField& densityField();

            const vectorField& velocityField() const;
            vectorField& velocityField();

            const scalarField& temperatureField() const;
            scalarField& temperatureField();

            const scalarField& pressureField() const;
            scalarField& pressureField();

            bool singleValueController() const;
            bool& singleValueController();

            bool fieldController() const;
            bool& fieldController();

            bool writeInTimeDir() const;
            bool writeInCase() const;

            scalar avReqDensity();
            vector avReqVelocity();
            scalar avReqTemperature();
            scalar avReqPressure();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
