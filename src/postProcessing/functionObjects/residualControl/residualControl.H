/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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
    Foam::residualControl

Description
    Function object that allows users to set target convergence criteria, and
    stop the run if the conditions are satisfied.

SourceFiles
    residualControl.C
    IOresidualControl.H

\*---------------------------------------------------------------------------*/

#ifndef residualControl_H
#define residualControl_H

#include "pointFieldFwd.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                       Class residualControl Declaration
\*---------------------------------------------------------------------------*/

class residualControl
{
protected:

    // Private data

        //- Name of this object
        word name_;

        //- Reference to object registry
        const objectRegistry& obr_;

        //- On/off switch - on if obr_ is an fvMesh object
        bool active_;

        //- List of variable name vs max residual
        List<Tuple2<word, scalar> > maxResiduals_;

        //- Flag to indicate whether convergence criteria have been met
        bool criteriaSatisfied_;


    // Protected Member Functions

        //- Perform residual control checks
        bool checkCriteria(const bool output) const;

        //- Disallow default bitwise copy construct
        residualControl(const residualControl&);

        //- Disallow default bitwise assignment
        void operator=(const residualControl&);


public:

    //- Runtime type information
    TypeName("residualControl");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        residualControl
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            const bool loadFromFilesUnused = false
        );


    //- Destructor
    virtual ~residualControl();


    // Member Functions

        //- Return name of the system call set
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the system calls
        virtual void read(const dictionary&);

        //- Execute the "executeCalls" at each time-step
        virtual void execute();

        //- Execute the "endCalls" at the final time-loop
        virtual void end();

        //- Write, execute the "writeCalls"
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const pointField&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
