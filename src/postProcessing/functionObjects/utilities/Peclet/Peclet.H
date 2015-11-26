/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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
    Foam::Peclet

Group
    grpUtilitiesFunctionObjects

Description
    This function object calculates and outputs the Peclet number as a
    surfaceScalarField.

    Example of function object specification to calculate the Peclet number:
    \verbatim
    Peclet1
    {
        type        Peclet;
        functionObjectLibs ("libutilityFunctionObjects.so");
        ...
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: Peclet       | yes         |
        phiName      | Name of flux field      | no          | phi
        resultName   | Name of Peclet field    | no          | <function name>
        log          | Log to standard output  | no          | yes
    \endtable

SourceFiles
    Peclet.C
    IOPeclet.H

\*---------------------------------------------------------------------------*/

#ifndef Peclet_H
#define Peclet_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "OFstream.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                          Class Peclet Declaration
\*---------------------------------------------------------------------------*/

class Peclet
{
    // Private data

        //- Name of this set of Peclet objects
        word name_;

        //- Reference to the database
        const objectRegistry& obr_;

        //- On/off switch
        bool active_;

        //- Name of flux field, default is "phi"
        word phiName_;

        //- Result name
        word resultName_;

        //- Switch to send output to Info as well as to file
        Switch log_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        Peclet(const Peclet&);

        //- Disallow default bitwise assignment
        void operator=(const Peclet&);


public:

    //- Runtime type information
    TypeName("Peclet");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        Peclet
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~Peclet();


    // Member Functions

        //- Return name of the set of Peclet
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the Peclet data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Calculate the Peclet and write
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
