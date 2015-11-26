/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd
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
    Foam::yPlus

Group
    grpUtilitiesFunctionObjects

Description
    This function object evaluates and outputs turbulence y+ for turbulence
    models. The field is stored on the mesh database so that it can be
    retrieved and used for other applications.

    Example of function object specification to calculate the y+ (LES):
    \verbatim
    yPlus1
    {
        type        yPlus;
        functionObjectLibs ("libutilityFunctionObjects.so");
        ...
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: yPlus        | yes         |
        phiName      | Name of flux field      | no          | phi
        resultName   | Name of y+ field        | no          | \<function name\>
        log          | Log to standard output  | no          | yes
    \endtable


SourceFiles
    yPlus.C
    IOyPlus.H

\*---------------------------------------------------------------------------*/

#ifndef yPlus_H
#define yPlus_H

#include "functionObjectFile.H"
#include "volFieldsFwd.H"
#include "Switch.H"
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
class fvMesh;

/*---------------------------------------------------------------------------*\
                            Class yPlus Declaration
\*---------------------------------------------------------------------------*/

class yPlus
:
    public functionObjectFile
{
    // Private data

        //- Name of this set of yPlus objects
        word name_;

        const objectRegistry& obr_;

        //- on/off switch
        bool active_;

        //- Name of mass/volume flux field (optional, default = phi)
        word phiName_;

        //- Result name
        word resultName_;

        //- Switch to send output to Info as well as to file
        Switch log_;


    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Calculate y+
        template<class TurbulenceModel>
        void calcYPlus
        (
            const TurbulenceModel& turbulenceModel,
            const fvMesh& mesh,
            volScalarField& yPlus
        );

        //- Disallow default bitwise copy construct
        yPlus(const yPlus&);

        //- Disallow default bitwise assignment
        void operator=(const yPlus&);


public:

    //- Runtime type information
    TypeName("yPlus");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        yPlus
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~yPlus();


    // Member Functions

        //- Return name of the set of yPlus
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the yPlus data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Calculate the yPlus and write
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

#ifdef NoRepository
    #include "yPlusTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
