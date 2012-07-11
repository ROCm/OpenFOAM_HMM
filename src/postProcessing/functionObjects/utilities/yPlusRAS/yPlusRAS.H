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

\addtogroup utilitiesFunctionObjects
@{

Class
    Foam::yPlusRAS

Description
    Evaluates and outputs turbulence y+ for RAS models.  Values written to
    time folders as field 'yPlusRAS'

SourceFiles
    yPlusRAS.C
    IOyPlusRAS.H
@}

\*---------------------------------------------------------------------------*/

#ifndef yPlusRAS_H
#define yPlusRAS_H

#include "volFieldsFwd.H"
#include "pointFieldFwd.H"
#include "Switch.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;
class fvMesh;

/*---------------------------------------------------------------------------*\
                          Class yPlusRAS Declaration
\*---------------------------------------------------------------------------*/

class yPlusRAS
{
    // Private data

        //- Name of this set of yPlusRAS objects
        word name_;

        const objectRegistry& obr_;

        //- on/off switch
        bool active_;

        //- Switch to send output to Info as well as to file
        Switch log_;

        //- Name of mass/volume flux field (optional, default = phi)
        word phiName_;

        //- Output file pointer
        autoPtr<OFstream> outputFilePtr_;


    // Private Member Functions

        //- Make the output file
        virtual void makeFile();

        //- Calculate incompressible form of y+
        void calcIncompressibleYPlus(const fvMesh& mesh, volScalarField& yPlus);

        //- Calculate compressible form of y+
        void calcCompressibleYPlus(const fvMesh& mesh, volScalarField& yPlus);

        //- Disallow default bitwise copy construct
        yPlusRAS(const yPlusRAS&);

        //- Disallow default bitwise assignment
        void operator=(const yPlusRAS&);


public:

    //- Runtime type information
    TypeName("yPlusRAS");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        yPlusRAS
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~yPlusRAS();


    // Member Functions

        //- Return name of the set of yPlusRAS
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the yPlusRAS data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Calculate the yPlusRAS and write
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
