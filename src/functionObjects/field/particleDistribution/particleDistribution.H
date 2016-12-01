/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::functionObjects::particleDistribution

Group
    grpFieldFunctionObjects

Description
    Generates a particle distrbution for lagrangian data.

Usage
    \verbatim
    particleDistribution1
    {
        type        particleDistribution;
        libs        ("libfieldFunctionObjects.so");
        ...
        cloud       "myCloud";
        field       "d"
        distributionBinWidth 0.1;
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property    | Description              |  Required  | Default value
        type        | Type name: particleDistribution | yes |
        cloud       | Name of cloud to process | Yes        |
        field       | Name of cloud field to process | Yes  |
        distributionBinWidth | Width of distribution bins | yes |
    \endtable

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile
    Foam::functionObjects::timeControl

SourceFiles
    particleDistribution.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_particleDistribution_H
#define functionObjects_particleDistribution_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "scalarField.H"
#include "cachedRandom.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                    Class particleDistribution Declaration
\*---------------------------------------------------------------------------*/

class particleDistribution
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected data

        //- Cloud name
        word cloudName_;

        //- Field name
        wordList fieldNames_;

        //- Tag field name - used to filter the particles into groups
        word tagFieldName_;

        //- Distribution bin width
        scalar distributionBinWidth_;

        //- Random number generator - used by distribution models
        cachedRandom rndGen_;


    // Private Member Functions

        //- Generate the distribution
        void generateDistribution
        (
            const word& fieldName,
            const scalarField& field,
            const label tag = -1
        );

        // Helper function to retrieve the particle data
        template<class Type>
        bool processField
        (
            const objectRegistry& obr,
            const word& fieldName,
            const List<DynamicList<label>>& addr
        );

        //- Disallow default bitwise copy construct
        particleDistribution(const particleDistribution&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const particleDistribution&) = delete;


public:

    //- Runtime type information
    TypeName("particleDistribution");


    // Constructors

        //- Construct from Time and dictionary
        particleDistribution
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~particleDistribution();


    // Member Functions

        //- Read the particleDistribution data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the particleDistribution
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "particleDistributionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
