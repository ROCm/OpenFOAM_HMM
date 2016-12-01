/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenCFD Ltd.
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
    Foam::functionObjects::extractEulerianParticleDistribution
Group
    grpFieldFunctionObjects

Description
    Generates particle size information from Eulerian calculations, e.g. VoF.

    Particle data is written to the directory:

    \verbatim
        $FOAM_CASE/postProcessing/<name>/particleDistribution
    \endverbatim

Usage
    extractEulerianParticleDistribution1
    {
        type            extractEulerianParticleDistribution;
        libs            ("libfieldFunctionObjects.so");
        ...
        faceZone        f0;
        nLocations      10;
        alphaName       alpha.water;
        UName           U;
        rhoName         rho;
        phiName         phi;

        writeRawData    yes;
    }
    \endverbatim

    where the entries comprise:
    \table
        Property     | Description             | Required | Default value
        type         | type name: extractEulerianParticleDistribution | yes |
        faceZone     | Name of faceZone used as collection surface | yes |
        nLocations   | Number of injection bins to generate  | yes |
        aplhaName    | Name of phase indicator field | yes   |
        rhoName      | Name of density field   |      yes |
        phiName      | Name of flux field      |      yes |
        distributionBinWidth | Binned distribution bin width| yes |
        writeRawData | Flag to write raw particle data | yes |
    \endtable

SourceFiles
    extractEulerianParticleDistribution.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_extractEulerianParticleDistribution_H
#define functionObjects_extractEulerianParticleDistribution_H

#include "extractEulerianParticles.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
/*---------------------------------------------------------------------------*\
     Class extractEulerianParticleDistributionFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class extractEulerianParticleDistribution
:
    public extractEulerianParticles
{

protected:

    // Protected data

        // Agglomeration

            //- Number of sample locations to generate
            label nInjectorLocations_;

            //- Agglomeration addressing from fine to coarse
            labelList fineToCoarseAddr_;

            //- Agglomeration addressing from coarse to fine
            labelListList coarseToFineAddr_;

            //- Coarse face positions
            vectorList coarsePosition_;

            //- Global coarse face addressing
            globalIndex globalCoarseFaces_;


        // Particle collection info

            //- Flag to reset the distribution on each write
            bool resetDistributionOnWrite_;

            //- Diameter distribution bin width
            scalar distributionBinWidth_;

            //- Random class needed by distribution models
            cachedRandom rndGen_;


    // Protected Member Functions

        //- Initialise the particle collection bins
        virtual void initialiseBins();

        //- Write agglomerated particle data to stream
        virtual void writeBinnedParticleData();

        //- Disallow default bitwise copy construct
        extractEulerianParticleDistribution
        (
            const extractEulerianParticleDistribution&
        );

        //- Disallow default bitwise assignment
        void operator=(const extractEulerianParticleDistribution&);


public:

    // Static data members

        //- Static data staticData
        TypeName("extractEulerianParticleDistribution");


    // Constructors

        //- Construct from components
        extractEulerianParticleDistribution
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool readFields = true
        );


    //- Destructor
    virtual ~extractEulerianParticleDistribution();


    // Member Functions

        //- Read the field min/max data
        virtual bool read(const dictionary&);

        //- Execute
        virtual bool execute();

        //- Write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "extractEulerianParticleDistributionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
