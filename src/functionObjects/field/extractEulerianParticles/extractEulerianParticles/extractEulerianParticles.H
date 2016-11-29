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
    Foam::functionObjects::extractEulerianParticles
Group
    grpFieldFunctionObjects

Description
    Generates particle size information from Eulerian calculations, e.g. VoF.

    Particle data is written to the directory:

    \verbatim
        $FOAM_CASE/postProcessing/<name>/rawParticlesDict
        $FOAM_CASE/postProcessing/<name>/particleDistribution
    \endverbatim

Usage
    extractEulerianParticles1
    {
        type            extractEulerianParticles;
        libs            ("libfieldFunctionObjects.so");
        ...
        faceZone        f0;
        nLocations      10;
        alphaName       alpha.water;
        UName           U;
        rhoName         rho;
        phiName         phi;

        createDistribution yes;
        distributionBinWidth 1e-4;
        resetDistributionOnWrite no;

        writeRawData    yes;
    }
    \endverbatim

    where the entries comprise:
    \table
        Property     | Description             | Required | Default value
        type         | type name: extractEulerianParticles | yes |
        faceZone     | Name of faceZone used as collection surface | yes |
        nLocations   | Number of injection bins to generate  | yes |
        aplhaName    | Name of phase indicator field | yes   |
        rhoName      | Name of density field   |      yes |
        phiNane      | Name of flux field      |      yes |
        createDistribution | Flag to create a binned distribution | yes |
        distributionBinWidth | Binned distribution bin width| yes |
        writeRawData | Flag to write raw particle data | yes |
    \endtable

SourceFiles
    extractEulerianParticles.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_extractEulerianParticles_H
#define functionObjects_extractEulerianParticles_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "runTimeSelectionTables.H"
#include "polyMesh.H"
#include "surfaceFieldsFwd.H"
#include "vectorList.H"
#include "globalIndex.H"
#include "cachedRandom.H"
#include "eulerianParticle.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
/*---------------------------------------------------------------------------*\
         Class extractEulerianParticlesFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class extractEulerianParticles
:
    public fvMeshFunctionObject,
    public writeFile
{

protected:

    // Protected data

        //- Flag to write the raw collected particle data
        bool writeRawData_;

        //- Raw particle output dictionary
        autoPtr<IOdictionary> rawParticlesDictPtr_;

        //- List to keep track of output times
        DynamicList<scalar> outputTimes_;


        // faceZone info

            //- Name of faceZone to sample
            word faceZoneName_;

            //- Index of the faceZone
            label zoneID_;

            //- Patch indices where faceZone face intersect patch
            labelList patchIDs_;

            //- Patch face indices where faceZone face intersect patch
            labelList patchFaceIDs_;

            //- Global face addressing
            globalIndex globalFaces_;


        // Field names

            //- Name of phase fraction field
            word alphaName_;

            //- Value of phase fraction used to identify particle boundaries
            scalar alphaThreshold_;

            //- Name of the velocity field, default = U
            word UName_;

            //- Name of the density field, defauls = rho
            word rhoName_;

            //- Name of the flux field, default ="rho"
            word phiName_;


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

            //- Time of last write step
            scalar lastOutputTime_;

            //- Region indices in faceZone faces from last iteration
            labelList regions0_;

            //- Number of regions from last iteration
            label nRegions0_;

            //- Particle properties (partial, being accumulated)
            List<eulerianParticle> particles_;

            //- Map from region to index in particles_ list
            Map<label> regionToParticleMap_;

            // Collected particles
            DynamicList<eulerianParticle> collectedParticles_;

            //- Minimum diameter (optional)
            //  Can be used to filter out 'small' particles
            scalar minDiameter_;

            //- Maximum diameter (optional)
            //  Can be used to filter out 'large' particles
            scalar maxDiameter_;

            //- Flag to create a distribution
            bool createDistribution_;

            //- Flag to reset the distribution on each write
            bool resetDistributionOnWrite_;

            //- Diameter distribution bin width
            scalar distributionBinWidth_;

            //- Random class needed by distribution models
            cachedRandom rndGen_;


        // Statistics

            //- Total number of collected particles
            label nCollectedParticles_;

            //- Total number of discarded particles
            label nDiscardedParticles_;

            //- Total discarded volume [m3]
            scalar discardedVolume_;


    // Protected Member Functions

        //- Return the base directory for dictionary output
        virtual fileName dictBaseFileDir() const;

        //- Check that the faceZone is valid
        virtual void checkFaceZone();

        //- Initialise the particle collection bins
        virtual void initialiseBins();

        //- Return the volumetric flux
        virtual tmp<surfaceScalarField> phiU() const;

        //- Set the blocked faces, i.e. where alpha > alpha threshold value
        virtual void setBlockedFaces
        (
            const surfaceScalarField& alphaf,
            const faceZone& fz,
            boolList& blockedFaces
        );

        //- Calculate the addressing between regions between iterations
        virtual void calculateAddressing
        (
            const label nRegionsOld,
            const label nRegionsNew,
            const scalar time,
            labelList& regionFaceIDs
        );

        //- Collect particles that have passed through the faceZone
        virtual void collectParticles
        (
            const scalar time,
            const boolList& collectParticle
        );

        //- Process latest region information
        virtual void accumulateParticleInfo
        (
            const surfaceScalarField& alphaf,
            const surfaceScalarField& phi,
            const labelList& regionFaceIDs,
            const faceZone& fz
        );

        //- Write agglomerated particle data to stream
        virtual void writeBinnedParticleData();

        template<class Type>
        inline Type faceValue
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>& field,
            const label localFaceI,
            const label globalFaceI
        ) const;

        //- Disallow default bitwise copy construct
        extractEulerianParticles(const extractEulerianParticles&);

        //- Disallow default bitwise assignment
        void operator=(const extractEulerianParticles&);


public:

    // Static data members

        //- Static data staticData
        TypeName("extractEulerianParticles");


    // Constructors

        //- Construct from components
        extractEulerianParticles
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~extractEulerianParticles();


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
    #include "extractEulerianParticlesTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
