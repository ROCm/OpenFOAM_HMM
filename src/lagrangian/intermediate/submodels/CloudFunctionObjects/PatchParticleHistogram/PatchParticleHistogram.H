/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::PatchParticleHistogram

Description
    Computes a histogram for the distribution of particle diameters
    and corresponding number of particles hitting on a given list of patches.

    Operands:
    \table
      Operand        | Type | Location
      input          | -    | -
      output file    | dat  | $FOAM_CASE/postProcessing/\<FO\>/\<time\>/\<file\>
      output field   | -    | -
    \endtable

    The output file contains two columns, the first is the bin edges of
    the particle diameter (i.e. \c d), and the second is the number of
    particles whose diameter falling into the corresponding bin
    (i.e. \c nParticles).

Usage
    Minimal example by using
    \c constant/reactingCloud1Properties.cloudFunctions:
    \verbatim
    patchParticleHistogram1
    {
        // Mandatory entries (unmodifiable)
        type             patchParticleHistogram;
        patches          (<patch1> <patch2> ... <patchN>);
        nBins            10;
        min              0.1;
        max              10.0;
        maxStoredParcels 20;
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type   | Reqd  | Dflt
      type         | Type name: patchParticleHistogram  | word   |  yes  | -
      patches      | Names of operand patches         | wordList |  yes  | -
      nBins        | Number of histogram bins           | label  |  yes  | -
      max          | Maximum value of histogram data    | scalar |  yes  | -
      min          | Minimum value of histogram data    | scalar |  yes  | -
      maxStoredParcels | Maximum number of parcels to process | label | yes | -
    \endtable

Note
  - The underlying type of \c maxStoredParcels is set as a scalar for I/O.

SourceFiles
    PatchParticleHistogram.C

\*---------------------------------------------------------------------------*/

#ifndef PatchParticleHistogram_H
#define PatchParticleHistogram_H

#include "CloudFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class PatchParticleHistogram Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PatchParticleHistogram
:
    public CloudFunctionObject<CloudType>
{
    // Private Data

        //- Convenience typedef for parcel type
        typedef typename CloudType::particleType parcelType;

        //- Number of data bins
        const label nBins_;

        //- Minimum value of histogram data
        const scalar min_;

        //- Maximum value of histogram data
        const scalar max_;

        //- Bin width of histogram
        const scalar delta_;

        //- Maximum number of parcels to store - set as a scalar for I/O
        const scalar maxStoredParcels_;

        //- Bin edges of histogram
        scalarField binEdges_;

        //- List of patch indices to post-process
        labelList patchIDs_;

        //- Accumulated number of particles per patch
        //- binned according to the histogram settings
        List<List<scalar>> nParticlesCumulative_;

        //- List of time for each data record
        List<DynamicList<scalar>> times_;

        // List of patch-hit particle diameters
        List<DynamicList<scalar>> patchDiameters_;

        // List of number of patch-hit particles
        List<DynamicList<scalar>> patchParticles_;


    // Private Member Functions

        //- Return local patchi if patch is in patchIds_ list
        label applyToPatch(const label globalPatchi) const;


protected:

    // Protected Member Functions

        //- Write post-processing info
        void write();


public:

    //- Runtime type information
    TypeName("patchParticleHistogram");


    // Constructors

        //- No default construct
        PatchParticleHistogram() = delete;

        //- Construct from dictionary
        PatchParticleHistogram
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Copy construct
        PatchParticleHistogram(const PatchParticleHistogram<CloudType>& ppm);

        //- No copy assignment
        void operator=(const PatchParticleHistogram<CloudType>&) = delete;

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new PatchParticleHistogram<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~PatchParticleHistogram() = default;


    // Member Functions

        // Evaluation

            //- Post-patch hook
            virtual void postPatch
            (
                const parcelType& p,
                const polyPatch& pp,
                bool& keepParticle
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PatchParticleHistogram.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
