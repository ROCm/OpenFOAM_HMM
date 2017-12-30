/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::multiLevelDecomp

Description
    Decompose given using consecutive application of decomposers.

SourceFiles
    multiLevelDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef multiLevelDecomp_H
#define multiLevelDecomp_H

#include "decompositionMethod.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class multiLevelDecomp Declaration
\*---------------------------------------------------------------------------*/

class multiLevelDecomp
:
    public decompositionMethod
{
    // Private data

        //- Original coefficients for this method
        const dictionary& coeffsDict_;

        //- Rewritten dictionary of individual methods
        dictionary methodsDict_;

        PtrList<decompositionMethod> methods_;


    // Private Member Functions

        //- Fill the methodsDict_
        void createMethodsDict();

        //- Set methods based on the contents of the methodsDict_
        void setMethods();


        //- Given connectivity across processors work out connectivity
        //  for a (consistent) subset
        void subsetGlobalCellCells
        (
            const label nDomains,
            const label domainI,
            const labelList& dist,

            const labelListList& cellCells,
            const labelList& set,
            labelListList& subCellCells,
            labelList& cutConnections
        ) const;

        //- Decompose at 'currLevel' without addressing
        void decompose
        (
            const labelListList& pointPoints,
            const pointField& points,
            const scalarField& pointWeights,
            const labelUList& pointMap,  // map back to original points
            const label currLevel,
            const label leafOffset,

            labelList& finalDecomp
        );

        //- Disallow default bitwise copy construct and assignment
        void operator=(const multiLevelDecomp&) = delete;
        multiLevelDecomp(const multiLevelDecomp&) = delete;


public:

    //- Runtime type information
    TypeName("multiLevel");


    // Constructors

        //- Construct given the decomposition dictionary
        multiLevelDecomp(const dictionary& decompDict);

        //- Construct given decomposition dictionary and region name
        multiLevelDecomp
        (
            const dictionary& decompDict,
            const word& regionName
        );


    //- Destructor
    virtual ~multiLevelDecomp()
    {}


    // Member Functions

        //- Is method parallel aware (i.e. does it synchronize domains across
        //  proc boundaries)
        virtual bool parallelAware() const;

        //- Inherit decompose from decompositionMethod
        using decompositionMethod::decompose;

        //- Return for every coordinate the wanted processor number. Use the
        //  mesh connectivity (if needed)
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const pointField& points,
            const scalarField& pointWeights
        );

        //- Return for every coordinate the wanted processor number. Explicitly
        //  provided connectivity - does not use mesh_.
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cc,
            const scalarField& cWeights
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
