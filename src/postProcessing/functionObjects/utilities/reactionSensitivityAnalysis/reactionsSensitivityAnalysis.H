/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd
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
    Foam::reactionsSensitivityAnalysis

Group
    grpUtilitiesFunctionObjects

Description
    This function object creates four data files named:

    "consumption"    :   consumption rate
    "production"     :   destruction rate
    "productionInt"  :   integral between dumps of the production rate
    "consumptionInt" :   integral between dumps of the consumption rate

    The function object indicates reaction rates of creation or destruction
    of spcecies in each reaction.


SourceFiles
    reactionsSensitivityAnalysis.C
    IOreactionsSensitivityAnalysis.H

\*---------------------------------------------------------------------------*/

#ifndef reactionsSensitivityAnalysis_H
#define reactionsSensitivityAnalysis_H

#include "functionObjectFile.H"
#include "volFields.H"
#include "basicChemistryModel.H"
#include "autoPtr.H"
#include "basicMultiComponentMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;

/*---------------------------------------------------------------------------*\
                    Class reactionsSensitivityAnalysis Declaration
\*---------------------------------------------------------------------------*/
template<class chemistryType>
class reactionsSensitivityAnalysis
:
    public functionObjectFile

{
    // Private data

        //- Name of this set of reactionsSensitivityAnalysis objects
        word name_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- On/off switch
        bool active_;

        //- List list for species production
        scalarListList production_;

        //- List list for species consumption
        scalarListList consumption_;

        //- List list for species production integral
        scalarListList productionInt_;

        //- List list for species consumption integral
        scalarListList consumptionInt_;

        //- Start time of integration
        scalar startTime_;

        //- End time of integration
        scalar endTime_;

        //- Word list of species
        wordList speciesNames_;

        //-Number of reactions
        label nReactions_;


         // File streams

            //- Integrated coefficients
            autoPtr<OFstream> prodFilePtr_;

            //- Moment coefficient
            autoPtr<OFstream> consFilePtr_;

            //- Drag coefficient
            autoPtr<OFstream> prodIntFilePtr_;

            //- Lift coefficient
            autoPtr<OFstream> consIntFilePtr_;



    // Private Member Functions


        //- Create file names for forces and bins
        void createFileNames();

        //- Output file header information
        void writeFileHeader(OFstream& os);

        //- Calculate production and destruction of each species
        void calculateSpeciesRR(const basicChemistryModel&);

        //- Write species production/consumption rates
        void writeSpeciesRR();


        //- Disallow default bitwise copy construct
        reactionsSensitivityAnalysis(const reactionsSensitivityAnalysis&);

        //- Disallow default bitwise assignment
        void operator=(const reactionsSensitivityAnalysis&);


public:

    //- Runtime type information
    TypeName("reactionsSensitivityAnalysis");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        reactionsSensitivityAnalysis
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~reactionsSensitivityAnalysis();


    // Member Functions

        //- Return name of the set of reactionsSensitivityAnalysis
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the reactionsSensitivityAnalysis data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Calculate the reactionsSensitivityAnalysis and write
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
#   include "reactionsSensitivityAnalysis.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
