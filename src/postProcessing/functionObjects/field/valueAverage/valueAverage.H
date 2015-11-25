/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::valueAverage

Group
    grpFieldFunctionObjects

Description
    This function object calculates the average value from the output of
    function objects that generate singular values.

    Example of function object specification:
    \verbatim
    valueAverage1
    {
        type        valueAverage;
        functionObjectLibs ("libfieldFunctionObjects.so");
        ...
        writeToFile yes;
        log         yes;
        functionObjectName forceCoeffs1;
        fields      (Cm Cd Cl);
        window      0.5;
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: valueAverage | yes         |
        writeToFile  | write average data to file |  no      | yes
        log          | write average data to standard output | no | yes
        fields       | list of fields to process | yes       |
    \endtable

    Output data is written to the file \<timeDir\>/valueAverage.dat

SeeAlso
    Foam::functionObject
    Foam::functionObjectFile
    Foam::functionObjectState
    Foam::OutputFilterFunctionObject

SourceFiles
    valueAverage.C
    valueAverageTemplates.C
    IOvalueAverage.H

\*---------------------------------------------------------------------------*/

#ifndef valueAverage_H
#define valueAverage_H

#include "functionObjectState.H"
#include "functionObjectFile.H"
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
                        Class valueAverage Declaration
\*---------------------------------------------------------------------------*/

class valueAverage
:
    public functionObjectState,
    public functionObjectFile
{
protected:

    // Protected data

        //- Reference to the database
        const objectRegistry& obr_;

        //- Name of function object to retrueve data from
        word functionObjectName_;

        //- List of fields on which to operate
        wordList fieldNames_;

        //- Averaging window
        const scalar window_;

        //- Average time per field
        List<scalar> totalTime_;

        //- Reset the averaging process on restart flag
        Switch resetOnRestart_;

        //- Switch to send output to Info as well
        Switch log_;


    // Protected Member Functions

        //- Templated function to calculate the average
        template<class Type>
        void calc
        (
            const word& fieldName,
            const word& meanName,
            const scalar alpha,
            const scalar beta,
            bool& processed
        );

        //- Output file header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Disallow default bitwise copy construct
        valueAverage(const valueAverage&);

        //- Disallow default bitwise assignment
        void operator=(const valueAverage&);


public:

    //- Runtime type information
    TypeName("valueAverage");

    //- Constructor
    valueAverage
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict,
       const bool loadFromFiles = false
    );

    //- Destructor
    virtual ~valueAverage();


    // Public Member Functions

        //- Read the field min/max data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Write the fieldMinMax
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
    #include "valueAverageTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
