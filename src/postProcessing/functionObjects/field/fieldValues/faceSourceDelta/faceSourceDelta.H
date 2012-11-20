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

Class
    Foam::fieldValues::faceSourceDelta

Group
    grpFieldFunctionObjects

Description
    This function object provides a differencing option between two 'face
    source' function objects.

    Example of function object specification:
    \verbatim
    faceSourceDelta1
    {
        type            faceSourceDelta;
        functionObjectLibs ("libfieldFunctionObjects.so");
        faceSource1
        {
            ...
        }
        faceSource2
        {
            ...
        }
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: faceSourceDelta   | yes    |
    \endtable

SeeAlso
    Foam::faceSource

SourceFiles
    faceSourceDelta.C

\*---------------------------------------------------------------------------*/

#ifndef faceSourceDelta_H
#define faceSourceDelta_H

#include "functionObjectFile.H"
#include "faceSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fieldValues
{

/*---------------------------------------------------------------------------*\
                       Class faceSourceDelta Declaration
\*---------------------------------------------------------------------------*/

class faceSourceDelta
:
    public functionObjectFile
{

private:

    // Private data

        //- Database this class is registered to
        const objectRegistry& obr_;

        //- Switch to send output to Info as well as to file
        Switch log_;

        //- Face source 1
        faceSource faceSource1_;

        //- Face source 2
        faceSource faceSource2_;


    // Private Member Functions

        //- Templated function to process common fields
        template<class Type>
        void processFields(bool& found);


protected:

    // Functions to be over-ridden from IOoutputFilter class

        //- Update mesh
        virtual void updateMesh(const mapPolyMesh&);

        //- Move points
        virtual void movePoints(const Field<point>&);

        //- Output file header information
        virtual void writeFileHeader(const label i);


public:

    //- Run-time type information
    TypeName("faceSourceDelta");


    //- Construct from components
    faceSourceDelta
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict,
        const bool loadFromFiles = false
    );


    //- Destructor
    virtual ~faceSourceDelta();


    // Public Member Functions

        // Function object functions

            //- Read from dictionary
            virtual void read(const dictionary&);

            //- Calculate and write
            virtual void write();

            //- Execute
            virtual void execute();

            //- Execute the at the final time-loop, currently does nothing
            virtual void end();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fieldValues
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "faceSourceDeltaTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
