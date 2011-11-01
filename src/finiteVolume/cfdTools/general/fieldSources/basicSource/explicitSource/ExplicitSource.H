/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::ExplicitSource

Description
    Time activated explicit source.

    Sources described by:

    {
        active          true;      // on/off switch
        timeStart       0.2;       // start time
        duration        2.0;       // duration
        selectionMode   points;    // cellSet/cellZone/all
        volumeMode      absolute;  // specific

        fieldData                  // field data - usage for multiple fields
        (
            (H2O 0.005)
        );

        fieldData       0.005;     // field data - usage for single field

        points                     // list of points when selectionMode = points
        (
            (2.75 0.5 0)
        );

        cellSet         c0;        // cellSet name when selectionMode=cellSet
        cellZone        c0;        // cellZone name when selectionMode=cellZone
    }

SourceFiles
    ExplicitSource.C

\*---------------------------------------------------------------------------*/

#ifndef ExplicitSource_H
#define ExplicitSource_H

#include "Tuple2.H"
#include "volFieldsFwd.H"
#include "DimensionedField.H"
#include "basicSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

class fvMesh;

template<class Type>
class ExplicitSource;

// Forward declaration of friend functions

template<class Type>
Ostream& operator<<
(
    Ostream&,
    const ExplicitSource<Type>&
);

/*---------------------------------------------------------------------------*\
               Class ExplicitSource Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class ExplicitSource
:
    public basicSource
{
public:

    // Public data

        //- Enumeration for volume types
        enum volumeModeType
        {
            vmAbsolute,
            vmSpecific
        };

        //- Word list of volume mode type names
        static const wordList volumeModeTypeNames_;


protected:

    // Protected data

        typedef Tuple2<word, Type> fieldNameValuePair;

        //- Volume mode
        volumeModeType volumeMode_;

        //- Volume normalisation
        scalar VDash_;

        //- Source value per field
        List<fieldNameValuePair> fieldData_;


    // Protected functions

        //- Helper function to convert from a word to a volumeModeType
        volumeModeType wordToVolumeModeType(const word& vtName) const;

        //- Helper function to convert from a volumeModeType to a word
        word volumeModeTypeToWord(const volumeModeType& vtType) const;

        //- Set the local field data
        void setFieldData(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("ExplicitSource");


    // Constructors

        //- Construct from components
        ExplicitSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Member Functions

        // Check

            //- Return index of field name if found in fieldNames list
            virtual label applyToField(const word& fieldName) const;


        // Access

            //- Return const access to the volume mode
            inline const volumeModeType& volumeMode() const;

            //- Return const access to the source field value
            inline const Type& fieldData() const;


        // Edit

            //- Return access to the volume mode
            inline volumeModeType& volumeMode();

            //- Return access to the source field value
            inline Type& fieldData();


        // Evaluation

            //- Add explicit contribution to equation
            virtual void addSup(fvMatrix<Type>& eqn, const label fieldI);


        // I-O

            //- Write the source properties
            virtual void writeData(Ostream&) const;

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ExplicitSource.C"
#   include "ExplicitSourceIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ExplicitSourceI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
