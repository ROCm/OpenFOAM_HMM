/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::functionObjects::fieldVisualisationBase

Description
    Base class for scene objects

SourceFiles
    fieldVisualisationBase.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fieldVisualisationBase_H
#define functionObjects_fieldVisualisationBase_H

#include "dictionary.H"
#include "Tuple2.H"
#include "Enum.H"
#include "vector.H"
#include "HashPtrTable.H"
#include "Function1.H"

#include "vtkSmartPointer.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class vtkActor;
class vtkLookupTable;
class vtkPolyData;
class vtkPolyDataMapper;
class vtkRenderer;


class vtkMapper;

namespace Foam
{
namespace functionObjects
{

class runTimePostProcessing;

/*---------------------------------------------------------------------------*\
                   Class fieldVisualisationBase Declaration
\*---------------------------------------------------------------------------*/

class fieldVisualisationBase
{
public:

    // Public enumerations

        enum colourByType
        {
            cbColour,
            cbField
        };

        static const Enum<colourByType> colourByTypeNames;

        enum colourMapType
        {
            cmRainbow,
            cmBlueWhiteRed,
            cmFire,
            cmGreyscale
        };

        static const Enum<colourMapType> colourMapTypeNames;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        fieldVisualisationBase(const fieldVisualisationBase&);

        //- Disallow default bitwise assignment
        void operator=(const fieldVisualisationBase&);


protected:

    // Protected data

        struct scalarBar
        {
            bool visible_;
            bool vertical_;
            Tuple2<scalar, scalar> position_;
            string title_;
            label fontSize_;
            string labelFormat_;
            label numberOfLabels_;
        };

        //- Colours
        const HashPtrTable<Function1<vector>, word>& colours_;

        //- Field name
        word fieldName_;

        //- Colour by type
        colourByType colourBy_;

        //- Colour map type
        colourMapType colourMap_;

        //- Range of values
        Tuple2<scalar, scalar> range_;

        //- Scalar bar
        scalarBar scalarBar_;


    // Protected Member Functions

        //- Set the colour map
        void setColourMap(vtkLookupTable* lut) const;

        //- Add scalar bar to renderer
        void addScalarBar
        (
            const scalar position,
            vtkRenderer* renderer,
            vtkLookupTable* lut
        ) const;

        //- Set field/configure mapper, add scalar bar
        void setField
        (
            const scalar position,
            const word& colourFieldName,
            vtkPolyDataMapper* mapper,
            vtkRenderer* renderer,
            vtkPolyData* pData
        ) const;

        //- Add glyphs
        void addGlyphs
        (
            const scalar position,
            const word& scaleFieldName,
            const word& colourFieldName,
            const scalar maxGlyphLength,
            vtkPolyData* data,
            vtkActor* actor,
            vtkRenderer* renderer
        ) const;


public:

    // Constructors

        //- Construct from dictionary
        fieldVisualisationBase
        (
            const dictionary& dict,
            const HashPtrTable<Function1<vector>, word>& colours
        );


    //- Destructor
    virtual ~fieldVisualisationBase();


    // Member Functions

        // Access

            //- Return the colours
            const HashPtrTable<Function1<vector>, word>& colours() const;

            //- Return the field name
            const word& fieldName() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
