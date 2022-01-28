/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2022 OpenCFD Ltd.
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
    Foam::gltfSetWriter

Description
    Writes point data in glTF v2 format

    Two files are generated:
    - filename.bin  : a binary file containing all scene entities
    - filename.gltf : a JSON file that ties fields to the binary data

    The output can contain both geometry and fields, with additional support
    for colours using a user-supplied colour map, and animation of particle
    tracks.

    Controls are provided via the optional formatOptions dictionary.

    For non-particle track data:

    \verbatim
    formatOptions
    {
        // Apply colours flag (yes | no ) [optional]
        colours     yes;

        // List of options per field
        fieldInfo
        {
            p
            {
                // Colour map [optional]
                colourMap       <colourMap>;

                // Colour map minimum and maximum limits [optional]
                // Uses field min and max if not specified
                min             0;
                max             1;

                // Alpha channel [optional] (<scalar>)
                alpha           0.5;
            }
        }
    }
    \verbatim

    For particle tracks:

    \verbatim
    formatOptions
    {
        // Apply colours flag (yes | no) [optional]
        colours     yes;

        // Animate tracks (yes | no) [optional]
        animate     yes;

        // Animation properties [optional]
        animationInfo
        {
            // Colour map [optional]
            colourMap       <colourMap>;

            // Colour [optional] (<vector> | uniform | field)
            colour          (1 0 0); // RGB in range [0-1]

            //colour          uniform;
            //colourValue     (1 0 0); // RGB in range [0-1]

            //colour          field;
            //colourField     d;

            // Colour map minimum and maximum limits [optional]
            // Note: for colour = field option
            // Uses field min and max if not specified
            min             0;
            max             1;

            // Alpha channel [optional] (<scalar>)
            alpha           0.5;
        }
    }
    \endverbatim

Note
    When writing particle animations, the particle field and colour properties
    correspond to initial particle state (first data point) and cannot be
    animated (limitation of the file format).

    For more information on the specification see
    https://www.khronos.org/registry/glTF/

SourceFiles
    gltfCoordSetWriter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_coordSetWriters_gltfWriter_H
#define Foam_coordSetWriters_gltfWriter_H

#include "writer.H"
#include "colourTable.H"
#include "MinMax.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class gltfSetWriter Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class gltfSetWriter
:
    public writer<Type>
{
public:

    // Enumerations

        //- Field option used for colours
        enum class fieldOption : char
        {
            NONE,       //!< Placeholder type (unnamed)
            UNIFORM,    //!< Uniform value
            FIELD       //!< Field value
        };


    //- Strings corresponding to the field options
    static const Enum<fieldOption> fieldOptionNames_;


private:

    // Private Data

        //- Flag to animate - for particle tracks only
        bool animate_;

        //- Flag to add field colours
        bool colour_;

        //- Local field information
        const dictionary fieldInfoDict_;

        //- Animation information
        const dictionary animationDict_;


    // Private Member Functions

        //- Return the colour map name
        word getColourMap(const dictionary& dict) const;

        //- Return the colour table corresponding to the colour map
        const colourTable& getColourTable(const dictionary& dict) const;

        //- Return the named min/max field limits (from sub-dictionary)
        scalarMinMax getFieldLimits(const word& fieldName) const;

        //- Return the alpha field for mesh values
        tmp<scalarField> getAlphaField(const dictionary& dict) const;

        //- Return the animation colour when animating tracks
        vector getTrackAnimationColour
        (
            const colourTable& colours,
            const wordList& valueSetNames,
            const List<List<Field<Type>>>& valueSets,
            const label tracki
        ) const;


public:

    //- Runtime type information
    TypeName("gltf");


    // Constructors

        //- Default construct
        gltfSetWriter();

        //- Construct from dictionary
        explicit gltfSetWriter(const dictionary& dict);


    //- Destructor
    virtual ~gltfSetWriter() = default;


    // Member Functions

        //- Return the file name
        virtual fileName getFileName
        (
            const coordSet&,
            const wordList&
        ) const;

        //- Write
        virtual void write
        (
            const coordSet&,
            const wordList&,
            const List<const Field<Type>*>&,
            Ostream&
        ) const;

        //- Write tracks (main entry point)
        virtual void write
        (
            const bool writeTracks,
            const List<scalarField>& times,
            const PtrList<coordSet>&,
            const wordList& valueSetNames,
            const List<List<Field<Type>>>&,
            Ostream&
        ) const;

        //- Write animated tracks
        virtual void writeAnimateTracks
        (
            const bool writeTracks,
            const List<scalarField>& times,
            const PtrList<coordSet>& tracks,
            const wordList& valueSetNames,
            const List<List<Field<Type>>>& valueSets,
            Ostream&
        ) const;

        //- Write static tracks
        virtual void writeStaticTracks
        (
            const bool writeTracks,
            const List<scalarField>& times,
            const PtrList<coordSet>& tracks,
            const wordList& valueSetNames,
            const List<List<Field<Type>>>& valueSets,
            Ostream&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "gltfCoordSetWriter.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
