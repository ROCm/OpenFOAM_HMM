/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
    Foam::functionObjects::histogram

Group
    grpFieldFunctionObjects

Description
    Computes the volume-weighted histogram of an input \c volScalarField.

    Operands:
    \table
      Operand        | Type              | Location
      input          | volScalarField    | $FOAM_CASE/\<time\>/\<inpField\>
      output file    | dat  | $FOAM_CASE/postProcessing/\<FO\>/\<time\>/\<file\>
      output field   | -                 | -
    \endtable

    The set written contains two columns, the first the volume averaged values,
    the second the raw bin count.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    histogram1
    {
        // Mandatory entries (unmodifiable)
        type        histogram;
        libs        (fieldFunctionObjects);

        // Mandatory (inherited) entries (runtime modifiable)
        field       p;
        nBins       100;
        setFormat   gnuplot;

        // Optional entries (runtime modifiable)
        max         5;
        min        -5;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: histogram               | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      field        | Name of operand field              | word |  yes  | -
      nBins        | Number of histogram bins           | label | yes  | -
      setFormat    | Output format                      | word |  yes  | -
      max          | Maximum value sampled       | scalar | no  | fieldMax
      min          | minimum value sampled       | scalar | no  | 0.0
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link writeFile.H \endlink

    Usage by the \c postProcess utility is not available.

Note
    If \c max is not provided it will use the field's min and max as the bin
    extremes. If \c max is provided but not \c min it will use 0.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::writeFile
    - ExtendedCodeGuide::functionObjects::field::histogram

SourceFiles
    histogram.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_histogram_H
#define functionObjects_histogram_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "writer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class histogram Declaration
\*---------------------------------------------------------------------------*/

class histogram
:
    public fvMeshFunctionObject,
    public writeFile
{
    // Private Data

        //- Number of bins
        label nBins_;

        //- Name of field
        word fieldName_;

        //- Maximum value
        scalar max_;

        //- Minimum value
        scalar min_;

        //- Output formatter to write
        autoPtr<writer<scalar>> formatterPtr_;


    // Private Member Functions

        void writeGraph
        (
            const coordSet& coords,
            const word& valueName,
            const scalarField& normalizedValues,
            const scalarField& absoluteValues
        ) const;


public:

    //- Runtime type information
    TypeName("histogram");


    // Constructors

        //- Construct from Time and dictionary
        histogram
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        histogram(const histogram&) = delete;

        //- No copy assignment
        void operator=(const histogram&) = delete;


    // Destructor
    virtual ~histogram() = default;


    // Member Functions

        //- Read the histogram data
        virtual bool read(const dictionary&);

        //- Execute (effectively no-op)
        virtual bool execute();

        //- Calculate the histogram and write
        //  postProcess overrides the usual writeControl behaviour and
        //  forces writing always (used in post-processing mode)
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
