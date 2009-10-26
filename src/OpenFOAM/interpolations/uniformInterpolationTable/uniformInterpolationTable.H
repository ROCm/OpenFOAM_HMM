/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::uniformInterpolationTable

Description
    Table with uniform interval in independant variable, with linear
    interpolation

    Values specified in a dictionary:

    {
        x0              0;      // lower limit
        dx              0.2;    // fixed interval
        log10           true;   // take log(10) when interpolating?
        data                    // list of dependent data values
        (
            7870                // value at x0
            7870                // value at x0 + dx
            ...
            7870                // value at x0 + n*dx
        );
    }

SourceFiles
    uniformInterpolationTable.C

\*---------------------------------------------------------------------------*/

#ifndef uniformInterpolationTable_H
#define uniformInterpolationTable_H

#include "List.H"
#include "Switch.H"
#include "IOobject.H"
#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class uniformInterpolationTable Declaration
\*---------------------------------------------------------------------------*/

class uniformInterpolationTable
:
    public IOobject,
    public List<scalar>
{
    // Private data

        // Control parameetrs

            //- Lower limit
            scalar x0_;

            //- Fixed interval
            scalar dx_;

            //- Flag to indicate that x data is given in log10(x) form
            Switch log10_;

            //- Bound x values
            Switch bound_;


    // Private Member Functions

        //- Check that the table is valid
        void checkTable() const;

        //- Disallow default bitwise assignment
        void operator=(const uniformInterpolationTable&);


public:

    //- Runtime type information
    TypeName("uniformInterpolationTable");


    // Constructors

        //- Construct from IOobject and readFields flag. Creates a null object
        //  if readFields = false
        uniformInterpolationTable(const IOobject& io, const bool readFields);

        //- Construct from name, objectRegistry and dictionary.
        //  If initialiseOnly flag is set, control parameters are read from
        //  the dictionary, but not the data table
        uniformInterpolationTable
        (
            const word& tableName,
            const objectRegistry& db,
            const dictionary& dict,
            const bool initialiseOnly = false
        );

        //- Construct as copy
        uniformInterpolationTable(const uniformInterpolationTable& uit);


    //- Destructor
    ~uniformInterpolationTable();


    // Member Functions

        // Access

            //- Return the lower limit
            inline scalar x0() const;

            //- Return the fixed interval
            inline scalar dx() const;

            //- Return the log10(x) flag
            inline const Switch& log10() const;

            //- Return the bound flag
            inline const Switch& bound() const;


        // Edit

            //- Return the lower limit
            inline scalar& x0();

            //- Return the fixed interval
            inline scalar& dx();

            //- Return the log10(x) flag
            inline Switch& log10();

            //- Return the bound flag
            inline Switch& bound();


        // Evaluation

            //- Return the minimum x value
            inline scalar xMin() const;

            //- Return the maximum x value
            inline scalar xMax() const;

            //- Interpolate
            scalar interpolate(scalar x) const;

            //- Interpolate - takes log10 flag into account
            scalar interpolateLog10(scalar x) const;


        // Override ancestor size() function and [] operator

            //- Return the size of the table
            label size() const
            {
                return List<scalar>::size();
            }

            //- Use List[] operator for read access
            scalar operator[](label x) const
            {
                return List<scalar>::operator[](x);
            }

            //- Use List[] operator for write access
            scalar& operator[](label x)
            {
                return List<scalar>::operator[](x);
            }


        // I-O

            //- Write
            void write() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "uniformInterpolationTableI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
