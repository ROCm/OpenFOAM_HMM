/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::fixedTemperatureSource

Description
    Fixed temperature equation constraint

    Sources described by:

        fixedTemperatureSourceCoeffs
        {
            fieldNames      (h e hs);   // names of fields to apply source
            mode            uniform;    // constant or lookup

            // uniform option
            temperature     constant 500; // fixed temperature [K]

            // loolup option
            // TName        T;          // optional temperature field name
        }

Note:
    The 'uniform' option allows the use of a time-varying uniform temperature
    by means of the DataEntry type.

SourceFiles
    basicSource.C

\*---------------------------------------------------------------------------*/

#ifndef fixedTemperatureSource_H
#define fixedTemperatureSource_H

#include "basicSource.H"
#include "NamedEnum.H"
#include "DataEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class fixedTemperatureSource Declaration
\*---------------------------------------------------------------------------*/

class fixedTemperatureSource
:
    public basicSource
{

public:

    //- Temperature mode
    enum temperatureMode
    {
        tmUniform,
        tmLookup
    };


    //- String representation of temperatureMode enums
    static const NamedEnum<temperatureMode, 2> temperatureModeNames_;


protected:

    // Protected data

        //- Operation mode
        temperatureMode mode_;

        //- Uniform temperature [K]
        autoPtr<DataEntry<scalar> > Tuniform_;

        //- Temperature field name
        word TName_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        fixedTemperatureSource(const fixedTemperatureSource&);

        //- Disallow default bitwise assignment
        void operator=(const fixedTemperatureSource&);


public:

    //- Runtime type information
    TypeName("fixedTemperatureSource");


    // Constructors

        //- Construct from components
        fixedTemperatureSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~fixedTemperatureSource()
    {}


    // Member Functions

        virtual bool alwaysApply() const;


        // Set values directly

            //- Scalar
            virtual void setValue(fvMatrix<scalar>& eqn, const label fieldI);


        // I-O

            //- Write data
            virtual void writeData(Ostream&) const;

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
