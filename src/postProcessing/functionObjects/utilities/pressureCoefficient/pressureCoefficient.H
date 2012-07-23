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
    Foam::pressureCoefficient

Group
    grpUtilitiesFunctionObjects

Description
    Calculates pressure coefficient, \f$c_p\f$

        \f[
            c_p = \frac{p}{p_{dyn,inf}}
        \f]
    where:
        \f[
            p_{dyn,inf} = 0.5 \rho |U_{inf}|^2
        \f]

    where
    \vartable
        c_p     | pressure coefficient
        p       | pressure [bar]
        \rho    | density [kg/m3]
        U       | velocity [m/s]
    \endvartable


SourceFiles
    pressureCoefficient.C
    IOpressureCoefficient.H

\*---------------------------------------------------------------------------*/

#ifndef pressureCoefficient_H
#define pressureCoefficient_H

#include "volFieldsFwd.H"
#include "pointFieldFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                  Class pressureCoefficient Declaration
\*---------------------------------------------------------------------------*/

class pressureCoefficient
{
    // Private data

        //- Name of this set of pressureCoefficient objects
        word name_;

        const objectRegistry& obr_;

        //- on/off switch
        bool active_;

        //- Name of pressure field, default is "p"
        word pName_;

        //- Name of density field (optional)
        word rhoName_;

        //- Free stream velocity magnitude [m/s]
        scalar magUinf_;


    // Private Member Functions

        //- Return 1 if the pressure field is kinematic, i.e. p/rho
        //  otherwise return rho from database
        tmp<volScalarField> rho(const volScalarField& p) const;

        //- Disallow default bitwise copy construct
        pressureCoefficient(const pressureCoefficient&);

        //- Disallow default bitwise assignment
        void operator=(const pressureCoefficient&);


public:

    //- Runtime type information
    TypeName("pressureCoefficient");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        pressureCoefficient
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~pressureCoefficient();


    // Member Functions

        //- Return name of the set of pressureCoefficient
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the pressureCoefficient data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Calculate the pressureCoefficient and write
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const pointField&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
