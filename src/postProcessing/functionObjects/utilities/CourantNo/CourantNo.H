/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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
    Foam::CourantNo

Group
    grpUtilitiesFunctionObjects

Description
    This function object calculates and outputs the Courant number as a
    volScalarField.  The field is stored on the mesh database so that it can
    be retrieved and used for other applications.

    Example of function object specification to calculate the Courant number:
    \verbatim
    CourantNo1
    {
        type        CourantNo;
        functionObjectLibs ("libutilityFunctionObjects.so");
        ...
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: CourantNo    | yes         |
        rhoName      | Name of density field   | no          | rho
        phiName      | Name of flux field      | no          | phi
        resultName   | Name of Courant number field | no     | <function name>
        log          | Log to standard output  | no          | yes
    \endtable


SourceFiles
    CourantNo.C
    IOCourantNo.H

\*---------------------------------------------------------------------------*/

#ifndef CourantNo_H
#define CourantNo_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                          Class CourantNo Declaration
\*---------------------------------------------------------------------------*/

class CourantNo
{
    // Private data

        //- Name of this set of CourantNo objects
        word name_;

        //- Reference to the database
        const objectRegistry& obr_;

        //- On/off switch
        bool active_;

        //- Name of flux field, default is "phi"
        word phiName_;

        //- Name of density field (optional)
        word rhoName_;

        //- Result name
        word resultName_;

        //- Switch to send output to Info as well as to file
        Switch log_;


    // Private Member Functions

        //- Divide Co by rho if required
        tmp<volScalarField::DimensionedInternalField> byRho
        (
            const tmp<volScalarField::DimensionedInternalField>& Co
        ) const;

        //- Disallow default bitwise copy construct
        CourantNo(const CourantNo&);

        //- Disallow default bitwise assignment
        void operator=(const CourantNo&);


public:

    //- Runtime type information
    TypeName("CourantNo");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        CourantNo
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~CourantNo();


    // Member Functions

        //- Return name of the set of CourantNo
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the CourantNo data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Calculate the CourantNo and write
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

#endif

// ************************************************************************* //
