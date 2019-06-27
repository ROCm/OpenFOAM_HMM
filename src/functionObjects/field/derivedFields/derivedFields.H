/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::functionObjects::derivedFields

Group
    grpFieldFunctionObjects

Description
    A limited set of predefined derived fields ("rhoU", "pTotal").

    \verbatim
    derived
    {
        type    derivedFields;
        libs    ("libfieldFunctionObjects.so");

        fields  (rhoU pTotal);

        // Optional: reference density for incompressible
        rhoRef  1.25;
    }
    \endverbatim

    Entries:
    \table
        Property | Description                              | Required | Default
        type     | derivedFields                            | yes |
        derived  | Derived fields (pTotal/rhoU)             | yes |
        rhoRef   | Reference density (incompressible)       | no  | 1
    \endtable

    The known derived fields
    \plaintable
       rhoU          | (rho * U)
       pTotal        | (p + 1/2 * rho * U)
    \endplaintable

SourceFiles
    derivedFields.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_derivedFields_H
#define functionObjects_derivedFields_H

#include "fvMeshFunctionObject.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class derivedFields Declaration
\*---------------------------------------------------------------------------*/

class derivedFields
:
    public fvMeshFunctionObject
{
public:

    // Public Enumerations

        //- Derived/calculated field type
        enum derivedType
        {
            NONE = 0,           //!< "none"
            MASS_FLUX,          //!< "rhoU"
            TOTAL_PRESSURE,     //!< "pTotal"
            UNKNOWN
        };

        //- Known derived field types
        static const Enum<derivedType> knownNames;


protected:

    // Read from dictionary

        //- List of derived field (types) to create
        List<derivedType> derivedTypes_;

        //- Reference density (to convert from kinematic to static pressure)
        scalar rhoRef_;


    // Private Member Functions

        //- Hard-coded derived field (rho * U)
        //  \return true if field did not previously exist
        bool add_rhoU(const word& derivedName);

        //- Hard-coded derived field (p + 1/2 * rho * U)
        //  \return true if field did not previously exist
        bool add_pTotal(const word& derivedName);


public:

    //- Run-time type information
    TypeName("derivedFields");


    // Constructors

        //- Construct from Time and dictionary
        derivedFields
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~derivedFields() = default;


    // Member Functions

        //- Remove (checkOut) derived fields from the object registry
        void removeDerivedFields();

        //- Read the data
        virtual bool read(const dictionary& dict);

        //- Calculate the derived fields
        virtual bool execute();

        //- Write derived fields
        virtual bool write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Update for mesh point-motion
        virtual void movePoints(const polyMesh& m);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
