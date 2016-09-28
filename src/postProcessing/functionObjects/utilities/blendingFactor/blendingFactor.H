/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2015 OpenFOAM Foundation
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
    Foam::blendingFactor

Group
    grpUtilitiesFunctionObjects

Description
    This function object provides information on the mode of operation of
    blended convection schemes.

    The weight of a blended scheme is given by a function of the blending
    factor, f:

        weight = f*scheme1 + (1 - f)*scheme2

    The factor is a face-based quantity, which is converted to a cell-based
    quantity by assigning the minimum blending factor for any cell face.

    An indicator (volume) field, named <functionObjectName>:<fieldName>, is
    generated that is set to (1 - f), i.e. values of:
    - 0 represent scheme1 as active, and
    - 1 represent scheme2 as active.
    - intermediate values show the contribution to scheme2

    Additional reporting is written to the standard output, providing
    statistics as to the number of cells used by each scheme.

    Example of function object specification to calculate the blending factor:
    \verbatim
    blendingFactor1
    {
        type        blendingFactor;
        functionObjectLibs ("libutilityFunctionObjects.so");

        ...

        // Name of field
        fieldName   U;
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: blendingFactor | yes       |
        phiName      | Name of flux field      | no          | phi
        fieldName    | Name of field to evaluate | yes       |
        tolerance    | Tolerance for number of blended cells | no | 0.001
        log          | Log to standard output  | no          | yes
    \endtable

SourceFiles
    blendingFactor.C
    IOblendingFactor.H

\*---------------------------------------------------------------------------*/

#ifndef blendingFactor_H
#define blendingFactor_H

#include "functionObjectState.H"
#include "functionObjectFile.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "OFstream.H"
#include "Switch.H"

#include "convectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                       Class blendingFactor Declaration
\*---------------------------------------------------------------------------*/

class blendingFactor
:
    public functionObjectState,
    public functionObjectFile
{
    // Private data

        //- Name
        const word name_;

        //- Reference to the database
        const objectRegistry& obr_;

        //- Name of flux field, default is "phi"
        word phiName_;

        //- Field name
        word fieldName_;

        //- Result field name
        word resultName_;

        //- Tolerance used when calculating the number of blended cells
        scalar tolerance_;

        //- Switch to send output to Info as well as to file
        Switch log_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        blendingFactor(const blendingFactor&);

        //- Disallow default bitwise assignment
        void operator=(const blendingFactor&);

        //- Helper function to calculate the blending factor for the scheme
        template<class Type>
        void calcScheme
        (
            const GeometricField<Type, fvPatchField, volMesh>& field,
            const typename fv::convectionScheme<Type>& cs
        );

        //- Calculate the blending factor
        template<class Type>
        void calc();


protected:

    // Protected Member Functions

        //- Write the file header
        virtual void writeFileHeader(Ostream& os) const;


public:

    //- Runtime type information
    TypeName("blendingFactor");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        blendingFactor
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~blendingFactor();


    // Member Functions

        //- Return name of the set of blendingFactor
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the blendingFactor data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Calculate the blendingFactor and write
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

#ifdef NoRepository
    #include "blendingFactorTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
