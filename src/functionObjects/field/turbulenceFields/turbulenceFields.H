/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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
    Foam::functionObjects::turbulenceFields

Group
    grpFieldFunctionObjects

Description
    Computes various turbulence-related quantities that are not typically
    output during calculations, and stores/writes them on the mesh database
    for further manipulation.

    Fields are stored as copies of the original, with the prefix
    "turbulenceModel:", e.g.:

    \verbatim
      turbulenceModel:R
    \endverbatim

    Operands:
    \table
      Operand           | Type             | Location
      input             | -                | -
      output file       | -                | -
      output field      | vol\<Type\>Field | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

    where \c \<Type\>=Scalar/SymmTensor.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    turbulenceFields1
    {
        // Mandatory entries (unmodifiable)
        type        turbulenceFields;
        libs        (fieldFunctionObjects);

        // Mandatory entries (runtime modifiable)
        // Either field or fields entries
        fields      (R devRhoReff);
        field       R;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                           | Type | Req'd | Dflt
      type         | Type name: turbulenceFields           | word |  yes  | -
      libs         | Library name: fieldFunctionObjects    | word |  yes  | -
      fields       | Names of fields to store (see below)  | wordList | yes | -
      field        | Name of a field to store (see below)  | word | yes | -
    \endtable

    where \c fields can include:
    \verbatim
      k           | turbulent kinetic energy
      epsilon     | turbulent kinetic energy dissipation rate
      omega       | specific dissipation rate
      nuTilda     | modified turbulent viscosity
      nut         | turbulent viscosity (incompressible)
      nuEff       | effective turbulent viscosity (incompressible)
      mut         | turbulent viscosity (compressible)
      muEff       | effective turbulent viscosity (compressible)
      alphat      | turbulence thermal diffusivity (compressible)
      alphaEff    | effective turbulence thermal diffusivity (compressible)
      R           | Reynolds stress tensor
      devReff     | deviatoric part of the effective Reynolds stress
      devRhoReff  | divergence of the Reynolds stress
      L           | integral-length scale / mixing-length scale
      I           | turbulence intensity
    \endverbatim

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        <solver> -postProcess -func turbulenceFields
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::turbulenceFields

SourceFiles
    turbulenceFields.C
    turbulenceFieldsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_turbulenceFields_H
#define functionObjects_turbulenceFields_H

#include "fvMeshFunctionObject.H"
#include "HashSet.H"
#include "Enum.H"
#include "volFieldsFwd.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                      Class turbulenceFields Declaration
\*---------------------------------------------------------------------------*/

class turbulenceFields
:
    public fvMeshFunctionObject
{
public:

    // Public Enumerations

        //- Options for the turbulence fields (compressible)
        enum compressibleField
        {
            cfK,            //!< "Turbulent kinetic energy"
            cfEpsilon,      //!< "Turbulent kinetic energy dissipation rate"
            cfOmega,        //!< "Specific dissipation rate"
            cfNuTilda,      //!< "Modified turbulent viscosity"
            cfMut,          //!< "Turbulent viscosity"
            cfMuEff,        //!< "Effective turbulent viscosity"
            cfAlphat,       //!< "Turbulence thermal diffusivity"
            cfAlphaEff,     //!< "Effective turbulence thermal diffusivity"
            cfR,            //!< "Reynolds stress tensor"
            cfDevRhoReff,   //!< "Divergence of the Reynolds stress"
            cfL,            //!< "Integral-length/Mixing-length scale"
            cfI             //!< "Turbulence intensity"
        };

        //- Names for compressibleField turbulence fields
        static const Enum<compressibleField> compressibleFieldNames_;

        //- Options for the turbulence fields (incompressible)
        enum incompressibleField
        {
            ifK,            //!< "Turbulent kinetic energy"
            ifEpsilon,      //!< "Turbulent kinetic energy dissipation rate"
            ifOmega,        //!< "Specific dissipation rate"
            ifNuTilda,      //!< "Modified turbulent viscosity"
            ifNut,          //!< "Turbulent viscosity"
            ifNuEff,        //!< "Effective turbulent viscosity"
            ifR,            //!< "Reynolds stress tensor"
            ifDevReff,   //!< "Deviatoric part of the effective Reynolds stress"
            ifL,            //!< "Integral-length/Mixing-length scale"
            ifI             //!< "Turbulence intensity"
        };

        //- Names for incompressibleField turbulence fields
        static const Enum<incompressibleField> incompressibleFieldNames_;

        //- Turbulence closure model name
        static const word modelName_;

protected:

    // Protected Data

        //- Fields to load
        wordHashSet fieldSet_;


    // Protected Member Functions

        //- Return true if compressible turbulence model is identified
        bool compressible();

        //- Process the turbulence field
        template<class Type>
        void processField
        (
            const word& fieldName,
            const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvalue
        );

        //- Return omega calculated from k and epsilon
        template<class Model>
        tmp<volScalarField> omega(const Model& model) const;

        //- Return nuTilda calculated from k and omega
        template<class Model>
        tmp<volScalarField> nuTilda(const Model& model) const;

        //- Return integral length scale, L, calculated from k and epsilon
        template<class Model>
        tmp<volScalarField> L(const Model& model) const;

        //- Return turbulence intensity, I, calculated from k and U
        template<class Model>
        tmp<volScalarField> I(const Model& model) const;


public:

    //- Runtime type information
    TypeName("turbulenceFields");


    // Constructors

        //- Construct from Time and dictionary
        turbulenceFields
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        turbulenceFields(const turbulenceFields&) = delete;

        //- No copy assignment
        void operator=(const turbulenceFields&) = delete;


    //- Destructor
    virtual ~turbulenceFields() = default;


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary&);

        //- Calculate turbulence fields
        virtual bool execute();

        //- Do nothing.
        //  The turbulence fields are registered and written automatically
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "turbulenceFieldsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
