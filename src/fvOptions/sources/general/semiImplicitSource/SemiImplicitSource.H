/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::fv::SemiImplicitSource

Group
    grpFvOptionsSources

Description
    Applies semi-implicit source within a specified region for \c Type,
    where \c \<Type\>=Scalar/Vector/SphericalTensor/SymmTensor/Tensor.
    The injection rate coefficients are specified
    as pairs of Su-Sp coefficients, i.e.:

    \f[
        S(x) = S_u + S_p x
    \f]

    where
    \vartable
      S(x)    | net source for field 'x'
      S_u     | explicit source contribution
      S_p     | linearised implicit contribution
    \endvartable

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    <Type>SemiImplicitSource1
    {
        // Mandatory entries (unmodifiable)
        type                <Type>SemiImplicitSource;

        // Mandatory entries (runtime modifiable)
        volumeMode          <volumeModeType>;
        injectionRateSuSp
        {
            k           (30.7 0);
            epsilon     (1.5  0);

            // The injectionRate can also be specified as a Function1
            // by having dictionaries for the field entries instead
            k
            {
                // Time-ramp from 0 to 30.7 at time 5
                Su table
                (
                    (0 0.0)
                    (5 30.7)
                );
                Sp      0.0;
            }
            epsilon
            {
                Su      1.5;
                Sp      0.0;
            }
        }

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                          | Type | Reqd | Dflt
      type       | Type name: \<Type\>SemiImplicitSource <!--
                -->                                     | word | yes  | -
      volumeMode | Volume mode type                     | word | yes  | -
      injectionRateSuSp | Injection rate settings  | dictionary | yes | -
    \endtable

    The inherited entries are elaborated in:
     - \link fvOption.H \endlink
     - \link cellSetOption.H \endlink

    Options for the \c volumeMode entry:
    \verbatim
      absolute  | Values are given as \<quantity\>
      specific  | Values are given as \<quantity\>/m3
    \endverbatim

See also
  - Foam::fvOption

SourceFiles
    SemiImplicitSource.C

\*---------------------------------------------------------------------------*/

#ifndef SemiImplicitSource_H
#define SemiImplicitSource_H

#include "Tuple2.H"
#include "cellSetOption.H"
#include "Enum.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

// Forward declarations
template<class Type> class SemiImplicitSource;

template<class Type>
Ostream& operator<<
(
    Ostream&,
    const SemiImplicitSource<Type>&
);


/*---------------------------------------------------------------------------*\
                     Class SemiImplicitSource Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class SemiImplicitSource
:
    public fv::cellSetOption
{
public:

    // Public Enumerations

        //- Options for the volume mode type
        enum volumeModeType
        {
            vmAbsolute,
            vmSpecific
        };

        //- Names for volumeModeType
        static const Enum<volumeModeType> volumeModeTypeNames_;


protected:

    // Protected Data

        //- Volume mode
        volumeModeType volumeMode_;

        //- Volume normalisation
        scalar VDash_;

        //- Source field values
        PtrList<Function1<Type>> Su_;
        PtrList<Function1<scalar>> Sp_;


    // Protected Functions

        //- Set the local field data
        void setFieldData(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("SemiImplicitSource");


    // Constructors

        //- Construct from components
        SemiImplicitSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Member Functions

        // Access

            //- Return const access to the volume mode
            inline const volumeModeType& volumeMode() const;


        // Edit

            //- Return access to the volume mode
            inline volumeModeType& volumeMode();


        // Evaluation

            //- Add explicit contribution to equation
            virtual void addSup
            (
                fvMatrix<Type>& eqn,
                const label fieldi
            );

            //- Add explicit contribution to compressible equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<Type>& eqn,
                const label fieldi
            );


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SemiImplicitSource.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SemiImplicitSourceI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
