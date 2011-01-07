/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2011 OpenCFD Ltd.
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
    Foam::ReactingCloud

Description
    Templated base class for reacting cloud

    - Adds to thermodynamic cloud
      - Variable composition (single phase)
      - Phase change

SourceFiles
    ReactingCloudI.H
    ReactingCloud.C

\*---------------------------------------------------------------------------*/

#ifndef ReactingCloud_H
#define ReactingCloud_H

#include "ThermoCloud.H"
#include "reactingCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class CompositionModel;

template<class CloudType>
class PhaseChangeModel;

/*---------------------------------------------------------------------------*\
                      Class ReactingCloud Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class ReactingCloud
:
    public ThermoCloud<ParcelType>,
    public reactingCloud
{

private:

    // Private data

        //- Cloud copy pointer
        autoPtr<ReactingCloud<ParcelType> > cloudCopyPtr_;


    // Private member functions

        //- Disallow default bitwise copy construct
        ReactingCloud(const ReactingCloud&);

        //- Disallow default bitwise assignment
        void operator=(const ReactingCloud&);


protected:

    // Protected data

        //- Parcel constant properties
        typename ParcelType::constantProperties constProps_;


        // References to the cloud sub-models

            //- Reacting composition model
            autoPtr<CompositionModel<ReactingCloud<ParcelType> > >
                compositionModel_;

            //- Reacting phase change model
            autoPtr<PhaseChangeModel<ReactingCloud<ParcelType> > >
                phaseChangeModel_;


        // Sources

            //- Mass transfer fields - one per carrier phase specie
            PtrList<DimensionedField<scalar, volMesh> > rhoTrans_;


        // Check

            //- Total mass transferred to continuous phase via phase change
            scalar dMassPhaseChange_;


    // Protected Member Functions

        // New parcel helper functions

            //- Check that size of a composition field is valid
            void checkSuppliedComposition
            (
                const scalarField& YSupplied,
                const scalarField& Y,
                const word& YName
            );


        // Cloud evolution functions

            //- Reset state of cloud
            void cloudReset(ReactingCloud<ParcelType>& c);


public:

    // Constructors

        //- Construct given carrier gas fields
        ReactingCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            const SLGThermo& thermo,
            bool readFields = true
        );

        //- Copy constructor with new name
        ReactingCloud(ReactingCloud<ParcelType>& c, const word& name);

        //- Copy constructor with new name - creates bare cloud
        ReactingCloud
        (
            const fvMesh& mesh,
            const word& name,
            const ReactingCloud<ParcelType>& c
        );

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<ParcelType> > clone(const word& name)
        {
            return autoPtr<Cloud<ParcelType> >
            (
                new ReactingCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<ParcelType> > cloneBare(const word& name) const
        {
            return autoPtr<Cloud<ParcelType> >
            (
                new ReactingCloud(this->mesh(), name, *this)
            );
        }


    //- Destructor
    virtual ~ReactingCloud();


    //- Type of parcel the cloud was instantiated for
    typedef ParcelType parcelType;


    // Member Functions

        // Access

            //- Return a reference to the cloud copy
            inline const ReactingCloud& cloudCopy() const;

            //- Return the constant properties
            inline const typename ParcelType::constantProperties&
                constProps() const;


            // Sub-models

                //- Return reference to reacting composition model
                inline const CompositionModel<ReactingCloud<ParcelType> >&
                    composition() const;

                //- Return reference to reacting phase change model
                inline const PhaseChangeModel<ReactingCloud<ParcelType> >&
                    phaseChange() const;


            // Sources

                //- Mass

                    //- Return reference to mass source for field i
                    inline DimensionedField<scalar, volMesh>&
                        rhoTrans(const label i);

                    //- Return const access to mass source fields
                    inline const PtrList<DimensionedField<scalar, volMesh> >&
                        rhoTrans() const;

                    //- Return reference to mass source fields
                    inline PtrList<DimensionedField<scalar, volMesh> >&
                        rhoTrans();

                    //- Return mass source term for specie i - specie eqn
                    inline tmp<fvScalarMatrix> SYi
                    (
                        const label i,
                        volScalarField& Yi
                    ) const;

                    //- Return tmp mass source for field i - fully explicit
                    inline tmp<DimensionedField<scalar, volMesh> >
                        Srho(const label i) const;

                    //- Return tmp total mass source for carrier phase
                    //  - fully explicit
                    inline tmp<DimensionedField<scalar, volMesh> > Srho() const;

                    //- Return total mass source term [kg/m3/s]
                    inline tmp<fvScalarMatrix> Srho(volScalarField& rho) const;


        // Check

            //- Add to cumulative phase change mass transfer
            void addToMassPhaseChange(const scalar dMass);


        // Cloud evolution functions

            //- Check parcel properties
            void checkParcelProperties
            (
                ParcelType& parcel,
                const scalar lagrangianDt,
                const bool fullyDescribed
            );

            //- Store the current cloud state
            void storeState();

            //- Reset the current cloud to the previously stored state
            void restoreState();

            //- Reset the cloud source terms
            void resetSourceTerms();

            //- Apply relaxation to (steady state) cloud sources
            void relaxSources(const ReactingCloud<ParcelType>& cloudOldTime);

            //- Evolve the cloud
            void evolve();

            //- Print cloud information
            void info() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ReactingCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ReactingCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
