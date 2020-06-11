/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::turbulentDigitalFilterFvPatchVectorField

Group
    grpInletBoundaryConditions

Description
    Digital-filter based boundary condition for velocity, i.e. \c U, to generate
    synthetic turbulence-alike time-series for LES and DES turbulent flow
    computations from input turbulence statistics.

    References:
    \verbatim
        Digital-filter method-based generator (DFM) (tag:KSJ):
            Klein, M., Sadiki, A., & Janicka, J. (2003).
            A digital filter based generation of inflow data for spatially
            developing direct numerical or large eddy simulations.
            Journal of computational Physics, 186(2), 652-665.
            DOI:10.1016/S0021-9991(03)00090-1

        Forward-stepwise method-based generator (FSM) (tag:XC)
            Xie, Z. T., & Castro, I. P. (2008).
            Efficient generation of inflow conditions for
            large eddy simulation of street-scale flows.
            Flow, turbulence and combustion, 81(3), 449-470.
            DOI:10.1007/s10494-008-9151-5

        Mass-inflow rate correction (tag:KCX):
            Kim, Y., Castro, I. P., & Xie, Z. T. (2013).
            Divergence-free turbulence inflow conditions for
            large-eddy simulations with incompressible flow solvers.
            Computers & Fluids, 84, 56-68.
            DOI:10.1016/j.compfluid.2013.06.001
    \endverbatim

    In \c DFM or \c FSM, a random number set (mostly white noise), and a group
    of target statistics (mostly mean flow, Reynolds stress tensor profiles and
    length-scale sets) are merged into a new number set (stochastic time-series,
    yet consisting of the statistics) by a chain of mathematical operations
    whose characteristics are designated by the target statistics, so that the
    realised statistics of the new sets could match the target.

    \verbatim
    Random number sets ---->-|
                             |
                         DFM or FSM ---> New stochastic time-series consisting
                             |           turbulence statistics
    Turbulence statistics ->-|
    \endverbatim

    The main difference between \c DFM and \c FSM is that \c FSM replaces
    the expensive-to-run streamwise convolution summation in \c DFM by a simpler
    and an almost-equivalent-in-effect numerical procedure in order to reduce
    computational costs. Accordingly, \c FSM potentially brings computational
    resource advantages for computations involving relatively large streamwise
    length-scale sets and small time-steps.

    Synthetic turbulence is generated on a virtual rectangular structured-mesh
    plane, which is parallel to the chosen patch, and is mapped onto this patch
    by the selected mapping method.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries (unmodifiable)
        type                turbulentDigitalFilterInlet;
        n                   (<nHeight> <nWidth>);
        L                   (<L1> <L2> ... <L9>);
        R                   uniform (<Rxx> <Rxy> <Rxz> <Ryy> <Ryz> <Rzz>);
        UMean               uniform (1 0 0);
        Ubulk               10.0;

        // Optional entries (unmodifiable)
        fsm                 false;
        Gaussian            true;     // always false for FSM
        fixSeed             true;
        continuous          false;
        correctFlowRate     true;
        mapMethod           nearestCell;
        perturb             1e-5;
        C1                  -1.5707;  //-0.5*PI;
        C1FSM               -0.7854   //-0.25*PI;
        C2FSM               -1.5707;  //-0.5*PI;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                            | Type  | Req'd | Dflt
      type     | Type name: turbulentDigitalFilterInlet | word  | yes   | -
      n        | Number of cells on turbulence generation plane <!--
                                           --> | tuple of labels | yes   | -
      L        | Integral length-scale set <!--
              --> (Lxu Lxv Lxw Lyu Lyv Lyw Lzu Lzv Lzw) [m] | tensor | yes | -
      R        | Reynolds stress tensor set <!--
              --> (xx xy xz yy yz zz) [m2/s2]      | symmTensorField | yes | -
      UMean    | Mean velocity profile [m/s]       | vectorField     | yes | -
      Ubulk    | Characteristic patch-normal bulk flow speed  [m/s] <!--
                                                        --> | scalar | yes | -
      fsm      | Flag to turn on the forward-stepwise method | bool | no | false
      Gaussian | Autocorrelation function form         | bool | no | true
      fixSeed  | Flag to fix random-number generator seed to 1234 <!--
             --> or generate a new seed based on clock-time per simulation <!--
                                                     --> | bool | no | true
      continuous | Flag to write random-number sets at output time, <!--
             --> and to read them on restart. Otherwise, generate <!--
             --> new random-number sets of restart       | bool | no  false
      correctFlowRate | Flag to correct mass-inflow rate on turbulence <!--
                    --> plane in (only) streamwise direction | bool | no | true
      mapMethod  | Interpolation-to-patch method      | word | no | nearestCell
      perturb    | Point perturbation for planarInterpolation mapMethod <!--
                                                  --> | scalar | no | 1e-5
      C1         | Model constant shaping autocorrelation function <!--
                                     --> (KSJ:Eq. 14) | scalar | no | -0.5*PI
      C1FSM | Model coefficient in FSM (XC:Eq. 14)    | scalar | no | -0.25*PI
      C2FSM | Model coefficient in FSM (XC:Eq. 14)    | scalar | no | -0.5*PI
    \endtable

    The inherited entries are elaborated in:
      - \link fixedValueFvPatchFields.H \endlink

    Options for the \c fsm entry:
    \verbatim
      false | Method due to (KSJ)
      true  | Method due to (XC)
    \endverbatim

    Options for the \c Gaussian entry:
    \verbatim
      true  | Gaussian function
      false | Exponential function (only option for FSM)
    \endverbatim

    Options for the \c mapMethod entry:
    \verbatim
      nearestCell         | One-to-one direct map, no interpolation
      planarInterpolation | Bilinear interpolation
    \endverbatim

    Patch-profile input is available for two entries:
    \verbatim
      R     | Reynolds stress tensor
      UMean | Mean velocity
    \endverbatim
    where the input profiles and profile coordinates are located in:
    \verbatim
      Coordinates | $FOAM_CASE/constant/boundaryData/\<patchName\>/points
      R/UMean     | $FOAM_CASE/constant/boundaryData/\<patchName\>/0/\{R/UMean\}
    \endverbatim

    \c points file contains a list of three-dimensional coordinates, and
    profile data files provide a value corresponding to each coordinate.

Note
    - \c mapMethod=planarInterpolation option needs point coordinates that can
    form a plane.
    - \c adjustTimeStep=true option is currently not fully supported.
    - In order to obtain Reynolds stress tensor information, experiments, RANS
    simulations or engineering relations can be used.
    - \c continuous=true means deterministic-statistical consistent restart
    (relatively more expensive), and \c continuous=false means deterministic
    discontinuity in synthetic turbulence time-series by keeping statistical
    consistency (relatively cheaper).
    - For \c L, the first three entries should always correspond to the
    length scales in association with the convective (streamwise) mean flow
    direction.
    - Streamwise integral length scales are converted to integral time scales
    by using Taylor's frozen turbulence hypothesis, and \c Ubulk.

See also
    - turbulentDFSEMInletFvPatchVectorField.C

SourceFiles
    turbulentDigitalFilterInletFvPatchVectorField.C
    turbulentDigitalFilterInletFvPatchVectorFieldTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentDigitalFilterInletFvPatchVectorField_H
#define turbulentDigitalFilterInletFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Random.H"
#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class pointToPointPlanarInterpolation;

/*---------------------------------------------------------------------------*\
        Class turbulentDigitalFilterInletFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class turbulentDigitalFilterInletFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Bilinear interpolation (for 'mapMethod=planarInterpolation')
        mutable autoPtr<pointToPointPlanarInterpolation> mapperPtr_;

        //- Flag to enable the forward-stepwise method
        const bool fsm_;

        //- Flag to select correlation function form: Gaussian or exponential
        const bool Gaussian_;

        //- Flag to fix the random-number generator seed to 1234 or
        //- generate a new seed based on clock-time per simulation
        const bool fixSeed_;

        //- Flag to write random-number sets at output time, and to read them
        //- on restart. Otherwise, generate new random-number sets on restart
        const bool continuous_;

        //- Flag to correct mass flow rate on turbulence plane
        const bool correctFlowRate_;

        //- Internal flag to read R from data files
        bool interpR_;

        //- Internal flag to read UMean from data files
        bool interpUMean_;

        //- Method for interpolation between a patch and turbulence plane
        const word mapMethod_;

        //- Current time index
        label curTimeIndex_;

        //- Characteristic patch-normal bulk flow speed [m/s]
        const scalar Ubulk_;

        //- Model constant shaping autocorrelation function (KSJ:Eq. 14)
        const scalar C1_;

        //- Fraction of perturbation (fraction of bounding box) to add
        const scalar perturb_;

        //- First time-step mass/volumetric flow rate
        scalar flowRate_;

        //- Random number generator
        Random rndGen_;

        //- Number of cells on turbulence plane (<nHeight> <nWidth>) [-]
        const Tuple2<label, label> n_;

        //- Uniform mesh size on turbulence plane (reversed) [1/m]
        Vector2D<scalar> invDelta_;

        //- Index pairs between patch and turbulence plane
        //- for the nearest cell mapping
        const List<Pair<label>> indexPairs_;

        //- Reynolds stress tensor profile field in global coordinates [m2/s2]
        const symmTensorField R_;

        //- Lund-Wu-Squires transformed R field (Cholesky decomp.) [m/s]
        //- Mapped onto actual mesh patch rather than turbulence plane
        //  (KSJ:Eq. 5)
        const symmTensorField Lund_;

        //- Mean inlet velocity profile field in global coordinates [m/s]
        vectorField UMean_;

        //- Integral length-scale set per turbulence plane section in local
        //- coordinates (e1u, e1v, e1w, e2u, e2v, e2w, e3u, e3v, e3w) [m]
        //  Backup of L_ for restart purposes
        const tensor Lbak_;

        //- Integral length-scale set in mesh units [cell]
        const tensor L_;

        //- One of the two model coefficients in FSM
        //  (XC:The argument of the first exp func in Eq. 14)
        const scalar C1FSM_;

        //- One of the two model coefficients in FSM
        //  (XC:The argument of the second exp func in Eq. 14)
        const scalar C2FSM_;

        //- One of the two exponential functions in FSM
        //  (XC:The first exponential function in Eq. 14)
        const List<scalar> coeffs1FSM_;

        //- One of the two exponential functions in FSM
        //  (XC:The first exponential function in Eq. 14)
        const List<scalar> coeffs2FSM_;

        //- Number of cells in random-number box [cell]
        //- Random-number sets within box are filtered with filterCoeffs_
        const List<label> szBox_;

        //- Convenience factors for two-dimensional random-number box [-]
        const List<label> boxFactors2D_;

        //- Convenience factors for three-dimensional random-number box [-]
        const List<label> boxFactors3D_;

        //- Index to the first elem of last plane of random-number box [-]
        const List<label> iNextToLastPlane_;

        //- Random-number sets distributed over three-dimensional box (u, v, w)
        List<List<scalar>> box_;

        //- Filter coefficients corresponding to L [-]
        const List<List<scalar>> filterCoeffs_;

        //- Filter-applied random-number sets [m/s], i.e. turbulence plane
        List<List<scalar>> filteredBox_;

        //- Filter-applied previous-time-step velocity field [m/s] used in FSM
        vectorField U0_;


    // Private Member Functions

        //- Return a reference to the patch mapper object
        const pointToPointPlanarInterpolation& patchMapper() const;

        //- Helper function to interpolate values from the boundary data or
        //- read from dictionary
        template<class Type>
        tmp<Field<Type>> interpolateOrRead
        (
            const word& fieldName,
            const dictionary& dict,
            bool& interpolateField
        ) const;

        //- Helper function to interpolate values from the boundary data
        template<class Type>
        tmp<Field<Type>> interpolateBoundaryData
        (
            const word& fieldName
        ) const;

        //- Generate random-number sets obeying the standard normal distribution
        template<class Form, class Type>
        Form randomSet(const label len);

        //- Compute nearest cell index-pairs between turbulence plane and patch
        List<Pair<label>> indexPairs();

        //- Check R on patch for mathematical domain errors
        void checkR() const;

        //- Compute Lund-Wu-Squires transformation
        symmTensorField calcLund() const;

        //- Compute the first time-step mass/
        //- volumetric flow rate based on UMean
        scalar calcFlowRate() const;

        //- Convert integral length scales in meters
        //- to turbulence plane cell-size units
        tensor meterToCell(const tensor& L) const;

        //- Resource allocation functions for the convenience factors
        List<label> initBox() const;
        List<label> initFactors2D() const;
        List<label> initFactors3D() const;
        List<label> initPlaneFactors() const;

        //- Compute various convenience factors for random-number box
        List<List<scalar>> fillBox();

        //- Compute filter coeffs once per simulation
        List<List<scalar>> calcFilterCoeffs() const;

        //- Discard current time-step random-box plane (closest to patch) by
        //- shifting from the back to the front, and add new plane to the back
        void shiftRefill();

        //- Map two-point correlated random-number sets
        //- on patch based on chosen mapping method
        void mapFilteredBox(vectorField& U);

        //- Embed one-point correlations, i.e. R, on patch
        void onePointCorrs(vectorField& U) const;

        //- Embed two-point correlations, i.e. L, on box
        //  Three-dimensional "valid"-type separable
        //  convolution summation algorithm
        //  (Based on Song Ho Ahn's two-dimensional "full"-type convolution)
        void twoPointCorrs();

        //- Compute coeffs1FSM_ once per simulation
        List<scalar> calcCoeffs1FSM() const;

        //- Compute coeffs2FSM_ once per simulation
        List<scalar> calcCoeffs2FSM() const;


public:

   //- Runtime type information
   TypeName("turbulentDigitalFilterInlet");


    // Constructors

        //- Construct from patch and internal field
        turbulentDigitalFilterInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentDigitalFilterInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- turbulentDigitalFilterInletFvPatchVectorField onto a new patch
        turbulentDigitalFilterInletFvPatchVectorField
        (
            const turbulentDigitalFilterInletFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentDigitalFilterInletFvPatchVectorField
        (
            const turbulentDigitalFilterInletFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new turbulentDigitalFilterInletFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        turbulentDigitalFilterInletFvPatchVectorField
        (
            const turbulentDigitalFilterInletFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new turbulentDigitalFilterInletFvPatchVectorField(*this, iF)
            );
        }


    //- Destructor
    virtual ~turbulentDigitalFilterInletFvPatchVectorField() = default;


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper& m);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField& ptf,
                const labelList& addr
            );


       //- Write
       virtual void write(Ostream&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "turbulentDigitalFilterInletFvPatchVectorFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
