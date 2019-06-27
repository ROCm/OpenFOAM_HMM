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
    Foam::turbulentDigitalFilterFvPatchVectorField

Group
    grpInletBoundaryConditions

Description
    Velocity boundary condition generating synthetic turbulence-alike
    time-series for LES and DES turbulent flow computations.

    To this end, two synthetic turbulence generators can be chosen:

    - Digital-filter method-based generator (DFM)

    \verbatim
    Klein, M., Sadiki, A., and Janicka, J.
        A digital filter based generation of inflow data for spatially
        developing direct numerical or large eddy simulations,
        Journal of Computational Physics (2003) 186(2):652-665.
        doi:10.1016/S0021-9991(03)00090-1
    \endverbatim

    - Forward-stepwise method-based generator (FSM)

    \verbatim
    Xie, Z.-T., and Castro, I.
        Efficient generation of inflow conditions for large eddy simulation of
        street-scale flows, Flow, Turbulence and Combustion (2008) 81(3):449-470
        doi:10.1007/s10494-008-9151-5
    \endverbatim

    In DFM or FSM, a random number set (mostly white noise), and a group
    of target statistics (mostly mean flow, Reynolds stress tensor profiles and
    length-scale sets) are fused into a new number set (stochastic time-series,
    yet consisting of the statistics) by a chain of mathematical operations
    whose characteristics are designated by the target statistics, so that the
    realised statistics of the new sets could match the target.

    Random number sets ---->-|
                             |
                         DFM or FSM ---> New stochastic time-series consisting
                             |           turbulence statistics
    Turbulence statistics ->-|

    The main difference between DFM and FSM is that the latter replaces the
    streamwise convolution summation in DFM by a simpler and a quantitatively
    justified equivalent procedure in order to reduce computational costs.
    Accordingly, the latter potentially brings resource advantages for
    computations involving relatively large length-scale sets and small
    time-steps.

    Synthetic turbulence is produced on a virtual rectangular structured-mesh
    turbulence plane, which is parallel to the actual patch, and is mapped onto
    the chosen patch by the selected mapping method.

Usage
    \table
     Property    | Description                     | Required    | Default value
     planeDivisions   | Number of nodes on turbulence plane (e2, e3) [-] | yes |
     L        | Integral length-scale set (9-comp):{e1,e2,e3}{u,v,w} [m] | yes |
     R         | Reynolds stress tensor set (xx xy xz yy yz zz) [m2/s2]  | yes |
     patchNormalSpeed    | Characteristic mean flow speed  [m/s] | yes |
     isGaussian  | Autocorrelation function form                 | no | true
     isFixedSeed | Flag to identify random-number seed is fixed  | no | true
     isContinuous   | Flag for random-number restart behaviour   | no | false
     isCorrectedFlowRate | Flag for mass flow rate correction    | no | true
     interpolateR     | Placeholder flag: interpolate R field    | no | false
     interpolateUMean | Placeholder flag: interpolate UMean field   | no | false
     isInsideMesh | Placeholder flag: TP is inside mesh or on patch | no | false
     isTaylorHypot | Placeholder flag: Taylor's hypothesis is on | no | true
     mapMethod   | Method to map reference values             | no | nearestCell
     threshold   | Threshold to avoid unintentional 'tiny' input | no | 1e-8
     modelConst  | Model constant (Klein et al., Eq. 14)         | no | -0.5*PI
     perturb     | Point perturbation for interpolation          | no | 1e-5
     const1FSM | A model coefficient in FSM (Xie-Castro, Eq. 14) | no | -0.25*PI
     const2FSM | A model coefficient in FSM (Xie-Castro, Eq. 14) | no | -0.5*PI
    \endtable

    Minimal example of the boundary condition specification with commented
    options:
    \verbatim
    <patchName>
    {
        type                turbulentDigitalFilter;
        variant             digitalFilter;          // reducedDigitalFilter;
        planeDivisions      (<planeDivisionsHeight> <planeDivisionsWidth>);
        L               (<Lxu> <Lxv> <Lxw> <Lyu> <Lyv> <Lyw> <Lzu> <Lzv> <Lzw>);
        R                   (<Rxx> <Rxy> <Rxz> <Ryy> <Ryz> <Rzz>);
        patchNormalSpeed    <characteristic flow speed>;
        value               uniform (0 0 0);        // mandatory placeholder

        // Optional entries with default input
        isGaussian          true;           // false    // always false for FSM
        isFixedSeed         true;           // false
        isContinuous        false;          // true
        isCorrectedFlowRate true;           // false
        interpolateR        false;          // placeholder
        interpolateUMean    false;          // placeholder
        isInsideMesh        false;          // placeholder
        isTaylorHypot       true;           // placeholder
        mapMethod           nearestCell;    // planarInterpolation
        threshold           1e-8;
        modelConst          -1.5707;        //-0.5*PI;
        perturb             1e-5;

        // Optional entries for only FSM with default input
        const1FSM           -0.7854         //-0.25*PI;
        const2FSM           -1.5707;        //-0.5*PI;
    }
    \endverbatim

    Among the dictionary entries, two entries can be input as patch profiles:

    \verbatim
     - Reynolds stress tensor, R
     - Mean velocity, UMean
    \endverbatim

    Profile data and corresponding coordinates are then input in the following
    directories:

    \verbatim
     - $FOAM_CASE/constant/boundaryData/\<patchName\>/points
     - $FOAM_CASE/constant/boundaryData/\<patchName\>/0/\{R|UMean\}
    \endverbatim

    The profile data and corresponding coordinates take the same form used by
    the \c timeVaryingMappedFixedValue and \c turbulentDFSEMInlet boundary
    conditions, consisting of a \c points file containing a list of 3-D
    coordinates, and profile data files providing a value per coordinate.

    Reynolds stress tensor and mean velocity input are in the global coordinate
    system whereas integral length scale set input is in the local patch
    coordinate system.

    It is assumed that the patch normal direction is \c e1, and the remaining
    patch plane directions are \c e2 and \c e3 following the right-handed
    coordinate system, which should be taken into consideration when the
    integral length scale set is input. The first three integral scale entries,
    i.e. L_e1u, Le1v, Le2w, should always correspond to the length scales
    that are in association with the convective mean flow direction.

Note
    - \c mapMethod \c planarInterpolation option requires point coordinates
    which can form a plane, thus input point coordinates varying only in a
    single direction will trigger error.
    - \c adjustTimeStep = true option is not fully supported at the moment.

SeeAlso
    turbulentDFSEMInletFvPatchVectorField.C

SourceFiles
    turbulentDigitalFilterInletFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentDigitalFilterInletFvPatchVectorField_H
#define turbulentDigitalFilterInletFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Random.H"
#include <functional>
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

    // Private Enumerations

        //- Options for the synthetic turbulence generator variant
        enum variantType : uint8_t
        {
            DIGITAL_FILTER = 1,    //!< Digital-filter method (Klein et al.)
            FORWARD_STEPWISE = 2,  //!< Forward-stepwise method (Xie-Castro)
        };

        //- Names for variant types
        static const Enum<variantType> variantNames;


    // Private Data

        //- 2D interpolation (for 'planarInterpolation' mapMethod)
        mutable autoPtr<pointToPointPlanarInterpolation> mapperPtr_;

        //- Selected option for the synthetic turbulence generator variant
        const enum variantType variant_;

        //- Flag: correlation function form is Gaussian or Exponential
        //  for variantType::DIGITAL_FILTER, default=Gaussian
        //  for variantType::FORWARD_STEPWISE, default=Exponential (only option)
        const bool isGaussian_;

        //- Flag: random-number generator seed is fixed (default=true) or
        //- generated pseudo-randomly based on clock-time per simulation
        const bool isFixedSeed_;

        //- Flag: write non-manipulated random-number sets at output time, and
        //- to read them on restart. Otherwise, generate new random-number sets
        //- on restart. (default=false)
        //  true: deterministic & statistically consistent, more expensive
        //  false: deterministic discontinuity & statistically consistent, cheaper
        const bool isContinuous_;

        //- Flag: mass flow rate is corrected on turbulence plane (default=true)
        const bool isCorrectedFlowRate_;

        //- Flag: interpolate R field (default=false)
        bool interpolateR_;

        //- Flag: interpolate UMean field (default=false)
        bool interpolateUMean_;

        //- Flag: turbulence plane is inside mesh or on a patch (default=false)
        //  Currently, true option is not available.
        const bool isInsideMesh_;

        //- Flag: convert streamwise (x) length scales to time scales by
        //- Taylor's 'frozen turbulence' hypothesis (default=true)
        //  Currently, false option is not available.
        const bool isTaylorHypot_;

        //- Method for interpolation between a patch and turbulence plane
        //- (default=nearestCell)
        //  Options:
        //   - nearestCell: one-to-one direct map, no interpolation
        //   - planarInterpolation: bilinear interpolation
        const word mapMethod_;

        //- Current time index
        label curTimeIndex_;

        //- Threshold to avoid unintentional 'tiny' input scalars (default=1e-8)
        const scalar tiny_;

        //- Characteristic (e.g. bulk) mean speed of flow in the patch normal
        //- direction [m/s]
        const scalar patchNormalSpeed_;

        //- Model constant shaping autocorr function [-] (default='-0.5*pi')
        //  (Klein et al., 2003, Eq. 14)
        const scalar modelConst_;

        //- Fraction of perturbation (fraction of bounding box) to add (for
        //- 'planarInterpolation' mapMethod)
        const scalar perturb_;

        //- Initial (first time-step) mass/vol flow rate [m^3/s]
        scalar initialFlowRate_;

        //- Random number generator
        Random rndGen_;

        //- Number of nodes on turbulence plane (e2 e3) [-]
        const Tuple2<label, label> planeDivisions_;

        //- Turbulence plane mesh size (reversed) (e2 e3) [1/m]
        Vector2D<scalar> invDelta_;

        //- Nearest cell mapping: Index pairs between patch and turbulence plane
        const List<Pair<label>> indexPairs_;

        //- Reynolds stress tensor profile (xx xy xz yy yz zz) in global
        //- coordinates [m^2/s^2]
        symmTensorField R_;

        //- Lund-Wu-Squires transformation (Cholesky decomp.) [m/s]
        //- Mapped onto actual mesh patch rather than turbulence plane
        //  (Klein et al., 2003, Eq. 5)
        symmTensorField LundWuSquires_;

        //- Mean inlet velocity profile in global coordinates [m/s]
        vectorField UMean_;

        //- Integral length-scale set per turbulence plane section in local
        //- coordinates (e1u, e1v, e1w, e2u, e2v, e2w, e3u, e3v, e3w) [m]
        //  First three entries should always correspond to the length scales
        //  in association with the convective mean flow direction
        //  Backup of L_ for restart purposes
        const tensor Lbak_;

        //- Integral length-scale set in mesh units [node]
        const tensor L_;

        //- One of the two model coefficients in FSM
        //  (Xie-Castro, 2008, the argument of the first exp func in Eq. 14)
        const scalar const1FSM_;

        //- One of the two model coefficients in FSM
        //  (Xie-Castro, 2008, the argument of the second exp func in Eq. 14)
        const scalar const2FSM_;

        //- One of the two exponential functions in FSM
        //  (Xie-Castro, 2008, the first exponential function in Eq. 14)
        const List<scalar> constList1FSM_;

        //- One of the two exponential functions in FSM
        //  (Xie-Castro, 2008, the first exponential function in Eq. 14)
        const List<scalar> constList2FSM_;

        //- Number of nodes in random-number box [node]
        //- Random-number sets within box are filtered with filterCoeffs_
        //- (e1u, e1v, e1w, e2u, e2v, e2w, e3u, e3v, e3w)
        const List<label> lenRandomBox_;

        //- Convenience factors for 2-D random-number box [-]
        const List<label> randomBoxFactors2D_;

        //- Convenience factors for 3-D randomNum box [-]
        const List<label> randomBoxFactors3D_;

        //- Index to the first elem of last plane of random-number box [-]
        const List<label> iNextToLastPlane_;

        //- Random-number sets distributed over a 3-D box (u, v, w)
        List<List<scalar>> randomBox_;

        //- Filter coefficients corresponding to L_ [-]
        const List<List<scalar>> filterCoeffs_;

        //- Filter-applied random-number sets [m/s] (effectively turb plane)
        List<List<scalar>> filteredRandomBox_;

        //- Filter-applied previous-time-step velocity field [m/s] used in FSM
        vectorField U0_;

        //- Run-time function selector between the original and reduced methods
        const std::function
        <
            void(turbulentDigitalFilterInletFvPatchVectorField*, vectorField&)
        > computeVariant;


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
        Form generateRandomSet(const label len);

        //- Compute nearest cell index-pairs between turbulence plane and patch
        List<Pair<label>> patchIndexPairs();

        //- Check R_ (mapped on actual mesh) for mathematical domain errors
        void checkRTensorRealisable() const;

        //- Compute Lund-Wu-Squires transformation
        //  (Klein et al., 2003, Eq. 5)
        symmTensorField computeLundWuSquires() const;

        //- Compute patch-normal into the domain
        vector computePatchNormal() const;

        //- Compute initial (first time-step) mass/vol flow rate based on UMean_
        scalar computeInitialFlowRate() const;

        //- Convert streamwise integral length scales to integral time scales
        //- via Taylor's frozen turbulence hypothesis
        void convertToTimeScale(tensor& L) const;

        //- Convert length scale phys. unit to turbulence plane mesh-size unit
        //- (Klein et al., 2003, Eq. 13)
        tensor convertScalesToGridUnits(const tensor& L) const;

        //- Resource allocation functions for the convenience factors
        List<label> initLenRandomBox() const;
        List<label> initBoxFactors2D() const;
        List<label> initBoxFactors3D() const;
        List<label> initBoxPlaneFactors() const;

        //- Compute various convenience factors for random-number box
        List<List<scalar>> fillRandomBox();

        //- Compute filter coeffs once per simulation
        //  (Klein et al., 2003, Eq. 14)
        List<List<scalar>> computeFilterCoeffs() const;

        //- Discard current time-step random-box plane (closest to patch) by
        //- shifting from the back to the front, and dd new plane to the back
        void rndShiftRefill();

        //- Map two-point correlated random-number sets on patch based on chosen
        //- mapping method
        void mapFilteredRandomBox(vectorField& U);

        //- Map R_ on patch
        void embedOnePointCorrs(vectorField& U) const;

        //- Map UMean_ on patch
        void embedMeanVelocity(vectorField& U) const;

        //- Correct mass/vol flow rate in (only) streamwise direction
        void correctFlowRate(vectorField& U) const;

        //- 3-D 'valid'-type 'separable' convolution summation algorithm
        //  'Inspired' from Song Ho Ahn's 2-D 'full'-type convolution algorithm
        //  with his permission
        void embedTwoPointCorrs();

        //- Compute the DFM (Klein et al., 2003)
        void computeDFM(vectorField& U);

        //- Compute the reduced DFM (i.e. hybrid FSM-DFM) (Xie-Castro, 2008)
        void computeReducedDFM(vectorField& U);

        //- Compute constList1FSM_ once per simulation
        List<scalar> computeConstList1FSM() const;

        //- Compute constList2FSM_ once per simulation
        List<scalar> computeConstList2FSM() const;

        //- Compute the forward-stepwise method
        //  (Xie-Castro, 2008, Eq. 14)
        void computeFSM(vectorField& U);


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
