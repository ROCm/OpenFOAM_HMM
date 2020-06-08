/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::functionObjects::STDMD

Group
    grpFieldFunctionObjects

Description
    (Beta Release) STDMD (i.e. Streaming Total Dynamic Mode Decomposition) is
    a variant of a data-driven dimensionality reduction method.

    STDMD is being used as a mathematical post-processing tool to compute
    a set of dominant modes out of a given flow (or dataset) each of which is
    associated with a constant frequency and decay rate, so that dynamic
    features of a given flow may become interpretable, and tractable.
    Among other Dynamic Mode Decomposition (DMD) variants, STDMD is presumed
    to provide the general DMD method capabilities alongside economised and
    feasible memory and CPU usage.

    The code implementation corresponds to Figs. 15-16 of the first citation
    below, more broadly to Section 2.4.

    References:
    \verbatim
        STDMD and mode-sorting algorithms (tags:K, HRDC, KZ, HWR):
            Kiewat, M. (2019).
            Streaming modal decomposition approaches for vehicle aerodynamics.
            PhD thesis. Munich: Technical University of Munich.
            URL:mediatum.ub.tum.de/doc/1482652/1482652.pdf

            Hemati, M. S., Rowley, C. W.,
            Deem, E. A., & Cattafesta, L. N. (2017).
            De-biasing the dynamic mode decomposition
            for applied Koopman spectral analysis of noisy datasets.
            Theoretical and Computational Fluid Dynamics, 31(4), 349-368.
            DOI:10.1007/s00162-017-0432-2

            Kou, J., & Zhang, W. (2017).
            An improved criterion to select
            dominant modes from dynamic mode decomposition.
            European Journal of Mechanics-B/Fluids, 62, 109-129.
            DOI:10.1016/j.euromechflu.2016.11.015

            Hemati, M. S., Williams, M. O., & Rowley, C. W. (2014).
            Dynamic mode decomposition for large and streaming datasets.
            Physics of Fluids, 26(11), 111701.
            DOI:10.1063/1.4901016

        Parallel classical Gram-Schmidt process (tag:Ka):
            Katagiri, T. (2003).
            Performance evaluation of parallel
            Gram-Schmidt re-orthogonalization methods.
            In: Palma J. M. L. M., Sousa A. A., Dongarra J., Hernández V. (eds)
            High Performance Computing for Computational Science — VECPAR 2002.
            Lecture Notes in Computer Science, vol 2565, p. 302-314.
            Berlin, Heidelberg: Springer.
            DOI:10.1007/3-540-36569-9_19

        Parallel direct tall-skinny QR decomposition (tags:BGD, DGHL):
            Benson, A. R., Gleich, D. F., & Demmel, J. (2013).
            Direct QR factorizations for
            tall-and-skinny matrices in MapReduce architectures.
            2013 IEEE International Conference on Big Data.
            DOI:10.1109/bigdata.2013.6691583

            Demmel, J., Grigori, L., Hoemmen, M., & Langou, J. (2012).
            Communication-optimal parallel
            and sequential QR and LU factorizations.
            SIAM Journal on Scientific Computing, 34(1), A206-A239.
            DOI:10.1137/080731992

        DMD properties:
            Brunton S. L. (2018).
            Dynamic mode decomposition overview.
            Seattle, Washington: University of Washington.
            youtu.be/sQvrK8AGCAo (Retrieved:24-04-20)
    \endverbatim

    Operands:
    \table
      Operand      | Type              | Location
      input        | {vol,surface}\<Type\>Field <!--
               --> | $FOAM_CASE/\<time\>/\<inpField\>
      output file  | dat                        <!--
               --> | $FOAM_CASE/postProcessing/\<FO\>/\<time\>/\<file\>(s)
      output field | volScalarField(s) | $FOAM_CASE/\<time\>/\<outField\>(s)
    \endtable

    where \c \<Type\>=Scalar/Vector/SphericalTensor/SymmTensor/Tensor.

    Output fields:
    \verbatim
      modeReal<modeIndex><field> | Real part of a mode field
      modeImag<modeIndex><field> | Imaginary part of a mode field
    \endverbatim

    Output files:
    \verbatim
      uSTDMD.dat                 | Unfiltered STDMD output
      STDMD.dat                  | Filtered STDMD output
    \endverbatim

    wherein for each mode, the following quantities are output into files:
    \verbatim
      Frequency
      Magnitude
      Amplitude (real part)
      Amplitude (imaginary part)
      Eigenvalue (real part)
      Eigenvalue (imaginary part)
    \endverbatim

Note
    Operations on boundary fields, e.g. \c wallShearStress, are currently not
    available.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    STDMD1
    {
        // Mandatory entries (unmodifiable)
        type                STDMD;
        libs                (fieldFunctionObjects);
        field               <inpField>;

        // Conditional mandatory entries (unmodifiable)
        // either of the below
        stdmdInterval       5.5;
        executeInterval     10;

        // Optional entries (unmodifiable)
        modeSorter          kiewat;
        nModes              50;
        maxRank             50;
        nGramSchmidt        5;
        fMin                0;
        fMax                1000000000;

        // Optional entries (run-time modifiable, yet not recommended)
        testEigen           false;
        dumpEigen           false;
        minBasis            0.00000001;
        minEVal             0.00000001;
        absTol              0.001;
        relTol              0.0001;
        sortLimiter         500;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: STDMD                   | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      field        | Name of the operand field          | word |  yes  | -
      stdmdInterval | STDMD time-step size [s]   | scalar| conditional | <!--
                  --> executeInterval*(current time-step of the simulation)
      modeSorter   | Mode-sorting algorithm variant     | word | no  | firstSnap
      nModes | Number of output modes in input freq range | label | no | GREAT
      maxRank  | Max columns in accumulated matrices      | label | no | GREAT
      nGramSchmidt | Number of Gram-Schmidt iterations    | label | no | 5
      fMin         | Min (non-negative) output frequency  | label | no | 0
      fMax         | Max output frequency                 | label | no | GREAT
      testEigen    | Flag to verify eigendecompositions   | bool  | no | false
      dumpEigen    | Flag to log operands of eigendecomps | bool  | no | false
      minBasis | Orthogonal basis expansion threshold     | scalar| no | 1e-8
      minEVal  | Min EVal for below EVals are omitted     | scalar| no | 1e-8
      absTol   | Min abs tol in eigendecomposition tests  | scalar| no | 1e-4
      relTol   | Relative tol in eigendecomposition tests | scalar| no | 1e-6
      sortLimiter | Maximum allowable magnitude for <!--
               -->  mag(eigenvalue)*(number of STDMD steps) to be used in <!--
               -->  modeSorter={kiewat,kouZhang} to avoid overflow errors <!--
               -->                                        | scalar | no | 500.0
    \endtable

    Options for the \c modeSorter entry:
    \verbatim
      kiewat    | Modified weighted amplitude scaling method
      kouZhang  | Weighted amplitude scaling method
      firstSnap | Method of first snapshot amplitude magnitude
    \endverbatim

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link writeFile.H \endlink

    Usage by \c postProcess utility is not available.

Note
    - To specify the STDMD time-step size (not necessarily equal to the
    time step of the simulation), entries of either \c stdmdInterval or
    \c executeInterval must be available (highest to lowest precedence). While
    \c stdmdInterval allows users to directly specify the STDMD time-step size
    in seconds, in absence of \c stdmdInterval, for convenience,
    \c executeInterval allows users to compute the STDMD time-step internally
    by multiplying itself with the current time-step size of the simulation.
    - Limitation: Restart is currently not available since intermediate writing
    of STDMD matrices are not supported.
    - Limitation: Non-physical input (e.g. full-zero fields) may upset STDMD.
    - Warning: DMD is an active research area at the time of writing; therefore,
    there would be cases whereat oddities might be encountered.
    - Warning: This STDMD release is the \b beta release; therefore,
    small-to-medium changes in input/output interfaces and internal structures
    should be expected in the next versions.
    - Warning: In \c modeSorter algorithms of \c kiewat and \c kouZhang, there
    exists a function of \f$power(x,y)\f$ where \c x is the magnitude of an
    eigenvalue, and \c y is the number of STDMD snapshots. This function poses
    a risk of overflow errors since, for example, if \c x ends up above 1.5 with
    500 snapshots or more, this function automatically throws the floating point
    error meaning overflow. Therefore, the domain-range of this function is
    heuristically constrained by the optional entry \c sortLimiter which skips
    the evaluation of this function for a given
    mag(eigenvalue)*(no. of STDMD steps), i.e. x*y, whose magnitude is larger
    than \c sortLimiter.

See also
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::writeFile
    - ExtendedCodeGuide::functionObjects::field::STDMD

SourceFiles
    STDMD.C
    STDMDTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_STDMD_H
#define functionObjects_STDMD_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "EigenMatrix.H"
#include "QRMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                              Class STDMD Declaration
\*---------------------------------------------------------------------------*/

class STDMD
:
    public fvMeshFunctionObject,
    public writeFile
{

    typedef SquareMatrix<scalar> SMatrix;
    typedef SquareMatrix<complex> SCMatrix;
    typedef RectangularMatrix<scalar> RMatrix;
    typedef RectangularMatrix<complex> RCMatrix;


    // Private Enumerations

        //- Options for the mode-sorting algorithm
        enum modeSorterType
        {
            KIEWAT,          //!< "Modified weighted amplitude scaling method"
            KOU_ZHANG,       //!< "Weighted amplitude scaling method"
            FIRST_SNAPSHOT   //!< "Method of first snapshot amplitude magnitude"
        };

        //- Names for modeSorterType
        static const Enum<modeSorterType> modeSorterTypeNames;


    // Private Data

        //- Mode-sorting algorithm (default=modeSorterType::FIRST_SNAPSHOT)
        const enum modeSorterType modeSorter_;

        //- Name of the operand volume or surface field
        const word fieldName_;

        //- Flag: First execution-step initialisation
        bool initialised_;

        //- Flag: Verify eigendecompositions by using theoretical relations
        bool testEigen_;

        //- Flag: Write operands of eigendecompositions to log
        //  To activate it, testEigen=true
        bool dumpEigen_;

        //- Number of output modes within input frequency range
        //- starting from the most energetic mode
        const label nModes_;

        //- Maximum allowable matrix column for Qz_ and Gz_
        //  Qz_ is assumed to always have full-rank, thus Qz_.n() = rank
        const label maxRank_;

        //- Number of maximum iterations of the classical Gram-Schmidt procedure
        const label nGramSchmidt_;

        //- Min frequency: Output only entries corresponding to freqs >= fMin
        const label fMin_;

        //- Max frequency: Output only entries corresponding to freqs <= fMax
        const label fMax_;

        //- Number of base components of input field, e.g. 3 for vector
        label nComps_;

        //- Number of elements in a snapshot
        //  A snapshot is an input dataset to be processed per execution-step
        label nSnap_;

        //- Current execution-step index of STDMD,
        //- not necessarily that of the simulation
        label currIndex_;

        //- Maximum allowable magnitude for mag(eigenvalue)*(number of
        //- STDMD steps) to be used in modeSorter={kiewat,kouZhang} to avoid
        //- overflow errors
        scalar sortLimiter_;

        //- Min absolute tolerance being used in eigendecomposition tests
        scalar absTol_;

        //- Relative tolerance being used in eigendecomposition tests
        scalar relTol_;

        //- Min value to execute orthogonal basis expansion of Qz_ and Gz_
        scalar minBasis_;

        //- STDMD time-step size that equals to
        //- (executeInterval of STDMD)*(deltaT of simulation) [s]
        scalar dt_;

        //- Min EVal magnitude threshold below where EVals are omitted
        scalar minMagEVal_;

        //- L2-norm of column vector z_
        scalar zNorm_;

        //- L2-norm of column vector ez_
        scalar ezNorm_;

        //- Augmented snapshot matrix (effectively a column vector) (K:Eq. 60)
        //  Upper half z_ = current-time snapshot slot
        //  Lower half z_ = previous-time snapshot slot
        RMatrix z_;

        //- Working copy of the augmented snapshot matrix z_
        //- being used in the classical Gram-Schmidt process
        RMatrix ez_;

        //- First-processed snapshot required by the mode-sorting
        //- algorithms at the final output computations (K:p. 43)
        RMatrix X1_;

        //- Accumulated-in-time unitary similarity matrix produced by the
        //- orthonormal decomposition of the augmented snapshot matrix z_
        //  (K:Eq. 60)
        RMatrix Qz_;

        //- Accumulated-in-time (squared) upper triangular matrix produced by
        //- the orthonormal decomposition of the augmented snapshot matrix z_
        //  (K:Eq. 64)
        SMatrix Gz_;

        //- Upper half of Qz_
        RMatrix QzUH_;

        //- Lower half of Qz_
        RMatrix QzLH_;

        //- Moore-Penrose pseudo-inverse of R produced by
        //- the QR decomposition of the last time-step QzUH_
        RMatrix RxInv_;

        //- Projected STDMD operator (K:Eq. 78)
        SMatrix Ap_;

        //- Output eigenvectors
        RCMatrix oEVecs_;

        //- Output eigenvalues
        List<complex> oEVals_;

        //- Output amplitudes
        List<complex> oAmps_;

        //- Output (non-negative) frequencies
        List<scalar> oFreqs_;

        //- Indices of oFreqs_ where freqs are
        //- non-negative and within [fMin_, fMax_]
        DynamicList<label> iFreqs_;

        //- Output (descending) magnitudes of (complex) amplitudes
        List<scalar> oMags_;

        //- Indices of oMags_
        List<label> iMags_;


    // Private Member Functions

        // Process

        //- Move previous-time snapshot into previous-time slot in z_
        //- and copy new current-time snapshot into current-time slot in z_
        template<class GeoFieldType>
        bool getSnapshot();

        //- Get the input field type to be processed by snapshot()
        template<class Type>
        bool getSnapshotType();

        //- Get the number of base components of input field
        template<class GeoFieldType>
        bool getComps();

        //- Return (parallel) L2-norm of a given column vector
        scalar parnorm(const RMatrix& colVector) const;

        //- Move the current-time snapshot into the previous-time snapshot in z_
        //- and copy the new field into the current-time snapshot
        void snapshot();

        //- Initialise all working matrices at the first execution-step
        void init();

        //- Initialise Qz_, Gz_ (both require the first two snapshots) and X1_
        void initBasis();

        //- Execute (parallel) classical Gram-Schmidt
        //- process to orthonormalise ez_ (Ka:Fig. 5)
        void GramSchmidt();

        //- Expand orthonormal bases Qz_ and Gz_ by stacking a column
        //- (ez_/ezNorm_) to Qz_, and a row (Zero) and column (Zero)
        //- to Gz_ if (minBasis_ < ezNorm_/zNorm_)
        void expandBasis();

        //- Update Gz_ before the potential orthonormal basis compression
        void updateGz();

        //- Compress orthonormal basis for Qz_ and Gz_ if (maxRank_ < Qz_.n())
        void compressBasis();


        // Postprocess

        //- Return a new List containing elems of List at indices
        template<class Type>
        void filterIndexed
        (
            List<Type>& lst,
            const UList<label>& indices
        );

        //- Return a new Matrix containing columns of Matrix at indices
        template<class MatrixType>
        void filterIndexed
        (
            MatrixType& lst,
            const UList<label>& indices
        );

        //- Compute global Ap (K:Eq. 78)
        void calcAp();

        //- Compute eigenvalues and eigenvectors
        void calcEigen();

        //- Weak-type floating-point comparison
        //  bit.ly/2Trdbgs (Eq. 2), and PEP-485
        bool close
        (
            const scalar s1,
            const scalar s2,
            const scalar absTol = 0,   //<! comparisons near zero
            const scalar relTol = 1e-8 //<! e.g. vals the same within 8 decimals
        ) const;

        //- Test real/complex eigenvalues by using
        //- the theoretical relation: (sum(eigenvalues) - trace(A) ~ 0)
        void testEigenvalues
        (
            const SquareMatrix<scalar>& A,
            const DiagonalMatrix<scalar>& EValsRe
        ) const;

        //- Test real eigenvectors by using the theoretical relation:
        //- ((A & EVec - EVal*EVec).norm() ~ 0)
        void testEigenvectors
        (
            const SquareMatrix<scalar>& A,
            const DiagonalMatrix<scalar>& EValsRe,
            const SquareMatrix<scalar>& EVecs
        ) const;

        //- Test complex eigenvectors by using the theoretical relation:
        //- ((A & EVec - EVal*EVec).norm() ~ 0)
        void testEigenvectors
        (
            const SquareMatrix<scalar>& A,
            const List<complex>& EVals,
            const RectangularMatrix<complex>& EVecs
        ) const;

        //- Remove any eigenvalues whose magnitude is smaller than
        //- minMagEVal_ while keeping the order of elements the same
        void filterEVals();

        //- Compute and filter frequencies (K:Eq. 81)
        void calcFreqs();

        //- Compute frequency indices
        //  Locate indices where oFreqs_ are
        //  in [fMin_, fMax_], and store them in iFreqs_ indices
        void calcFreqI();

        //- Compute amplitudes
        void calcAmps();

        //- Compute magnitudes
        //  Includes advanced sorting methods using
        //  the chosen weighted amplitude scaling method
        void calcMags();

        //- Compute the ordered magnitude indices
        void calcMagI();

        //- Compute modes
        void calcModes();

        //- Eigenvalue weighted amplitude scaling (KZ:Eq. 33)
        //- Modified eigenvalue weighted amplitude scaling (K)
        scalar sorter
        (
            const List<scalar>& weight,
            const complex& amplitude,
            const complex& eval,
            const scalar modeNorm
        ) const;

        //- Output file header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Filter objects according to iFreqs_ and iMags_
        void filterOutput();

        //- Write unfiltered/filtered data
        void writeOutput(OFstream& os) const;

        //- Compute STDMD output
        void calcOutput();


public:

    //- Runtime type information
    TypeName("STDMD");


    // Constructors

        //- No default construct
        STDMD() = delete;

        //- Construct from Time and dictionary
        STDMD
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        STDMD(const STDMD&) = delete;

        //- No copy assignment
        void operator=(const STDMD&) = delete;


    //- Destructor
    virtual ~STDMD() = default;


    // Member Functions

        //- Read STDMD settings
        virtual bool read(const dictionary&);

        //- Execute STDMD
        virtual bool execute();

        //- Write STDMD output
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "STDMDTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
