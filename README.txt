Notes from merging Integration-TUD-addOns and plus.

0)
- Changed
dynamicMesh/polyTopoChange/polyTopoChange/hexRef8/hexRef8.C
to created internal faces out-of-nothing.
- Added the mapNewInternalFaces to dynamicRefineFvMesh
- unset FOAM_SETNAN, FOAM_SIGFPE and run

    testCases/testAMRandLoadBalancing/damBreakWithObstacle

- However when writing V0 we noticed a nan in the V0 field.
  (this doesn't get written anymore but it might indicate a bug)


1) Tried moving additional mapping of surface fields to

    virtual dynamicRefineFvMesh::mapFields(const mapPolyMesh&);

so it would get called from the fvMesh::updateMesh. However
this mapping (explicitly) gets done using unadapted addressing
and it causes the cell addressing to be wrong:


[1] --> FOAM FATAL ERROR:
[1] index 826015320 out of range 0 ... 71079
[1]
[1]     From function void Foam::UList<T>::checkIndex(Foam::label) const [with T = Foam::cell; Foam::label = int]
[1]     in file /home/preston2/mattijs/OpenFOAM/work/OpenFOAM-plus.integration-TUD/src/OpenFOAM/lnInclude/UListI.H at line 106.


From dynamicRefineFvMesh::mapNewInternalFaces : (I think)

                    cell faceOwner = this->cells()[owner[facei]];

So this adaptation should be done in a separate pass


2) Mapping new internal faces:
- currently done by averaging the values of 'properly' mapped
  faces from owner and neighbour.
- this would not work if a cell gets split into 3x3x3
- instead this should be done by some geometric interpolation
  from point values?
- have selection mechanism instead or use oriented flag on surfaceFields

3) Current averaging is valid only for internal faces & only scalar/vector
