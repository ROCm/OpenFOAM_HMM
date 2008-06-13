To run all               : Allrun
To run all one iteration : Alltest
To clean all             : Allclean

The structure now is that a case that only requires blockMesh
and application does not need an Allrun or Allclean script.
Only if running is special it requires an Allrun. Similarly
if cleaning is non-standard.


                                run clean
boundaryFoam
    boundaryLaunderSharma       ok
    boundaryWallFunctions       ok
bubbleFoam
    bubbleColumn                ok
buoyantFoam
    hotRoom                     ok
buoyantSimpleFoam
    hotRoom                     ok
buoyantSimpleRadiationFoam
    hotRadiationRoom            ok
cavitatingFoam
    nozzle2D                    No cavitatingFoam solver.
channelOodles
    channel395                  ok
coodles
    pitzDaily                   ok
dieselFoam
    aachenBomb
dnsFoam
    boxTurb16
electrostaticFoam
    chargedWire                 ok
engineFoam
    kivaTest
financialFoam
    europeanCall                ok
icoDyMFoam
    movingCone                  ok  ok
icoFoam
    cavity                      ok  ok
    cavityFine                  ok  ok
    cavityGrade                 ok  ok
    cavityHighRe                ok  ok
    cavityClipped               ok  ok
    elbow                       ok  ok
interDyMFoam
    damBreakWithObstacle        ok  ok
interFoam
    damBreak                    ok  ok
    damBreakFine                ok  ok
laplacianFoam
    flange                      ok  ok
lesInterFoam
    nozzleFlow2D
mdEquilibrationFoam
    periodicCube
mhdFoam                         ok
    hartmann
MRFSimpleFoam
    mixerVessel2D
multiphaseInterFoam
    damBreak4phase
    damBreak4phaseFine
nonNewtonianIcoFoam
    offsetCylinder
oodles
    pitzDaily
    pitzDailyDirectMapped
potentialFoam
    cylinder
    pitzDaily
rasInterFoam
    damBreak
    damBreakFine
rhoPimpleFoam
    angledDuct
rhoPorousSimpleFoam
    angledDuctExplicit
    angledDuctImplicit
rhopSonicFoam
    shockTube
    wedge15Ma5
rhoSonicFoam
    forwardStep
    shockTube               No setShock application
rhoTurbFoam
    cavity
rhoTurbTwinParcelFoam
    simplifiedSiwek
scalarTransportFoam
    pitzDaily
settlingFoam
    dahl
    tank3D
simpleFoam
    pitzDaily
    pitzDaily3Blocks
    pitzDailyExptInlet
simpleSRFFoam
    mixer
    simpleSRFFoam
solidDisplacementFoam
    plateHole
solidEquilibriumDisplacementFoam
    beamEndLoad
sonicFoam
    forwardStep
    shockTube
sonicLiquidFoam
    decompressionTank
    decompressionTankFine
sonicTurbFoam
    nacaAirfoil
    prism
turbFoam
    cavity
twoPhaseEulerFoam
    bed
    bed2
    bubbleColumn
XiFoam
    moriyoshiHomogeneous
    moriyoshiHomogeneousPart2
Xoodles
    pitzDaily
    pitzDaily3D

