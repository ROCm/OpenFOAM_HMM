2013-01

- to UPstream added allocation of communicators
- added communicators to lduMesh,lduInterface and (indirectly)
  to lduInterfaceField
- gamg agglomeration allocates new communicator for every level
- in all linear solvers/smoothers/preconditioners make they use
  communicator
- added lots of warnings if using unexpected communicator (UPstream::warnComm)
- did LUScalarMatrix for 'directSolveCoarsest'
