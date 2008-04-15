int CheckEntityParent( CCMIOID parent, CCMIOEntity childType )
{
    if ((childType == kCCMIOProcessor && parent.type != kCCMIOState) ||
	(childType == kCCMIOBoundaryFaces && parent.type != kCCMIOTopology) ||
	(childType == kCCMIOCellType && parent.type != kCCMIOProblemDescription) ||
	(childType == kCCMIOBoundaryRegion && parent.type != kCCMIOProblemDescription) ||
	(childType == kCCMIOInterfaces && parent.type != kCCMIOCells))
	return(0);
    return(1);
}


int CheckEntityParent (CCMIOID parent, CCMIOEntity childType)
{
  switch (childType)
  {
    case kCCMIOProcessor:
      return (parent.type != kCCMIOState) ? 0 : 1;

    case kCCMIOBoundaryFaces:
      return (parent.type != kCCMIOTopology) ? 0 : 1;

    case kCCMIOCellType:
      return (parent.type != kCCMIOProblemDescription) ? 0 : 1;

    case kCCMIOBoundaryRegion:
      return (parent.type != kCCMIOProblemDescription) ? 0 : 1;

    case kCCMIOInterfaces:
      return (parent.type != kCCMIOCells &&
              parent.type != kCCMIOTopology) ? 0 : 1;

    default:
      return 1;
  }
  return 1;
}
