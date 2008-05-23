forAll (pointsOnThisSegment, pS)
{
    const point& ptP = mesh_.points()[pointsOnThisSegment[pS]];

    // Assessing real cells in range is only required on the 1st iteration
    // because they do not change from iteration to iteration.

    if (iterationNo == 0)
    {
        forAll (mesh_.points(), pointIIndex)
        {

            // Compare separation of ptP to all other points in the mesh,
            // add unique reference to cell with any point within rCutMax_
            // to realCellsFoundInRange.

            const point& ptI(mesh_.points()[pointIIndex]);

            if (magSqr(ptP - ptI) <= rCutMaxSqr)
            {
                const labelList& ptICells(mesh_.pointCells()[pointIIndex]);

                forAll(ptICells, pIC)
                {
                    const label cellI(ptICells[pIC]);

                    if(findIndex(realCellsFoundInRange, cellI) == -1)
                    {
                        realCellsFoundInRange.append(cellI);
                    }
                }
            }
        }
    }

    forAll(referredInteractionList, rIL)
    {
        const vectorList& refCellPoints =
            referredInteractionList[rIL].vertexPositions();

        forAll(refCellPoints, rCP)
        {
            if (magSqr(ptP - refCellPoints[rCP]) <= rCutMaxSqr)
            {
                if(findIndex(referredCellsFoundInRange, rIL) == -1)
                {
                    referredCellsFoundInRange.append(rIL);
                }
            }
        }
    }
}
