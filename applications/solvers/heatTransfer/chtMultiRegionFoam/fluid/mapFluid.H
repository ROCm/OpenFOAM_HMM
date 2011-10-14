h = allToProcMappers[procI].decomposeField(allh(), true);
h.oldTime().timeIndex() = allh().oldTime().timeIndex();
h.correctBoundaryConditions();
