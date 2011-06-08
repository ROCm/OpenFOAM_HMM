Info<< "Creating pyrolysis model" << endl;

autoPtr<regionModels::pyrolysisModels::pyrolysisModel> pyrolysis
(
    regionModels::pyrolysisModels::pyrolysisModel::New(mesh)
);
