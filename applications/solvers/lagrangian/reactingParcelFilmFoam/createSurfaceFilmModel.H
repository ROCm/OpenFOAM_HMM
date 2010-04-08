Info<< "\nConstructing surface film model" << endl;

autoPtr<surfaceFilmModels::surfaceFilmModel>
    tsurfaceFilm(surfaceFilmModels::surfaceFilmModel::New(mesh, g));
surfaceFilmModels::surfaceFilmModel& surfaceFilm = tsurfaceFilm();

