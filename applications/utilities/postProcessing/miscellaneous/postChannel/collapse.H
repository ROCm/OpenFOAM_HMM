    scalarField UmeanXvalues = channelIndexing.collapse
    (
        Umean.component(vector::X)
    );

    scalarField UmeanYvalues = channelIndexing.collapse
    (
        Umean.component(vector::Y)
    );

    scalarField UmeanZvalues = channelIndexing.collapse
    (
        Umean.component(vector::Z)
    );

    scalarField RxxValues = channelIndexing.collapse(Rxx);
    scalarField RyyValues = channelIndexing.collapse(Ryy);
    scalarField RzzValues = channelIndexing.collapse(Rzz);
    scalarField RxyValues = channelIndexing.collapse(Rxy, true);

    scalarField pPrime2MeanValues = channelIndexing.collapse(pPrime2Mean);

    scalarField epsilonValues = channelIndexing.collapse(epsilonMean);

    /*
    scalarField nuMeanValues = channelIndexing.collapse(nuMean);
    scalarField nuPrimeValues = channelIndexing.collapse(nuPrime);

    scalarField gammaDotMeanValues = channelIndexing.collapse(gammaDotMean);
    scalarField gammaDotPrimeValues = channelIndexing.collapse(gammaDotPrime);
    */

    scalarField urmsValues = sqrt(mag(RxxValues));
    scalarField vrmsValues = sqrt(mag(RyyValues));
    scalarField wrmsValues = sqrt(mag(RzzValues));

    scalarField kValues
        = 0.5*(sqr(urmsValues) + sqr(vrmsValues) + sqr(wrmsValues));



    scalarField y = channelIndexing.y(mesh.C());

    makeGraph(y, UmeanXvalues, "Uf", R.path(), gFormat);
    makeGraph(y, urmsValues, "u", R.path(), gFormat);
    makeGraph(y, vrmsValues, "v", R.path(), gFormat);
    makeGraph(y, wrmsValues, "w", R.path(), gFormat);
    makeGraph(y, RxyValues, "uv", R.path(), gFormat);
    makeGraph(y, kValues, "k", R.path(), gFormat);

    makeGraph(y, pPrime2MeanValues, "pPrime2Mean", R.path(), gFormat);

    makeGraph(y, epsilonValues, "epsilon", R.path(), gFormat);

    /*
    makeGraph(y, nuMeanValues, "nu", R.path(), gFormat);
    makeGraph(y, nuPrimeValues, "nuPrime", R.path(), gFormat);
    makeGraph(y, gammaDotMeanValues, "gammaDot", R.path(), gFormat);
    makeGraph(y, gammaDotPrimeValues, "gammaDotPrime", R.path(), gFormat);
    */
