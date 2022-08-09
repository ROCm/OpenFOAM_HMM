## Implicit mass conservation test

The implicit mass conservation is switched on in the field (here T) overset BC.

The keyword is used as:

    free
    {
        type            overset;
        massCorrection  true;
        value           uniform 300;
    }

A special net fringe flux is output in this tutorial.

Switch `massCorrection` to `false` to evaluate the mass correction.

