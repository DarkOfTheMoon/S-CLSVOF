    // solve Level-Set function as the re-initialization equation
    Info<< "solve the reinitialization equation"
        << nl << endl;

    psi == psi0;

    for(int corr=0; corr<int(epsilon.value()/deltaTau.value()); corr++)
    {
        psi = psi + psi0/mag(psi0)*(1-mag(fvc::grad(psi)*dimChange))*deltaTau;
        psi.correctBoundaryConditions();
    }

    // update Dirac function
    forAll(mesh.cells(),celli)
    {
        if(mag(psi[celli]) > epsilon.value())
        {
            delta[celli] = 0;
        }
        else
        {
            delta[celli] = 0.5/epsilon.value()*(1.0+Foam::cos(M_PI*psi[celli]/epsilon.value()));
        }
    };

    // update Heaviside function
    forAll(mesh.cells(),celli)
    {
        if(psi[celli] < -epsilon.value())
        {
            H[celli] = 0.;
        }
        else if(epsilon.value() < psi[celli])
        {
            H[celli] = 1.;
        }
        else
        {
            H[celli] = 0.5*(1+psi[celli]/epsilon.value()+Foam::sin(M_PI*psi[celli]/epsilon.value())/M_PI);
        }
    };

