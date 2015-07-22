{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phic(mag(phi/mesh.magSf()));
    phic = min(interface.cAlpha()*phic, max(phic));
//    surfaceScalarField phir(phic*interface.nHatf());
    surfaceScalarField phir(phic*nVecf);

    surfaceScalarField phiAlpha
    (
        fvc::flux
        (
            phi,
            alpha1,
            alphaScheme
        )
      + fvc::flux
        (
            -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
            alpha1,
            alpharScheme
        )
    );

    MULES::explicitSolve(alpha1, phi, phiAlpha, 1, 0);

    rhoPhi = phiAlpha*(rho1 - rho2) + phi*rho2;
}
