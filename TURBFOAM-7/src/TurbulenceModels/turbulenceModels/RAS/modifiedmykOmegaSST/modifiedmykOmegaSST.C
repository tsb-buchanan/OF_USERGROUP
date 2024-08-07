/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "modifiedmykOmegaSST.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
modifiedmykOmegaSST<BasicTurbulenceModel>::modifiedmykOmegaSST
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    Foam::kOmegaSST
    <
        eddyViscosity<RASModel<BasicTurbulenceModel>>,
        BasicTurbulenceModel
    >
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    //====================== Propagation parameters ======================
    // Set these to 0 to turn off influence of model corrections, 1 to include fully.
    // The kDeficit and bDelta fields come from models further down.  Values in the
    // interval [0,1] make sense too.
    
    useRST_   // Use the LES RST (via a_ij) in the omega-eqn production
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "useRST",
            this->coeffDict_,
            1.0
        )
    ),
    usekDeficit_  // Use the k-eqn residual to correct the omega production
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "usekDeficit",
            this->coeffDict_,
            1.0
        )
    ),
    // Ramping gradually introduce corrections (both kDeficit and bDelta) to
    // aid solver stability.  Before `rampStartTime` corrections are zero,
    // after `rampEndTime` they are 1.0.  Linear in between.  Default is
    // full correction from beginning.
    rampStartTime_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "rampStartTime",
            this->coeffDict_,
	    dimTime,
            0
        )
    ),
    rampEndTime_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "rampEndTime",
            this->coeffDict_,
	    dimTime,
            1
        )
    ),
    xi_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "xi_ramp",
            this->coeffDict_,
	    dimless,
            1
        )
    ),

    //========================== Fields from frozen solver ================
    epsilon_
    (
        IOobject(
	    "epsilon",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("epsilon", dimensionSet(0,2,-3,0,0,0,0), 0.0) 
    ),   
    kDeficit_
    (
        IOobject(
	    "kDeficit",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("kDeficit", dimensionSet(0,2,-3,0,0,0,0), 0.0) 
    ),
    
    
    bijDelta_
    (
        IOobject
        (
            "bijDelta",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor
        (
        "bijDelta",
        dimensionSet(0,0,0,0,0,0,0),
        symmTensor(0,0,0,0,0,0)
        )
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}

template<class BasicTurbulenceModel>
tmp<fvVectorMatrix> modifiedmykOmegaSST<BasicTurbulenceModel>::divDevReff
(
    volVectorField& U
) const
{
    Info << "In: modifiedmykOmegaSST::divDevReff()" << endl;
    return
    (
      - fvc::div((this->alpha_*this->rho_*this->nuEff())*dev2(T(fvc::grad(U))))
      - fvm::laplacian(this->alpha_*this->rho_*this->nuEff(), U)    // linear part
      + fvc::div(dev(2.*this->k_*this->bijDelta_) * useRST_ * xi_)  // non-linear part
    );
}
  

template<class BasicTurbulenceModel>
void modifiedmykOmegaSST<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }
    Info << "In: modifiedmykOmegaSST.correct()" << endl;
    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    // Update omega and G at the wall
    volScalarField& omega_ = this->omega_;
    volScalarField& k_ = this->k_;

    BasicTurbulenceModel::correct();
    // This is the current iteration, 1000, etc.
    const dimensionedScalar time = this->runTime_;
    xi_ = (time < rampStartTime_)? 0.0:
          (time > rampEndTime_)? 1.0:
          (time - rampStartTime_) / (rampEndTime_ - rampStartTime_);
    Info << "Corrections: xi = " << xi_.value() <<
          ", kDeficit factor = " << (xi_*usekDeficit_).value() <<
          ", bijDelta factor = " << (xi_*useRST_).value() << endl;

    volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))()()
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    volSymmTensorField S(symm(tgradU()) / omega_);
    volTensorField W(skew(tgradU()) / omega_);
    volScalarField S2(2*magSqr(symm(tgradU())));


/*******************************Feature Addition**************************/
    //add features for models
    volSymmTensorField T1(S);
    volSymmTensorField T2((symm((S & W) - (W & S))));
    volSymmTensorField T3((symm((S & S) - 0.3333*tr(S & S)*I)));
    volSymmTensorField T4((symm((W & W) - 0.3333*tr(W & W)*I)));
    
    
    
    epsilon_ = (omega_*k_);

/************************************************************************/


/*******************************Model Addition**************************/
    //Model kDeficit_
    kDeficit_ = 0.043 * epsilon_; //ask about constant = 1

    
    //Model bij
    bijDelta_ =  (0.174 * T2) + (0.357 * T3);
    
/************************************************************************/

    volScalarField GbyNu(dev(twoSymm(tgradU())) && tgradU());
    volScalarField::Internal G(this->GName(), nut()*GbyNu);
    volScalarField G2
    (
         "G2",
	 nut*GbyNu - xi_ * useRST_ * (2*(this->k_)*bijDelta_ && tgradU())
    );


    tgradU.clear();

    volScalarField CDkOmega
    (
        (2*this->alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    volScalarField F1(this->F1(CDkOmega));
    volScalarField F23(this->F23());

    {
        volScalarField::Internal gamma(this->gamma(F1));
        volScalarField::Internal beta(this->beta(F1));

        // Turbulent frequency equation
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, omega_)
          + fvm::div(alphaRhoPhi, omega_)
          - fvm::laplacian(alpha*rho*this->DomegaEff(F1), omega_)
         ==
            alpha()*rho()*gamma
           *min
            (
   	     // Production modified due to RST correction
	     G2 / nut(),
                (this->c1_/this->a1_)*this->betaStar_*omega_()
               *max(this->a1_*omega_(), this->b1_*F23()*sqrt(S2()))
            )
	    // Production modified due to k-equation correction
	    + alpha()*rho()*gamma*kDeficit_*omega_()/k_()*(xi_ * usekDeficit_) 
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, omega_)
          - fvm::Sp(alpha()*rho()*beta*omega_(), omega_)
          - fvm::SuSp
            (
                alpha()*rho()*(F1() - scalar(1))*CDkOmega()/omega_(),
                omega_
            )
          + this->Qsas(S2(), gamma, beta)
          + this->omegaSource()
          + fvOptions(alpha, rho, omega_)
        );
	omega_.boundaryFieldRef().updateCoeffs();

        omegaEqn.ref().relax();
        fvOptions.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
        solve(omegaEqn);
        fvOptions.correct(omega_);
        bound(omega_, this->omegaMin_);
    }

    // --------------- Turbulent kinetic energy equation ----------------
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*this->DkEff(F1), k_)
     ==
	// Production modified due to RST correction
        alpha()*rho()*G2  //this->Pk(G)
	// Production modified due to k-equation correction
	+ alpha()*rho()*kDeficit_*(xi_ * usekDeficit_)
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*this->epsilonByk(F1, F23), k_)
      + this->kSource()
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    this->correctNut(S2, F23);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
