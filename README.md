# OF_USERGROUP
This is a repository to create your own custom turbulence model in openfoam7. Tutourial is based mostly on https://www.tfd.chalmers.se/~hani/kurser/OS_CFD_2019/lectureNotes/ImplementTurbulenceModel.pdf

By Tyler Buchanan

# Instruction on How to Create Custom Turbulence Models:

# Prerequisites
First git clone needed files:
git clone https://github.com/tsb-buchanan/OF_USERGROUP.git


# Docker
It is assuned that openfoam7 is already installed otherwise you can use the docker image and use the following lines:
```bash
Steps to run OpenFOAM-7 on docker 
Make a directory for the openfoam
mkdir openfoam
Copy the files into the ~/openfoam
// Go to the ~/openfoam
$cd ~/openfoam
// Build docker image
sh docker_image_builder.sh
// Start the docker container
sh docker_starter.sh
// Now you are in the container!
// To exit the container (run in the container)
exit
```

# Premade turbulence model (quick install)
In this git is all the files premade. Therefore only commands to be preformed to compile custom turbulence model is:
```bash
cd OF_USERGROUP/TURBFOAM-7/src/TurbulenceModels/incompressible/
wclean
wmakeLnInclude -u ../turbulenceModels
wmake
```
# Implementing Turbulence model from scratch

In this document we will implement a modified version of the kOmegaSST model, named modifiedmykOmegaSST.Before that however,  we will first practice a bit by making our own versions of kEpsilon and kOmegaSST. In fact, there is a problem when only trying to make a new version of kOmegaSST. That is likely due to the special design with the kOmegaSSTBase class. However, the problem disappears  if the kEpsilon model is first implemented. We first start by copying and renaming the kEpsilon and kOmegaSST folders, files and classes:
```bash
foam 
cp -r --parents src/TurbulenceModels/turbulenceModels/RAS/kEpsilon $WM_PROJECT_USER_DIR 
cd $WM_PROJECT_USER_DIR/src/TurbulenceModels/turbulenceModels/RAS 
mv kEpsilon mykEpsilon 
cd mykEpsilon 
mv kEpsilon.C mykEpsilon.C 
mv kEpsilon.H mykEpsilon.H 
sed -i s/kEpsilon/mykEpsilon/g mykEpsilon.* 

foam 
cp -r --parents src/TurbulenceModels/turbulenceModels/RAS/kOmegaSST $WM_PROJECT_USER_DIR 
cd $WM_PROJECT_USER_DIR/src/TurbulenceModels/turbulenceModels/RAS 
mv kOmegaSST mykOmegaSST 
cd mykOmegaSST 
mv kOmegaSST.C mykOmegaSST.C 
mv kOmegaSST.H mykOmegaSST.H
```
At this point we want to change the name of the class in those files using the sed command. The problem is that there is a string 'kOmegaSSTBase' occurring in the files, so we need to use a sed command that will not modify those strings. Here is a sed command that omits lines that contain the word kOmegaSSTBase (and gives an example of how to also omit lines with the words 'mykOmegaSST' and 'dummy'):

```bash
sed -Ei '/(kOmegaSSTBase|mykOmegaSST|dummy)/!s/kOmegaSST/mykOmegaSST/g' mykOmegaSST.*
```
Two important notes about the above command: 1: I experienced that it is important to have the flags in the order Ei and not iE, since otherwise the omission will not work. 2: According to comments in forums the command will not omit consecutive lines with those words. However, in this particular case it works.

Additionally, you want to keep just kOmegaSST for following line in mykOmegaSST.C that is commented below:
```cpp
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
mykOmegaSST<BasicTurbulenceModel>::mykOmegaSST
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
    Foam::kOmegaSST //Keep as kOmegaSST
    <
        eddyViscosity<RASModel<BasicTurbulenceModel>>,
        BasicTurbulenceModel
    >
```
Now we need to copy a file with macros that tell the compiler which instances of turbulence models to compile, as one of the available template options:
```bash
foam 
cp --parents src/TurbulenceModels/incompressible/turbulentTransportModels/turbulentTransportModels.C $WM_PROJECT_USER_DIR 
cp -r --parents src/TurbulenceModels/incompressible/Make $WM_PROJECT_USER_DIR 
cd $WM_PROJECT_USER_DIR/src/TurbulenceModels/incompressible/turbulentTransportModels 
mv turbulentTransportModels.C myTurbulentTransportModels.C
```
Open myTurbulentTransportModels.C and make sure that the active lines are:
```cpp
#include "turbulentTransportModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// -------------------------------------------------------------------------- //
// Laminar models
// -------------------------------------------------------------------------- //


// -------------------------------------------------------------------------- //
// RAS models
// -------------------------------------------------------------------------- //

#include "mykEpsilon.H" 
makeRASModel(mykEpsilon);

#include "mykOmegaSST.H"
makeRASModel(mykOmegaSST);

//#include "modifiedmykOmegaSST.H"
//makeRASModel(modifiedmykOmegaSST);

// -------------------------------------------------------------------------- //
// LES models
// -------------------------------------------------------------------------- //


// ************************************************************************* //
```

Now we need to make sure that we can compile by modifying Make/files and Make/options. First go there:
```bash
cd $WM_PROJECT_USER_DIR/src/TurbulenceModels/incompressible
```
Then make sure that Make/files contains only:
```bash
turbulentTransportModels/myTurbulentTransportModels.C 
LIB = $(FOAM_USER_LIBBIN)/libmyIncompressibleTurbulenceModels
```
Then add two entries under EXE_INC and one entry under LIB_LIBS in Make/options to connect to the original library:
```cpp
EXE_INC = \
    -I../turbulenceModels/lnInclude \
    -I$(LIB_SRC)/transportModels \
    -I$(LIB_SRC)/finiteVolume/lnInclude \
    -I$(LIB_SRC)/meshTools/lnInclude \
    -I$(LIB_SRC)/TurbulenceModels/turbulenceModels/lnInclude \
    -I$(LIB_SRC)/TurbulenceModels/incompressible/lnInclude 

LIB_LIBS = \
    -lincompressibleTransportModels \
    -lturbulenceModels \
    -lfiniteVolume \
    -lmeshTools \
    -lincompressibleTurbulenceModels
```
In order to compile we need to first create an lnInclude directory for the available turbulence models using the wmakeLnInclude command (this must be repeated if more turbulence models are added later), and the we can compile with wmake as usual:
```bash
wmakeLnInclude -u ../turbulenceModels 
wmake
```
# Test mykOmegaSST

In OF_USERGROUP/platforms/Case_Studies/Periodic_Hill/00Baseline/constant/turbulenceProperties change kOmegaSST to mykOmegaSST.

Then at the begining of system/controlDict add:
```cpp
libs ("libmyIncompressibleTurbulenceModels.so");
```
Run (You can just run ./run.sh in 00Baseline folder) with the two new turbulence models (Or just check mykOmegaSST) and check the log-file that they are selected. Try also with a dummy entry and see that both the original and new turbulence models are available for simpleFoam in this particular case (which has the libs-entry in controlDict).


# Implement the modifiedmykOmegaSST model

We will in the rest of the document go through the basic steps of implementing a new kOmegaSST model, named modifiedmykOmegaSST. The model is data-driven turbulence model that corrects the anisotropy reynolds stress (bijDelta) and turbulent production (kDeficit) in the kOmegaSST turbulence model. The Explict Algebraic Stress models are obtained using SpaRTA described in the following papers:

SpaRTA: https://doi.org/10.1007/s10494-019-00089-x

Models used in Tutorial: http://resolver.tudelft.nl/uuid:324d4b2d-bf58-40a0-b60d-4e2e0b992797

Now to begin implementing the modifiedmykOmegaSST, we copy mykOmegaSST, rename directory, files and class names, and finally update the lnInclude directory:
```bash
cd $WM_PROJECT_USER_DIR/src/TurbulenceModels/turbulenceModels/RAS 
cp -r mykOmegaSST modifiedmykOmegaSST
mv modifiedmykOmegaSST/mykOmegaSST.C modifiedmykOmegaSST/modifiedmykOmegaSST.C 
mv modifiedmykOmegaSST/mykOmegaSST.H modifiedmykOmegaSST/modifiedmykOmegaSST.H 
sed -i s/mykOmegaSST/modifiedmykOmegaSST/g modifiedmykOmegaSST/modifiedmykOmegaSST.* 
```


Additionally, you want to keep just kOmegaSST for following line in modifiedmykOmegaSST.C that is commented below:
```cpp
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
mykOmegaSST<BasicTurbulenceModel>::modifiedmykOmegaSST
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
    Foam::kOmegaSST //Keep as kOmegaSST
    <
        eddyViscosity<RASModel<BasicTurbulenceModel>>,
        BasicTurbulenceModel
    >
```

Now add modifiedmykOmegaSST the same way as mykOmegaSST in myTurbulentTransportModels.C (see previous section). Before we compile, we need to make the compiler aware that we have done modifications. The reason for this is that we did not modify any file that is listed in Make/files. We use the touch command to change the time-stamp of that file, so that the compiler will compile.
```bash
cd $WM_PROJECT_USER_DIR/src/TurbulenceModels/incompressible 
wmakeLnInclude -u ../turbulenceModels
touch turbulentTransportModels/myTurbulentTransportModels.C 
wmake
```
Notice that the modifiedmykOmegaSST model is not shown in the compilation output, but it is. The compilation will take more time than expected when we only added one additional turbulence model, and that is because all the templated alternatives mentioned in the above file must be recompiled. Right now, it means that kEpsilon, kOmegaSST, and modifiedmykOmegaSST are all compiled.

Now we need to modifiy modifiedmykOmegaSST with correction models. To do this we need to add the flowing code to modifiedmykOmegaSST.C. In this section, we add new IO objects for model corrections and propagation parameters. 
```cpp
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
```

After adding the IOobjects we now work on the turbulence model. First we will add a modified divDevReff fuction. divDevReff is a term that represents the divergence of the deviatoric part of the effective stress tensor. We modify this term by adding a nonlinear part to original term in kOmegaSSTBASE define below: 

\[ 
fvc::div\left( \text{dev}\left(2 \cdot k \cdot \mathbf{b_{ij\Delta}}\right) \cdot \text{useRST} \cdot \xi \right)
\]


Add the following lines to the bottom of modifiedmykOmegaSST.C:
```cpp
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
```

With the modifed divDevReff() function added, we now will add the modified k and omega equations allong with the Explict Alegbraic Stress Models defined in http://resolver.tudelft.nl/uuid:324d4b2d-bf58-40a0-b60d-4e2e0b992797. The following code can be added to modifiedmykOmegaSST.C:

```cpp
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
    kDeficit_ = 0.043 * epsilon_; 

    
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
```

With all the code snippits added to modifiedmykOmegaSST.C we also need to make a few small additions to modifiedmykOmegaSST.H. The first additions is initializing the added IOobjects described in the .C file: 

```cpp
protected: //Need to add before public:

        // Control of corrections
            dimensionedScalar useRST_;
            dimensionedScalar usekDeficit_;
            dimensionedScalar rampStartTime_;
            dimensionedScalar rampEndTime_;
            dimensionedScalar xi_;  // The ramp value - proportion of corrections to include

        // Frozen fields
            volScalarField epsilon_;
            volScalarField kDeficit_;
            volSymmTensorField bijDelta_;

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;
```

Once these IOojects are added, in modifiedmykOmegaSST.H, add at the end of the class declaration (before the destructor and before the curly bracket that is followed by a semicolon):

```cpp
    //- Solve the turbulence equations and correct the turbulence viscosity
    virtual void correct();

In addtion add divDevReff call after the destructor shown below:
    
    //- Destructor
    virtual ~modifiedmykOmegaSST()
    {}

    tmp<Foam::fvVectorMatrix> divDevReff(volVectorField& U) const;
```

We can now save both the .C and .H file and compile the newly completed turbulence model by using the following lines:
```bash
cd $WM_PROJECT_USER_DIR/src/TurbulenceModels/incompressible 
wclean
wmakeLnInclude -u ../turbulenceModels
wmake
```
After this we can test the new model by going to 01Modified in Case_Studies folder and running the run.sh script. When going into the constant/turbulenceProperties you will see the follwong additons to the file below. This is added to determine when the correction models are added and how much of the correction is being added. (Note case may already be ran so deleted *00 folders)
```cpp
simulationType RAS;

RAS
{
    RASModel        modifiedmykOmegaSST; //mykOmegaSST
    turbulence      on;
    printCoeffs     on;


    useRST          1.0; //Amount of bijDelta correction 1 = 100% and 0 = 0%
    usekDeficit     1.0; //Amount of kDeficit correction 1 = 100% and 0 = 0%
    rampStartTime   2000; // iteration to start adding correction
    rampEndTime     3000; // iteration to reach full corrections defined by useRST and usekDeficit

}
```

After running the case you can visualize both the baseline, modified case and compare to LES data saved as /home/tbuchanan/OF_USERGROUP/Case_Studies/Periodic_Hill/00Baseline/Reference_solution.SST.2000 using the paraFoam function. Make sure to cp -r Reference_solution.SST.2000 2000 so paraFoam can read it. 
