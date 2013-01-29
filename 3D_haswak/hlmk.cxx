/*******************************************************************************
 * 2-fluid equations
 * Same as Maxim's version of BOUT - simplified 2-fluid for benchmarking
 *******************************************************************************/

#include <bout.hxx>
#include <boutmain.hxx>

#include <initialprofiles.hxx>
#include <derivs.hxx>
#include <interpolation.hxx>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// 2D initial profiles
Field2D Ni0, Ti0, Te0, Vi0, phi0, Ve0, rho0, Ajpar0,gradNi0;
Vector2D b0xcv, b0,B0; /// for curvature terms

// 3D evolving fields
Field3D rho, Te, Ni, Ajpar, Vi, Ti;

//Flux                                                                         
Vector3D Gamma,vEB;

// Derived 3D variables
Field3D phi, Apar, Ve, jpar, phi_filt;

// Non-linear coefficients
Field3D nu, mu_i, kapa_Te, kapa_Ti;

// 3D total values
Field3D Nit, Tit, Tet, Vit,rhot;

// pressures
Field3D pei, pe;
Field2D pei0, pe0;

// Metric coefficients
Field2D Rxy, Bpxy, Btxy, hthe,DXRxy;

// parameters
BoutReal Te_x, Ti_x, Ni_x, Vi_x, bmag, rho_s, fmei, AA, ZZ;
BoutReal lambda_ei, lambda_ii;
BoutReal nu_hat, mui_hat, wci, nueix, nuiix;
BoutReal beta_p;

// settings
bool estatic, ZeroElMass; // Switch for electrostatic operation (true = no Apar)

bool noDC,plusDC;
bool nonlinear, haswak,par_damp,transport;
BoutReal zlowpass;
int nzpass;
int MZ;


BoutReal zeff, nu_perp;
bool evolve_rho, evolve_te, evolve_ni, evolve_ajpar, evolve_vi, evolve_ti;
BoutReal ShearFactor;

int phi_flags, apar_flags; // Inversion flags

// Field routines
int solve_phi_tridag(Field3D &r, Field3D &p, int flags);
int solve_apar_tridag(Field3D &aj, Field3D &ap, int flags);


FieldGroup comms; // Group of variables for communications

int physics_init(bool restarting)
{
  Field2D I; // Shear factor 
  
  output.write("Solving 6-variable 2-fluid equations\n");

  /************* LOAD DATA FROM GRID FILE ****************/

  // Load 2D profiles (set to zero if not found)
  GRID_LOAD(Ni0);
  GRID_LOAD(Ti0);
  GRID_LOAD(Te0);
  GRID_LOAD(Vi0);
  GRID_LOAD(Ve0);
  GRID_LOAD(phi0);
  GRID_LOAD(rho0);
  GRID_LOAD(Ajpar0);

  // Load magnetic curvature term
  b0xcv.covariant = false; // Read contravariant components
  mesh->get(b0xcv, "bxcv"); // b0xkappa terms

  // Load metrics
  GRID_LOAD(Rxy);
  GRID_LOAD(Bpxy);
  GRID_LOAD(Btxy);
  GRID_LOAD(hthe);
  mesh->get(mesh->dx,   "dpsi");
  mesh->get(I,    "sinty");
  mesh->get(mesh->zShift, "qinty");

  // Load normalisation values
  GRID_LOAD(Te_x);
  GRID_LOAD(Ti_x);
  GRID_LOAD(Ni_x);
  GRID_LOAD(bmag);

  Ni_x *= 1.0e14;
  bmag *= 1.0e4;

  /*************** READ OPTIONS *************************/

  // Read some parameters
  Options *globalOptions = Options::getRoot();
  Options *options = globalOptions->getSection("2fluid");
  OPTION(options, AA, 2.0);
  OPTION(options, ZZ, 1.0);

  OPTION(options, estatic,     false);
  OPTION(options, ZeroElMass,  false);
  OPTION(options, zeff,        1.0);
  OPTION(options, nu_perp,     0.0);
  OPTION(options, ShearFactor, 1.0);
  
  OPTION(options, noDC, false);
  OPTION(options, plusDC, false);
  OPTION(options,nonlinear,false);
  OPTION(options,haswak,false);
  OPTION(options,par_damp,false);
  OPTION(options,zlowpass,0);
  OPTION(options,transport,true);

  OPTION(options, phi_flags,   0);
  OPTION(options, apar_flags,  0);
  
  (globalOptions->getSection("Ni"))->get("evolve", evolve_ni,    true);
  (globalOptions->getSection("rho"))->get("evolve", evolve_rho,   true);
  (globalOptions->getSection("vi"))->get("evolve", evolve_vi,   true);
  (globalOptions->getSection("te"))->get("evolve", evolve_te,   true);
  (globalOptions->getSection("ti"))->get("evolve", evolve_ti,   true);
  (globalOptions->getSection("Ajpar"))->get("evolve", evolve_ajpar, true);
  
  OPTION(options,MZ,33);
  
  if (zlowpass != 0)
    nzpass = MZ*zlowpass;
  
  if(ZeroElMass)
    evolve_ajpar = false; // Don't need ajpar - calculated from ohm's law

  /************* SHIFTED RADIAL COORDINATES ************/

  if(mesh->ShiftXderivs) {
    ShearFactor = 0.0;  // I disappears from metric
    b0xcv.z += I*b0xcv.x;
  }

  /************** CALCULATE PARAMETERS *****************/

  rho_s = 1.02*sqrt(AA*Te_x)/ZZ/bmag;
  fmei  = 1./1836.2/AA;

  lambda_ei = 24.-log(sqrt(Ni_x)/Te_x);
  lambda_ii = 23.-log(ZZ*ZZ*ZZ*sqrt(2.*Ni_x)/pow(Ti_x, 1.5));
  wci       = (1.0)*9.58e3*ZZ*bmag/AA;
  nueix     = 2.91e-6*Ni_x*lambda_ei/pow(Te_x, 1.5);
  nuiix     = 4.78e-8*pow(ZZ,4.)*Ni_x*lambda_ii/pow(Ti_x, 1.5)/sqrt(AA);
  nu_hat    = zeff*nueix/wci;

  if(nu_perp < 1.e-10) {
    mui_hat      = (3./10.)*nuiix/wci;
  } else
    mui_hat      = nu_perp;

  if(estatic) {
    beta_p    = 1.e-29;
  }else
    beta_p    = 4.03e-11*Ni_x*Te_x/bmag/bmag;

  Vi_x = wci * rho_s;

  /************** PRINT Z INFORMTION ******************/
  
  BoutReal hthe0;
  if(mesh->get(hthe0, "hthe0") == 0) {
    output.write("    ****NOTE: input from BOUT, Z length needs to be divided by %e\n", hthe0/rho_s);
  }

  /************** SHIFTED GRIDS LOCATION ***************/

  // Velocities defined on cell boundaries
  Vi.setLocation(CELL_YLOW);
  Ajpar.setLocation(CELL_YLOW);

  // Apar and jpar too
  Apar.setLocation(CELL_YLOW); 
  jpar.setLocation(CELL_YLOW);

  /************** NORMALISE QUANTITIES *****************/

  output.write("\tNormalising to rho_s = %e\n", rho_s);

  // Normalise profiles
  Ni0 /= Ni_x/1.0e14;
  Ti0 /= Te_x;
  Te0 /= Te_x;
  phi0 /= Te_x;
  Vi0 /= Vi_x;
  

  // Normalise curvature term
  b0xcv.x /= (bmag/1e4);
  b0xcv.y *= rho_s*rho_s;
  b0xcv.z *= rho_s*rho_s;
  
  // Normalise geometry 
  Rxy /= rho_s;
  hthe /= rho_s;
  I *= rho_s*rho_s*(bmag/1e4)*ShearFactor;
  output.write("mesh->dx = %e\n", mesh->dx[0]);
  mesh->dx /= rho_s*rho_s*(bmag/1e4);
  output.write("mesh->dx = %e\n", mesh->dx[0]);
  // Normalise magnetic field
  Bpxy /= (bmag/1.e4);
  Btxy /= (bmag/1.e4);
  mesh->Bxy  /= (bmag/1.e4);
  
  // calculate pressures
  pei0 = (Ti0 + Te0)*Ni0;
  pe0 = Te0*Ni0;

  /**************** CALCULATE METRICS ******************/

  mesh->g11 = (Rxy*Bpxy)^2;
  mesh->g22 = 1.0 / (hthe^2);
  mesh->g33 = (I^2)*mesh->g11 + (mesh->Bxy^2)/mesh->g11;
  mesh->g12 = 0.0;
  mesh->g13 = -I*mesh->g11;
  mesh->g23 = -Btxy/(hthe*Bpxy*Rxy);
  
  mesh->J = hthe / Bpxy;
  
  mesh->g_11 = 1.0/mesh->g11 + ((I*Rxy)^2);
  mesh->g_22 = (mesh->Bxy*hthe/Bpxy)^2;
  mesh->g_33 = Rxy*Rxy;
  mesh->g_12 = Btxy*hthe*I*Rxy/Bpxy;
  mesh->g_13 = I*Rxy*Rxy;
  mesh->g_23 = Btxy*hthe*Rxy/Bpxy;

  mesh->geometry();
 
  B0.y = 1.0/mesh->J;
  B0.z = 0;
  B0.x = 0;
  B0.covariant = false;
  
  b0 = B0/abs(B0);

  b0xcv = b0 ^ V_dot_Grad(b0,b0);
  DXRxy = DDX(Rxy);
  gradNi0 = DDX(Ni0);
  
  /**************** SET EVOLVING VARIABLES *************/

  // Tell BOUT++ which variables to evolve
  // add evolving variables to the communication object
  if(evolve_rho) {
    bout_solve(rho, "rho");
    comms.add(rho);
    output.write("rho\n");
  }else
    initial_profile("rho", rho);

  if(evolve_ni) {
    bout_solve(Ni, "Ni");
    comms.add(Ni);
    output.write("ni\n");
  }else
    initial_profile("Ni", Ni);

  if(evolve_te) {
    bout_solve(Te, "Te");
    comms.add(Te);
    output.write("te\n");
  }else
    initial_profile("Te", Te);

  if(evolve_ajpar) {
    bout_solve(Ajpar, "Ajpar");
    comms.add(Ajpar);
    output.write("ajpar\n");
  }else {
    initial_profile("Ajpar", Ajpar);
    if(ZeroElMass)
      dump.add(Ajpar, "Ajpar", 1); // output calculated Ajpar
  }

  if(evolve_vi) {
    bout_solve(Vi, "Vi");
    comms.add(Vi);
    output.write("vi\n");
  }else
    initial_profile("Vi", Vi);

  if(evolve_ti) {
    bout_solve(Ti, "Ti");
    comms.add(Ti);
    output.write("ti\n");
  }else
    initial_profile("Ti", Ti);
  
  // Set boundary conditions
  jpar.setBoundary("jpar");
  phi.setBoundary("phi");
  /************** SETUP COMMUNICATIONS **************/

  // add extra variables to communication
  comms.add(phi);
  comms.add(Apar);

  if(transport) {
    dump.add(Gamma.x,"Gammax",1);
    dump.add(Gamma.y,"Gammay",1);
    dump.add(Gamma.z,"Gammaz",1);
  }

  // Add any other variables to be dumped to file
  dump.add(phi,  "phi",  1);
  dump.add(Apar, "Apar", 1);
  dump.add(jpar, "jpar", 1);
  
  dump.add(ddt(Ni),"ddtNi",1);
  dump.add(ddt(rho),"ddtrho",1);

  dump.add(Ni0, "Ni0", 0);
  
  dump.add(gradNi0,"gradNi0",0);
  dump.add(DXRxy,"DXRxy",0);
  dump.add(Te0, "Te0", 0);
  dump.add(Ti0, "Ti0", 0);

  dump.add(Te_x,  "Te_x", 0);
  dump.add(Ti_x,  "Ti_x", 0);
  dump.add(Ni_x,  "Ni_x", 0);
  dump.add(rho_s, "rho_s", 0);
  dump.add(wci,   "wci", 0);

  dump.add(Nit, "Nit", 0);
  
  return(0);
}

// just define a macro for V_E dot Grad
#define vE_Grad(f, p) ( b0xGrad_dot_Grad(p, f) / mesh->Bxy )
#define ParLaplac(f) (mesh->G2*DDY(f)+mesh->g22*D2DY2(f))


int physics_run(BoutReal t)
{
  // Solve EM fields
  int ncalls = solver->rhs_ncalls;
  rho.applyBoundary();
  //rho and phi_flags come in, phi is set
  solve_phi_tridag(rho, phi, phi_flags); //is this causing issues
  

  //let's just grab the value from the middle
  // phi_filt.allocate();
  // for(int jx=1;jx<mesh->ngx-1;jx++)
  //   for(int jy=0;jy<mesh->ngy;jy++)
  //     for(int jz=0;jz<mesh->ngz;jz++) 
  // 	phi_filt[jx][jy][jz] =phi[2][jy][jz];
  
  // phi = phi_filt;
  phi.applyBoundary();

  if(estatic || ZeroElMass) {
    // Electrostatic operation
    Apar = 0.0;
  }else {
    solve_apar_tridag(Ajpar, Apar, apar_flags); // Linear Apar solver
  }

  // Communicate variables
  mesh->communicate(comms);


  // Update profiles
  Nit = Ni0;// + Ni.DC();
  Tit = Ti0; // + Ti.DC();
  Tet = Te0; // + Te.DC();
  Vit = Vi0;// + Vi.DC();
  rhot = rho0;
  
  if(plusDC) {
    Nit += Ni.DC();
    Tit += Ti.DC();
    Tet += Te.DC();
    Vit += Vi.DC();
    rhot += rho.DC();
  }
  

  // Update non-linear coefficients on the mesh
  nu      = nu_hat * Nit / (Tet^1.5);
  mu_i    = mui_hat * Nit / (Tit^0.5);
  kapa_Te = 3.2*(1./fmei)*(wci/nueix)*(Tet^2.5);
  kapa_Ti = 3.9*(wci/nuiix)*(Tit^2.5);
  
  // note: nonlinear terms are not here
  pei = (Te0+Ti0)*Ni + (Te + Ti)*Ni0;
  pe  = Te0*Ni + Te*Ni0;
  
  if(ZeroElMass) {
    // Set jpar,Ve,Ajpar neglecting the electron inertia term
    //jpar = ((Te0*Grad_par(Ni, CELL_YLOW)) - (Ni0*Grad_par(phi, CELL_YLOW)))/(fmei*0.51*nu);
    //jpar = ((Tet*Grad_par_LtoC(Ni)) - (Nit*Grad_par_LtoC(phi)))/(fmei*0.51*nu);
    jpar = ((Te0*Grad_par_LtoC(Ni)) - (Ni0*Grad_par_LtoC(phi)))/(fmei*0.51*nu);
  
    if (zlowpass != 0)
      jpar = lowPass(jpar,nzpass);

    // Set boundary conditions on jpar (in BOUT.inp)
    jpar.applyBoundary();
    
    // Need to communicate jpar
    mesh->communicate(jpar);

    //Ve = Vi - jpar/Ni0;
    Ve = jpar/Ni0;
    Ajpar = Ve;
  }else {
    Ve = Ajpar + Apar;
    jpar = Ni0*(Vi - Ve);
  }
 
  //Flux                                                                             
  if(transport) {
    vEB = (B0^Grad(phi))/(mesh->Bxy^2);
    Gamma = vEB*Ni;
  }
  

  // DENSITY EQUATION

  ddt(Ni) = 0.0;
  if(evolve_ni) {
 
    ddt(Ni) -=  vE_Grad(Ni0, phi);// + vE_Grad(Ni, phi0);// + vE_Grad(Ni, phi);
    //ddt(Ni) -= Vpar_Grad_par(Ve, Ni0) + Vpar_Grad_par(Ve0, Ni);// + Vpar_Grad_par(Ve, Ni);

    if (nonlinear)
       ddt(Ni) -= vE_Grad(Ni, phi);
    
    //ddt(Ni) += Ni0*Div_par_CtoL(Ve);// + Ni*Div_par_CtoL(Ve0);// + Ni*Div_par(Vi);
    // ddt(Ni) -= Ni0*Div_par(Ve) + Ni*Div_par(Ve0);
    if (haswak)
      ddt(Ni) +=  Grad_par_CtoL(jpar);
    //ddt(Ni) += Grad_par(jpar);
    //ddt(Ni) += (2.0)*V_dot_Grad(b0xcv, pe);
    /*
    ddt(Ni) -= (2.0)*(Ni0*V_dot_Grad(b0xcv, phi) + Ni*V_dot_Grad(b0xcv, phi0));
    */
    
    //ddt(Ni) += .001*(1.0/(mesh->ngz)) *Laplacian(Ni);
    //if(minusDC) 
     // REMOVE TOROIDAL AVERAGE DENSITY
    if (noDC)
      Ni -= Ni.DC();
    
    if (zlowpass != 0) {
      ddt(Ni) = lowPass(ddt(Ni),nzpass);
      Ni = lowPass(Ni,nzpass);
    }
    if (par_damp)
      ddt(Ni) = lowPass_Y(ddt(Ni),1);
    


    //ddt(Ni) = smooth_y(ddt(Ni));
  }
  
  // ION VELOCITY

  ddt(Vi) = 0.0;
  if(evolve_vi) {
    ddt(Vi) -= vE_Grad(Vi0, phi) + vE_Grad(Vi, phi0);
    if (nonlinear)
      ddt(Vi) -= vE_Grad(Vi, phi);
    
    ddt(Vi) -= Vpar_Grad_par(Vi0, Vi) + Vpar_Grad_par(Vi, Vi0);
     if (nonlinear)
      ddt(Vi) -= vE_Grad(Vi,Vi);
     

    ddt(Vi) -= Grad_par(pei)/Ni0;
     /*ddt(Vi) -= Vpar_Grad_par(Ve0, Vi) + Vpar_Grad_par(Ve, Vi0);// + Vpar_Grad_par(Vi, Vi);
   
    */
    //ddt(Vi) += .001*(1.0/(mesh->ngz)) *Laplacian(Vi);
    if(noDC) 
      ddt(Vi) -= ddt(Vi).DC();
 
    if (zlowpass !=0 ) {
      ddt(Vi) = lowPass(ddt(Vi),nzpass);
      Vi = lowPass(Vi,nzpass);
    }
    if (par_damp)
      ddt(Ni) = lowPass_Y(ddt(Ni),1);
    
  }

  // ELECTRON TEMPERATURE

  ddt(Te) = 0.0;
  if(evolve_te) {
    ddt(Te) -= vE_Grad(Te0, phi) + vE_Grad(Te, phi0) + vE_Grad(Te, phi);
    ddt(Te) -= Vpar_Grad_par(Ve, Te0) + Vpar_Grad_par(Ve0, Te) + Vpar_Grad_par(Ve, Te);
    ddt(Te) += 1.333*Te0*( V_dot_Grad(b0xcv, pe)/Ni0 - V_dot_Grad(b0xcv, phi) );
    ddt(Te) += 3.333*Te0*V_dot_Grad(b0xcv, Te);
    ddt(Te) += (0.6666667/Ni0)*Div_par_K_Grad_par(kapa_Te, Te);
  }

  // ION TEMPERATURE

  ddt(Ti) = 0.0;
  if(evolve_ti) {
    ddt(Ti) -= vE_Grad(Ti0, phi) + vE_Grad(Ti, phi0) + vE_Grad(Ti, phi);
    ddt(Ti) -= Vpar_Grad_par(Vi, Ti0) + Vpar_Grad_par(Vi0, Ti) + Vpar_Grad_par(Vi, Ti);
    ddt(Ti) += 1.333*( Ti0*V_dot_Grad(b0xcv, pe)/Ni0 - Ti*V_dot_Grad(b0xcv, phi) );
    ddt(Ti) -= 3.333*Ti0*V_dot_Grad(b0xcv, Ti);
    ddt(Ti) += (0.6666667/Ni0)*Div_par_K_Grad_par(kapa_Ti, Ti);
  }

  // VORTICITY

  ddt(rho) = 0.0;
  
  if(evolve_rho) {
      
    // ddt(rho) -= vE_Grad(rho0, phi);// + vE_Grad(rho, phi0);//+ vE_Grad(rho, phi);
    //ddt(rho) -= vE_Grad(Ni0, phi);
    //ddt(rho) += mesh->Bxy*mesh->Bxy*Div_par(jpar, CELL_CENTRE);
    
    ddt(rho) += mesh->Bxy*mesh->Bxy*Grad_par_CtoL(jpar); 
    //ddt(rho) += Grad_par_CtoL(jpar); //for 2D the above line is equivalent
    // ddt(rho) += 2.0*mesh->Bxy*V_dot_Grad(b0xcv, pei);
    
    //ddt(rho) -= Vpar_Grad_par(Vi, rho0) + Vpar_Grad_par(Vi0, rho);// + Vpar_Grad_par(Vi, rho);
    
    
    if (nonlinear)
       ddt(rho) -= vE_Grad(rho, phi);  
    
    
    /*
    for(int jx=MXG;jx<mesh->ngx-MXG;jx++) {
      for(int jy=MYG;jy<mesh->ngy-MYG;jy++) {
	for(int jz=0;jz<mesh->ngz;jz++) {
	  ddt(rho)[jx][jy][jz] = Bxy[jx][jy]*Bxy[jx][jy] * (jpar[jx][jy+1][jz] - jpar[jx][jy][jz]) / (dy[jx][jy] * sqrt(mesh->g_22[jx][jy]));
	}
      }
    }
    */
    //output.write("mesh->dz)^2: %e \n",(mesh->dz)^2);
    //ddt(rho) += 10*(1.0/(mesh->ngz)) * (1.0/(mesh->ngz))*Laplacian(rho);
    //ddt(rho) += Laplacian(rho);
    //ddt(rho) += 1.0/100000000.0 *mu_i*Delp2(rho,-1.0); //no smoothing, check the value of mu_i
    //ddt(rho) += 1.0/100000000.0 *Delp2(rho); 
    if(noDC) 
      rho -= rho.DC();
    
    //ddt(rho) += 1e-4 * mu_i * Laplacian(rho);
    if (zlowpass != 0) {
      ddt(rho) = lowPass(ddt(rho),nzpass);
      rho = lowPass(rho,nzpass);
    }
//ddt(rho) = smooth_y(ddt(rho));
    if (par_damp)
      ddt(rho) = lowPass_Y(ddt(rho),1);


  }
  

  // AJPAR
  

  ddt(Ajpar) = 0.0;
  if(evolve_ajpar) {
    //ddt(Ajpar) -= vE_Grad(Ajpar0, phi) + vE_Grad(Ajpar, phi0) + vE_Grad(Ajpar, phi);

    /*
    for(int jx=MXG;jx<mesh->ngx-MXG;jx++) {
      for(int jy=MYG;jy<mesh->ngy-MYG;jy++) {
	for(int jz=0;jz<mesh->ngz;jz++) {
	  ddt(Ajpar)[jx][jy][jz] += (1./fmei) * (phi[jx][jy][jz] - phi[jx][jy-1][jz]) / (dy[jx][jy] * sqrt(mesh->g_22[jx][jy]));
	  ddt(Ajpar)[jx][jy][jz] -= (1./fmei)*(Te0[jx][jy]/Ni0[jx][jy])*(Ni[jx][jy][jz] - Ni[jx][jy-1][jz]) / (dy[jx][jy] * sqrt(mesh->g_22[jx][jy]));
	}
      }
    }
    */
    //ddt(Ajpar) -= (1./fmei)*1.71*Grad_par(Te);

    // ddt(Ajpar) += (1./fmei)*Grad_par(phi, CELL_YLOW);
    // ddt(Ajpar) -= (1./fmei)*(Te0/Ni0)*Grad_par(Ni, CELL_YLOW);
    // ddt(Ajpar) += 0.51*interp_to(nu, CELL_YLOW)*jpar/Ni0;

    ddt(Ajpar) += (1./fmei)*Grad_par_LtoC(phi);
    ddt(Ajpar) -= (1./fmei)*(Te0/Ni0)*Grad_par_LtoC(Ni);
    ddt(Ajpar) += 0.51*interp_to(nu, CELL_YLOW)*jpar/Ni0;
    
  }

  return(0);
}

/*******************************************************************************
 *                       FAST LINEAR FIELD SOLVERS
 *******************************************************************************/

#include <invert_laplace.hxx>

// Performs inversion of rho (r) to get phi (p)
int solve_phi_tridag(Field3D &r, Field3D &p, int flags)
{
  //output.write("Solving phi: %e, %e -> %e\n", max(abs(r)), min(Ni0), max(abs(r/Ni0)));

  if(invert_laplace(r/Ni0, p, flags, NULL)) {
    return 1;
  }

  //Field3D pertPi = Ti*Ni0 + Ni*Ti0;
  //p -= pertPi/Ni0;
  return(0);
}

int solve_apar_tridag(Field3D &aj, Field3D &ap, int flags)
{
  static Field2D a;
  static int set = 0;

  if(set == 0) {
    // calculate a
    a = (-0.5*beta_p/fmei)*Ni0;
    set = 1;
  }

  if(invert_laplace(a*(Vi - aj), ap, flags, &a))
    return 1;

  return(0);
}



