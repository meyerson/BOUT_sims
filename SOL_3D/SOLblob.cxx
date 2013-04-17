/*******************************************************************
 * Bare Bones Bout with Python
 *
 * D. SMeyerson 2013
 *******************************************************************/

#include <bout.hxx>
#include <boutmain.hxx>
#include <mpi.h>


#include <typeinfo>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#include <derivs.hxx>
#include <initialprofiles.hxx>
#include <invert_laplace.hxx>
#include <invert_parderiv.hxx>
#include <invert_laplace_gmres.hxx>
#include <boutexception.hxx>

//#include "./StandardAlpha.h"
//#include <python2.6/Python.h>

// Evolving variables 
Field3D u, n; //vorticity, density

//derived variables
Field3D phi,brkt;


//other fields
Field3D test1, test2, ReyN, ReyU;
Field2D n0;

//Constrained 
Field3D C_phi;

//other params
BoutReal nu, mu,gam, beta,alpha_c;

Field3D alpha;
//solver options
bool use_jacobian, use_precon;

//experimental
bool use_constraint;
bool chaosalpha;

int MZ;
// 2D initial profiles
Field2D Ni0, Ti0, Te0, Vi0, phi0, Ve0, rho0, Ajpar0,gradNi0;
Vector2D b0xcv, b0,B0,B0xcv,b0xcv_cn; /// for curvature terms

Field2D pei0,pe0;
// Metric coefficients
Field2D Rxy, Bpxy, Btxy, hthe,I;

// parameters
BoutReal Te_x, Ti_x, Ni_x, Vi_x, bmag, rho_s, fmei, AA, ZZ;
BoutReal lambda_ei, lambda_ii;
BoutReal nu_hat, mui_hat, wci, nueix, nuiix;
BoutReal beta_p;

// settings
bool estatic, ZeroElMass; // Switch for electrostatic operation (true = no Apar)
BoutReal zeff, nu_perp;
bool evolve_rho, evolve_te, evolve_ni, evolve_ajpar, evolve_vi, evolve_ti;
BoutReal ShearFactor;

int phi_flags, apar_flags; // Inversion flags



FieldGroup comms; // Group of variables for communications

const Field3D mybracket(const Field3D &phi, const Field3D &A);
int jacobian(BoutReal t); // Jacobian-vector multiply
int precon(BoutReal t, BoutReal cj, BoutReal delta); // Preconditioner


//BoutReal alphamap(BoutReal x,BoutReal z);
BoutReal alphamap(double x, double Lx, double y,double Ly,
		  double k=.50,double q0=5.0,double R=100,
		  int max_orbit =1000,double period=10.0);
//int alphamapPy();
//int precon_phi(BoutReal t, BoutReal cj, BoutReal delta);
//int jacobian_constrain(BoutReal t); // Jacobian-vector multiply

int physics_init(bool restarting)
{
  
  GRID_LOAD(Ni0);
  GRID_LOAD(Ti0);
  GRID_LOAD(Te0);
  GRID_LOAD(phi0);

  // Load magnetic curvature term
  b0xcv.covariant = false; // Read contravariant components
  mesh->get(b0xcv, "bxcv"); // b0xkappa terms

  // Load metrics
  GRID_LOAD(Rxy);
  GRID_LOAD(Bpxy);
  GRID_LOAD(Btxy);
  GRID_LOAD(hthe);
  mesh->get(mesh->dx,   "dx"); //dpsi really
  mesh->get(I,    "sinty");
  //mesh->get(mesh->zShift, "qinty");


  // // Load normalisation values
  GRID_LOAD(Te_x);
  GRID_LOAD(Ti_x);
  GRID_LOAD(Ni_x);
  GRID_LOAD(bmag);

  //Ni_x *= 1.0e14;
  //bmag *= 1.0e4;
 
  /*************** READ OPTIONS *************************/
  Options *globaloptions = Options::getRoot();
  Options *options = globaloptions->getSection("physics");
  Options *solveropts = globaloptions->getSection("solver");


  OPTION(options, AA, 2.0);
  OPTION(options, ZZ, 1.0);

  OPTION(options, phi_flags, 0);
  //OPTION(options, alpha,3e-5);
  OPTION(options, nu, 2e-3);
  //OPTION(options, mu, 0.040);
  OPTION(options,chaosalpha,false);
  OPTION(options, mu, 2e-3);
  OPTION(options, gam, 1e1);
  OPTION(options, beta, 6e-4);

  OPTION(globaloptions,MZ,33);

  OPTION(solveropts,use_precon,false);
  OPTION(solveropts,use_jacobian,true);
  OPTION(solveropts,use_constraint,false);


 /************** CALCULATE PARAMETERS *****************/
  
  Ni_x *= 1.0e14; //old uedge idl routine spits out in units of 10^20 m^-3, we want per cm^3
  Ni0 *= 1.0e14;
  
  bmag *= 1.0e4; //NRL formulas are in cgs the grid produces b in SI

  rho_s = 1.02*sqrt(AA*Te_x)/ZZ/bmag; //in METERS!!
  fmei  = 1./1836.2/AA;

  lambda_ei = 24.-log(sqrt(Ni_x)/Te_x);
  lambda_ii = 23.-log(ZZ*ZZ*ZZ*sqrt(2.*Ni_x)/pow(Ti_x, 1.5));
  wci       = 9.58e3*ZZ*bmag/AA;
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
  phi.setLocation(CELL_CENTRE);
  //jpar.setLocation(CELL_YLOW);

  /************** NORMALISE QUANTITIES *****************/
  Ni0 /= Ni_x;
  Ti0 /= Te_x;
  Te0 /= Te_x;
 
  //again we already changed bmag to cgs, but many other uedge generrated quanitites are SI
  // Normalise curvature term
  b0xcv.x /= (bmag/1e4);
  b0xcv.y *= rho_s*rho_s;
  b0xcv.z *= rho_s*rho_s;
  
  // Normalise geometry, all original lengths are in meters
  Rxy /= rho_s;
  hthe /= rho_s;
  I *= rho_s*rho_s*(bmag/1e4)*ShearFactor;
  output.write("mesh->dx = %e\n", mesh->dx[0]);
  mesh->dx /= rho_s*rho_s*(bmag/1e4);
  //mesh->dy /= rho_s; //nope, a dy is already unitless

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

  B0xcv = b0 ^ V_dot_Grad(b0,b0);
  dump.add(B0xcv.x,"B0xcvx",0);
  dump.add(B0xcv.y,"B0xcvy",0);
  dump.add(B0xcv.z,"B0xcvz",0);
  
  //dump the original b0xcv 
  dump.add(b0xcv.x,"b0xcvx",0);
  dump.add(b0xcv.y,"b0xcvy",0);
  dump.add(b0xcv.z,"b0xcvz",0);

  bout_solve(u, "u");
  comms.add(u);
  //phi = invert_laplace(u, phi_flags);
  static Field2D A = 0.0;
  static Field2D C = 1e-12;
  static Field2D D = 1.0;
  
  phi = invert_laplace(u, phi_flags,&A,&C,&D);
  //Laplacian *lap = Laplacian::create();
  
  bout_solve(n, "n");
  comms.add(n);

  //brute force way to set alpha
  
  if (chaosalpha){
    alpha.allocate();
    BoutReal ***a = alpha.getData();

    for(int jz=0;jz<mesh->ngz-1;jz++) {
      for(int jx=0;jx<mesh->ngx;jx++)
	for(int jy=0;jy<mesh->ngy;jy++){
	  a[jx][jy][jz]=alphamap(mesh->GlobalX(jx),1.0,mesh->dz*jz,mesh->zlength);
	}
    }
    alpha = (2.0*.1)/alpha;
    dump.add(alpha,"alpha",0);
  } else{
    OPTION(options, alpha_c,3e-5);
    alpha = alpha_c;
  }
  
  dump.add(brkt,"brkt",1);
  dump.add(test1,"test1",1);
  dump.add(test2,"test2",1);
  dump.add(ReyN,"ReyN",1);
  dump.add(ReyU,"ReyU",1);


    // Add any other variables to be dumped to file
  //dump.add(phi,  "phi",  1);
  // dump.add(Apar, "Apar", 1);
  // dump.add(jpar, "jpar", 1);
  
  // dump.add(ddt(Ni),"ddtNi",1);
  // dump.add(ddt(rho),"ddtrho",1);

  dump.add(Ni0, "Ni0", 0);
  
  // dump.add(gradNi0,"gradNi0",0);
  // dump.add(DXRxy,"DXRxy",0);
  dump.add(Te0, "Te0", 0);
  dump.add(Ti0, "Ti0", 0);

  dump.add(Te_x,  "Te_x", 0);
  dump.add(Ti_x,  "Ti_x", 0);
  dump.add(Ni_x,  "Ni_x", 0);
  dump.add(rho_s, "rho_s", 0);
  dump.add(wci,   "wci", 0);
  
  //dump.add(mesh->g11,   "g11", 0);

  

  if (use_constraint){
    //solver->setPrecon(precon_phi);
    //solver->setJacobian(jacobian_constrain);
    phi.setBoundary("phi");
    //bout_solve(phi,"phi");
  }else
    dump.add(phi,"phi",1);

  comms.add(phi); //super duper important 

  if (use_jacobian)
    solver->setJacobian(jacobian);

  // if (use_precon)
  //   solver->setPrecon(precon);
    
  output.write("use jacobian %i \n",use_jacobian);
  output.write("use precon %i \n",use_precon);
  output<< "mesh->hthe0 "<<hthe0<<endl;
  output.write("DONE WITH PHYSICS_INIT\n");

  n0 = 1.0;
  return 0;
}

#define bracket3D(f, g) ( b0xGrad_dot_Grad(f, g) / mesh->Bxy)
#define vE_Grad(f, p) ( b0xGrad_dot_Grad(p, f) / mesh->Bxy )
#define LapXZ(f)(mesh->g11*D2DX2(f) + mesh->g33*D2DZ2(f))

int physics_run(BoutReal t)
{
  // Run communications
  mesh->communicate(comms);
  //phi = invert_laplace(u, phi_flags);
  
  static Field2D A = 0.0;
  static Field2D C = 1e-24;
  static Field2D D = 1.0;
  
  phi = invert_laplace(u, phi_flags,&A,&C,&D);

  phi.applyBoundary("dirichlet");
  // Density
  //f = lowPass(f,1);
  //f = lowPass(g,1);
  mesh->communicate(comms);
  //mesh->communicate(phi);
  ddt(u)=0;
  ddt(n)=0;
  
  ReyU = bracket3D(phi,u)/(nu*LapXZ(u)+1e-5);

 
  ddt(u) -= bracket3D(phi,u);
  ddt(u) += alpha * phi;
  ddt(u) += nu * LapXZ(u);
  //ddt(u) -= beta* DDZ(n+n0)/(n+n0);
  ddt(u) += 2.0*mesh->Bxy*V_dot_Grad(b0xcv, n);
  //ddt(u) += 2.0*mesh->Bxy*b0xcv*Grad(n);
  
  ReyN = bracket3D(phi,n)/(mu * LapXZ(n)+1e-5);
  
  ddt(n)  -= bracket3D(phi,n+n0);
  ddt(n) += mu * LapXZ(n+n0);

 
  return 0;
}


const Field3D mybracket(const Field3D &phi, const Field3D &A)
{
  Field3D dpdx, dpdy, dpdz;
  Field3D vx, vy, vz;
  Field3D result;

  //output.write("mesh->Bxy = %e\n", (mesh->J*sqrt(mesh->g_22))[2][2]);

#ifdef CHECK
  int msg_pos = msg_stack.push("b0xGrad_dot_Grad( Field3D , Field3D )");
#endif

  // Calculate phi derivatives
  #pragma omp parallel sections
  {
    #pragma omp section
    dpdx = DDX(phi); 
    
    #pragma omp section
    dpdy = DDY(phi);
    
    #pragma omp section
    dpdz = 0;
  }
  
  // Calculate advection velocity
  #pragma omp parallel sections
  {

    #pragma omp section
    vx = mesh->g_23*dpdz - mesh->g_33*dpdy; 
    //vx = mesh->g_22*dpdz - mesh->g_23*dpdy;
    
    #pragma omp section
    vy = mesh->g_33*dpdx - mesh->g_13*dpdz;
      //vy = mesh->g_23*dpdx - mesh->g_12*dpdz;
    
    #pragma omp section
    //vz = mesh->g_12*dpdy - mesh->g_22*dpdx;
    vz =0;
  }


  // Upwind A using these velocities
  
  Field3D ry, rz;
  #pragma omp parallel sections
  {
    #pragma omp section
    result = VDDX(vx, A);
    
    #pragma omp section
    ry = VDDY(vy, A);
    
    #pragma omp section
    rz = VDDZ(vz, A);
  }
  //output.write("mesh->g_22: %g  \n" ,vx[4][4]);
  //result = (ry + rz); //usually  convention
  result = (result + ry) / (mesh->J*sqrt(mesh->g_22));

#ifdef TRACK
  result.name = "b0xGrad_dot_Grad("+phi.name+","+A.name+")";
#endif
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif
  return result;
}



/* computes Jv, where ddt() is holding v and the system state holds Jv */ 
int jacobian(BoutReal t) {
  mesh->communicate(ddt(u),ddt(n));
  
  static Field2D A = 0.0;
  static Field2D C = 1e-12;
  static Field2D D = 1.0;
  
  ddt(phi) = invert_laplace(ddt(u), phi_flags,&A,&C,&D);
  //ddt(phi) = invert_laplace(ddt(u), phi_flags); 

  mesh->communicate(ddt(phi));

  u=0;
  n=0;

  //u -= mybracket(ddt(phi),ddt(u));
  u += bracket3D(ddt(phi),ddt(u));
  //ddt(u) += alpha * phi;
  u += nu * Delp2(ddt(u));
  //ddt(u) -= beta * DDY(n)/n; 
  //ddt(u) -= beta* Grad_par(n)/n; 
  //u -= Grad_par(ddt(n)); 
  //ddt(u).applyBoundary("dirichlet");
  u -= beta* DDZ(n); 
  //mesh->communicate(comms); no don't do this here
  //.applyBoundary();
  //brkt = VDDY(DDY(phi), n) +  VDDZ(DDZ(phi), n) ;
 
  n += bracket3D(ddt(phi),ddt(n));
  //n -= mybracket(ddt(phi),ddt(n));
  n += mu * Delp2(ddt(n));
  //n -= alpha* n;
  
  n.applyBoundary();
  u.applyBoundary();
  return 0;
}

/* computes P^-1 r, where ddt() holds (-r) and the system state hold P^-1 r*/

int precon(BoutReal t, BoutReal gamma, BoutReal delta) {
  // mesh->communicate(rhscomms);
  mesh->communicate(ddt(n),ddt(u));

  n= 0;
  u =0;

  n += ddt(n);
  // mesh->communicate(n);
  // Ni -= (mesh->Bxy*mesh->Bxy*ddt(Ni) - ddt(rho))/(mesh->Bxy*mesh->Bxy);

  u += ddt(u);
  u -= gamma * Grad_par(ddt(n)); 
 
return 0;
 
}

BoutReal alphamap(double x, double Lx,double y,double Ly,
		  double k,double q0 ,double R,int max_orbit,double period){ 
  
  //rescale
  // x = jx*dx+x0;
  // u = jy*dy+y0;

  
  // x = (x-x0)*(2*M_PI/Lx)-M_PI;
  // y = (y-y0)*(2*M_PI/Ly);
  
  x = x*(2*M_PI/Lx)-M_PI;
  y = y*(2*M_PI/Ly);
  //output << "[" << x  << "], "<<endl;
  int count = 0;
  bool hit_divert = false;
  bool inSOL;
  double q;

  double L = 0.0;
  
  while(count < max_orbit and not hit_divert){
    inSOL = x>0; //are we in SOL now?
    
    //update the field line
    x = x+k*sin(period * y);
    y = fmod((x+y),(2*M_PI)); //one can argue that x = 2*M_PI*fmod(q(R),1)

    q =q0 +x/(2*M_PI);

    L = L+q*R*2*M_PI; //one can argue that for x>0, it should be L+q*R*M_PI

    hit_divert = (inSOL and x>0) or (x>M_PI); //did the field line line stay in SOL?
    //cout <<count<<" L:" << L<<endl;
    count++;
    }
  
  return L;

}
