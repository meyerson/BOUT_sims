/*******************************************************************
 * Bare Bones Bout with Python
 *
 * D. Meyerson 2013
 *******************************************************************/

#include <bout.hxx>
#include <boutmain.hxx>
#include <mpi.h>
#include <math.h>

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
#include <field_factory.hxx>

//#include "./StandardAlpha.h"
//#include <python2.6/Python.h>

// Evolving variables 
Field3D u, n,n_prev, u_prev,Te, Te_prev; //vorticity, density

//derived variables
Field3D phi,brkt;
int phi_flags;

//other fields
Field3D test1, test2, ReyN, ReyU;
Field2D n0;

//Constrained 
Field3D C_phi;

//other params
BoutReal nu, mu,gam, beta,alpha_c, eps,fmei,kpar,AA,ZZ,a_dw;
BoutReal TIMESTEP;
Field3D alpha, temp,edgefld,alpha_s, alpha_j,source,sink,nave,uave,uDC,nDC,n_rms,u_rms;

Field3D alpha_mask,div_jpar;
BoutReal Te0;
//solver options
bool use_jacobian, use_precon;
bool evolve_te;
//experimental

bool withsource,wave_bc,diff_bc,withsink;
bool use_constraint;
string chaosalpha;
bool inc_jpar;
bool log_n;
int smoother_a;
BoutReal n_sol;
int max_orbit;

int m;

int MZ;

FieldGroup comms; // Group of variables for communications

const Field3D mybracket(const Field3D &phi, const Field3D &A);
int jacobian(BoutReal t); // Jacobian-vector multiply
int precon(BoutReal t, BoutReal cj, BoutReal delta); // Preconditioner


//BoutReal alphamap(BoutReal x,BoutReal z);
BoutReal alphamap(double x, double Lx, double y,double Ly,
		  double k=1.20,double q0=3.0,double R=100,
		  int max_orbit =4000,double period=1.0,
		  bool count_turn = 0);

BoutReal Ullmann(double x, double Lx, double y,double Ly, double x_sol, BoutReal eps, double m);
BoutReal Newton_root(double x_in,double y_in,double b = 50.0, double C=.01, double m = 3.0);

const Field3D smooth_xz(const Field3D &f); 
//const Field3D GT(const Field3D &f,Field3D &g);

//int alphamapPy();
//int precon_phi(BoutReal t, BoutReal cj, BoutReal delta);
//int jacobian_constrain(BoutReal t); // Jacobian-vector multiply

int physics_init(bool restarting)
{
  


  Options *globaloptions = Options::getRoot();
  Options *options = globaloptions->getSection("physics");
  Options *solveropts = globaloptions->getSection("solver");
  Options *laplaceopts = globaloptions->getSection("laplaceopts");

  OPTION(laplaceopts, phi_flags, 0);
  //OPTION(options, alpha,3e-5);
  OPTION(options, nu, 2e-3);
  //OPTION(options, mu, 0.040);
  OPTION(options,chaosalpha,"chaos");
  OPTION(options, mu, 2e-3);
  OPTION(options, gam, 1e1);
  OPTION(options, beta, 6e-4);
  OPTION(options, inc_jpar,false);
  OPTION(options, smoother_a,0);
  OPTION(options, eps, 2e-1);
  OPTION(options, m, 3);
  OPTION(options, max_orbit, 100);

  OPTION(globaloptions,MZ,33);

  (globaloptions->getSection("te"))->get("evolve", evolve_te,   false);

  OPTION(solveropts,use_precon,false);
  OPTION(solveropts,use_jacobian,true);
  OPTION(solveropts,use_constraint,false);

  OPTION(options,withsource,false);
  OPTION(options,withsink,true);
  OPTION(options,wave_bc,true);
  OPTION(options,diff_bc,false);
  OPTION(options,log_n,true);

  OPTION(options,Te0,1e0);
  //OPTION(options,n0,1e0);
  OPTION(options, AA, 2.0); //deutrium ?
  OPTION(options, ZZ, 1.0); //singly ionized
  OPTION(options, alpha_c,3e-5);
  OPTION(globaloptions,TIMESTEP,1.0);

  bout_solve(u, "u");
  comms.add(u);
  //phi = invert_laplace(u, phi_flags);
  static Field2D A = 0.0;
  static Field2D C = 1e-12;
  static Field2D D = 1.0;
  
  phi = invert_laplace(u, phi_flags,&A,&C,&D);
  //phi = u*0.0;
  //Laplacian *lap = Laplacian::create();
  
  bout_solve(n, "n");
  comms.add(n);
 
  phi = phi + n.DC();
  FieldFactory f(mesh);
  if(withsource){
    //initial_profile("source", v);
    source = f.create3D("gauss(x-0.0,0.02)");
    //source = f.create3D("h(.05-x)");

    dump.add(source,"source",0);
    
  }

  if(withsink){
    //sink = 1.0 - f.create3D("h(x-.9)");
    sink = f.create3D("gauss(x-1.0,.02)");
    dump.add(sink,"sink",0);
    
  }

  if(inc_jpar){
    dump.add(div_jpar,"div_jpar",1);
    div_jpar.setBoundary("phi");
    if (evolve_te == false)
      Te = Te0;
    // comms.add(div_jpar);
  }
  
  if(evolve_te) {
    bout_solve(Te, "Te");
    comms.add(Te);
    //rhscomms.add(ddt(Te));
    output.write("te\n");
    Te_prev = Te;
    
  }

  /************** CALCULATE PARAMETERS *****************/

  //rho_s = 1.02*sqrt(AA*Te_x)/ZZ/bmag;
  //AA is the mass of the ion species in ion masses , default is 2
  // ZZ is the charge state , default is 1
  // zeff
 
  // wci       = (1.0)*9.58e3*ZZ*bmag/AA;
  // lambda_ei = 24.-log(sqrt(Ni_x)/Te_x);
  // nueix     = 2.91e-6*Ni_x*lambda_ei/pow(Te_x, 1.5);
  // nu_hat    = nue
  //brute force way to set alpha
  

  fmei  = 1./1836.2/AA;
  //kpar = alpha_c;  //simplest possible 
  //a_dw = pow(kpar,2.0)/(fmei*.51*.1);
  a_dw = .20;

  alpha.allocate();
  alpha_j.allocate();
  alpha_mask.allocate();

  BoutReal ***a = alpha.getData();
  BoutReal ***a_j = alpha_j.getData();
  BoutReal ***a_m = alpha_mask.getData();

  BoutReal edge[mesh->ngz];
  BoutReal zoomfactor = 3.0;
  BoutReal lowR = .55;
  BoutReal Lxz = 0;
  BoutReal x_sol = .3;
  BoutReal rho_s = .5;

  for(int jz=0;jz<mesh->ngz;jz++) 
    for(int jx=0;jx<mesh->ngx;jx++){
      Lxz = Ullmann(mesh->GlobalX(jx),1.0,mesh->dz*jz,mesh->zlength,x_sol,eps,m);
      
      // for(int jy=0;jy<mesh->ngy;jy++){
      // 	a[jx][jy][jz]=(Lxz>0)*(1.0/Lxz);
      // 	a_m[jx][jy][jz]=(Lxz<0);
      // 	if (mesh->firstX()){
      // 	  a_m[jx][jy][jz]=0;
      // 	  a[jx][jy][jz]= 0;
      // 	}
      // 	a_j[jx][jy][jz]= alpha_c*double(mesh->GlobalX(jx) > x_sol);
      // }

      for(int jy=0;jy<mesh->ngy;jy++){
	if  ("jump" == chaosalpha) {
	  a[jx][jy][jz]=alpha_c*double(mesh->GlobalX(jx) > x_sol);
	  a_m[jx][jy][jz]= double(a[jx][jy][jz] == 0.); //double(mesh->GlobalX(jx) <= x_sol);
	}
	else {
	  a[jx][jy][jz]=(Lxz>0)*(1.0/Lxz);
	  a_m[jx][jy][jz]=(Lxz<0);
	}
	
	if (mesh->firstX()){
      	  a_m[0][jy][jz]=0;
      	  a[0][jy][jz]= 0;
	}


      }
    }
    


  //alpha = rho_s/alpha;
  alpha = rho_s * alpha;
  //alpha = double(alpha > 0)*alpha;
    //alpha = alpha * alpha_c/alpha.max(1);

  alpha_s = lowPass(alpha,0);
  
  //normalize
  alpha = alpha * alpha_c/alpha_s.max(1);
  output.write("smoothing %g \n",smoother_a);
  for (int a_i=smoother_a;a_i>=0;a_i--){
    alpha = smooth_xz(alpha);
    output.write("smoothing %i \n",a_i);
  }

  alpha_s = alpha_s * alpha_c/alpha_s.max(1);

  
  // if ("jump" == chaosalpha){
  //   alpha_mask = ;
  //   alpha = alpha_j; }
  if ("smooth" == chaosalpha){
    alpha = alpha_s;
    alpha_mask = lowPass(alpha_mask,0);
  }
  
   
  dump.add(alpha,"alpha",0);
  dump.add(alpha_mask,"alpha_mask",0);
  dump.add(alpha_s,"alpha_smooth",0);
  dump.add(eps,"eps",0);
  dump.add(m,"m",0);
  dump.add(a_dw,"a_dw",0);
  

  //dump.add(brkt,"brkt",1);
  //dump.add(test1,"test1",1);
  //dump.add(test2,"test2",1);
  dump.add(ReyN,"ReyN",1);
  dump.add(ReyU,"ReyU",1);
  dump.add(nave,"nave",0);
  dump.add(uave,"uave",0);

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
  output.write("DONE WITH PHYSICS_INIT\n");

  n0 = 1e-4;
  n_prev =n;
  u_prev = u;
  nave = n;
  uave = u;


  if (log_n){
    n = log(n+n0);
    n0 = log(n0);
    n_prev =n;
  }



  ddt(n).setBoundary ("ddt[n]") ;
  phi.setBoundary("phi");
  ddt(u).setBoundary ("ddt[u]") ;
 
  mesh->communicate(comms);
  return 0;
}

#define bracket3D(f, g) ( b0xGrad_dot_Grad(f, g) )
#define LapXZ(f)(mesh->g11*D2DX2(f) + mesh->g33*D2DZ2(f))
//#define LapXZ(f)(D2DX2(f));// + mesh->g33*D2DZ2(f))

int physics_run(BoutReal t)
{
  // Run communications
  mesh->communicate(comms);
  //phi = invert_laplace(u, phi_flags);
  

  //mesh->communicate(comms);

  //n.applyBoundary("dirichlet(4.6)");

  if (wave_bc){

    uDC = u.DC();
    nDC = n.DC();

    Field3D nprevDC = n_prev.DC();    
    Field3D uprevDC = u_prev.DC();


    if (mesh->firstX())
      for(int i=1;i>=0;i--)
    	for(int j =0;j< mesh->ngy;j++)
    	  for(int k=0;k < mesh->ngz; k++){
	    n[i][j][k] =(.5*n[i+1][j][k] + n_prev[i][j][k])/(1.0+.5);
	    //n[i][j][k] =(.5*(nDC[i+1][j][k]) + (nprevDC)[i][j][k])/(1.0+.5);
	    //u[i][j][k] =(.5*u[i+1][j][k] + u[i][j][k])/(1.0+.5);
	    //n[i][j][k] = (n.DC)[i][j][k];
	    //u[i][j][k] = (u.DC)[i][j][k];
	  }
    

    
    
  } else{
    //n.applyBoundary();
    
    if (evolve_te)
      Te.applyBoundary();
  }
  
  static Field2D A = 0.0;
  static Field2D C = 1e-24;
  static Field2D D = 1.0;
  FieldFactory f(mesh);

  u.applyBoundary(); //BIG speed up
  phi = invert_laplace(u, phi_flags,&A,&C,&D);
  //phi = u*0.0;
  //phi.applyBoundary("dirichlet"); //SLLOW and unstable
  //phi.applyBoundary();

  
  mesh->communicate(comms);
  //mesh->communicate(phi);
  ddt(u)=0;
  ddt(n)=0;
 


  
  //ReyU = bracket3D(phi,u)/(nu*LapXZ(u)+1e-5);

 
  ddt(u) -= bracket3D(phi,u);
  ddt(u) += alpha * phi;
  ddt(u) += nu * LapXZ(u);



  if (log_n){
    ddt(u) += beta* DDZ(n);
 
    ddt(n) -= bracket3D(phi,n);
   
    ddt(n) += mu * (LapXZ(n) + Grad(n)*Grad(n)) ;
  
    ddt(n) -= alpha;	    
   
  } else {
    ddt(u) += beta* DDZ(n+n0)/(n+n0);
   
    //output.write ("no log_n \n");
    // ReyN = bracket3D(phi,n)/(mu * LapXZ(n)+1e-5);
    
    ddt(n) -= bracket3D(phi,n);
    ddt(n) += mu * (LapXZ(n)) ;
    ddt(n) -= alpha* n;
  }
  
  //ddt(n) = n*0;
  //ddt(u) = u*0;
  
  
 
  if(withsource ){
    if (log_n)
      //ddt(n) += (5.0e0 * alpha_c * source)/exp(n);
      ddt(n) += (5.0e-1 * alpha_c * source)/exp(n);
    else
      ddt(n) += (1.0e0 * alpha_c * source);
  }

  if(withsink){
    BoutReal target_val;
    //target_profile = smooth_x(smooth_x(n.DC()));
    
    Field3D region_select = f.create3D("h(.9-x)");
    //Field3D region_select = f.create3D("gauss(x-.95,.02)");// + f.create3D("h(.95-x)") - 1.0  ;

    //target_val = (region_select.DC()*n.DC()).mean(true)/((region_select.DC()).mean(true)); //slow


    target_val = min((n* region_select).DC(),true)-.01;
    //output.write("tarval_val  %g \n",target_val);
    //target_val = -6.0;
    ddt(n) -= (1e1*alpha_c *sink)*(1.0-exp(target_val)/exp(n));
    ddt(u) -= (3e0*alpha_c *sink)*(u - u.DC());


  }
  // if(withsink){
  //   ddt(n) -= (2.0e-2 * n * sink);
  // }

  u.applyBoundary();
  //n.applyBoundary();


  if(inc_jpar){
      // Update non-linear coefficients on the mesh
    //nu      = nu_hat * n/ (Te0^1.5);
    //lambda_ei = 24.-log(sqrt(n0)/Te_x);
    //nueix     = 2.91e-6*Ni_x*lambda_ei/pow(Te_x, 1.5);
    
    //jpar = ((Te0*Grad_par_LtoC(n)) - (n0*Grad_par_LtoC(phi)));///(fmei*0.51*nu);
    //div_jpar = -pow(kpar,2.0)*(lazy_log(n)*Te0 -phi)/(fmei*.51*.1);//*(log(n)*Te - phi)/(fmei*.51*nu);
    
    if (log_n)
      div_jpar = (n*Te0 - phi);//*(log(n)*
    else {
      //phi = phi - lazy_log(n[0][0][0])*Te0;
      div_jpar = (log(abs(n+n0))*Te0 - phi);//*(log(n)

    }
    div_jpar = div_jpar - div_jpar.DC();
    div_jpar.applyBoundary();
    // div_jpar.applyBoundary();
    // //for values where alpha  = min
    
    
    if (log_n)
      {
	ddt(u) -= a_dw*smooth_xz(smooth_xz(smooth_xz(alpha_mask)))*div_jpar/(exp(n).DC());
	ddt(n) -= a_dw*smooth_xz(smooth_xz(smooth_xz(alpha_mask)))*div_jpar/(exp(n).DC());
      }
    else
      {
	ddt(u) -= a_dw*smooth_xz(smooth_xz(smooth_xz(alpha_mask)))*div_jpar/n.DC();
	ddt(n) -= a_dw*smooth_xz(smooth_xz(smooth_xz(alpha_mask)))*div_jpar;
      }
  }

  ddt(Te) = 0.0;
  if(evolve_te) {
    ddt(Te)  -= bracket3D(phi,Te);
    ddt(Te) += (mu/10.) * LapXZ(Te);
    ddt(Te) -= alpha *Te;
    Te_prev = Te;
  }
   

  ddt(n).applyBoundary(); //extermely important 
  ddt(u).applyBoundary();



  mesh->communicate(comms);
  n_prev = n;
  u_prev = u;

  //n_rms = n-n.DC();
  
  nave = nave + n;
  uave = uave + u;
  
  //ddt(u) = uDC;

  return 0;
}
// const Field3D GT(const Field3D &f, const Field3D &g)
// {

// }

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
  // ddt(u).applyBoundary("dirichlet");
  // ddt(n).applyBoundary("dirichlet");
  // ddt(phi).applyBoundary("dirichlet");
  mesh->communicate(ddt(phi),ddt(u),ddt(n));

  u=0;
  n=0;

  //u -= mybracket(ddt(phi),ddt(u));
  //u -= bracket3D(ddt(phi),ddt(u));
  u += alpha * ddt(phi);

  //u += nu * LapXZ(ddt(u));
 
  //ddt(u).applyBoundary("dirichlet");


  
  if (log_n){
    u += beta* DDZ(ddt(n));
    n -= bracket3D(ddt(phi),ddt(n));

    if(inc_jpar){
      div_jpar = (ddt(n) - ddt(phi));
      div_jpar = div_jpar - div_jpar.DC();
      div_jpar.applyBoundary();
      u -= a_dw*smooth_xz(smooth_xz(smooth_xz(alpha_mask)))*div_jpar/(exp(ddt(n)).DC());
      n -= a_dw*smooth_xz(smooth_xz(smooth_xz(alpha_mask)))*div_jpar/(exp(ddt(n)).DC());
    }
    //n += mu * (LapXZ(ddt(n))+ Grad(ddt(n))*Grad(ddt(n))) ;
    //n -= alpha; //very very slow term
   
  }
  else {
    u += beta* DDZ(ddt(n)+n0)/(n0); 
    n -= bracket3D(ddt(phi),ddt(n));

    //n += mu * LapXZ(ddt(n)) ;
    n -= alpha*ddt(n);

  }
 
  // n -= bracket3D(ddt(phi),ddt(n));
  // //n -= mybracket(ddt(phi),ddt(n));
  // n += mu * LapXZ(ddt(n));
  // n -= alpha* ddt(n);
  

  if(inc_jpar){

  // //   mesh->communicate(ddt(Te));
  // //   Te = 0;
    if (log_n)
      div_jpar = -a_dw*(ddt(n)*Te0 - ddt(phi));//*(log(n)*
    else {
      ddt(n).applyBoundary("dirichlet");
      div_jpar = -a_dw*(log(ddt(n) + n0)*Te0 - ddt(phi));
      //div_jpar = -pow(kpar,2.0)*(log(ddt(n) + n0)*Te0 - ddt(phi))/(fmei*.51*.1);
    }
    div_jpar = div_jpar - div_jpar.DC();
    div_jpar.applyBoundary();
   

    u += smooth_xz(alpha_mask)*div_jpar;
    n += smooth_xz(alpha_mask)*div_jpar;
  }

  u.applyBoundary("dirichlet");
  n.applyBoundary("dirichlet");
  phi.applyBoundary("dirichlet");
  //n.applyBoundary();
  //u.applyBoundary();
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
const Field3D smooth_xz(const Field3D &f){
  Field3D result;
  result.allocate();
  result = f;
  for(int x=2;x<mesh->ngx-2;x++)
    for(int y=0;y<mesh->ngy;y++)
      for(int z=2;z<mesh->ngz-2;z++) {
        result[x][y][z] = 0.5*f[x][y][z] + 0.125*( 0.5*f[x+1][y][z] + 0.125*(f[x+2][y][z] + f[x][y][z] + f[x+1][y][z-1] + f[x+1][y][z+1]) +
                                                   0.5*f[x-1][y][z] + 0.125*(f[x][y][z] + f[x-2][y][z] + f[x-1][y][z-1] + f[x-1][y][z+1]) +
                                                   0.5*f[x][y][z-1] + 0.125*(f[x+1][y][z-1] + f[x-1][y][z-1] + f[x][y][z-2] + f[x][y][z]) +
                                                   0.5*f[x][y][z+1] + 0.125*(f[x+1][y][z+1] + f[x-1][y][z+1] + f[x][y][z] + f[x][y][z+2]));
      }

  mesh->communicate(result);
  return result;
}

BoutReal Ullmann(double x, double Lx, double y,double Ly,double x_sol,double eps,
		 double m){
  
  
  int count = 0;
  bool hit_divert = false;
  bool inSOL;
  double q,qmax;

  //int max_orbit = 100;
  
  double L = 0.0;
  //double eps = .5;
  double aa = -.01;
  //double m = 3.0;
  double l = 10.0;
  double R = 85.0;
  double q0 = 3.0;
  double b = 68.0;
  double a= 60.0;
  
  double nu = 2.0;
 
  //q = q0;

  //double width = eps*3./5.; //very rough, .12 for eps = .2
  //double width = eps*(3./5.)*(3./m);
  double width = eps*30.;  
  double offset = x_sol * .2;

  //will cover from b(1-offset) to b(1- offset + .4)
  //x = a*(x_sol*(x/Lx)/2. + 1.-x_sol);
  //x = b*(a/b + 3.*(b - a)/b * (x/Lx));
  //x = a + 3.*(b - a) * (x/Lx);
  //x = b - 3.*b*width + 8*b*width*(x/Lx); //the chaotic region should be between 1/4 of the total domain size witht this setup
  //x = b*(40./55. + 2.*(55. - 40)/55. * (x/Lx));
  x =b - 3*width +  (3*(b-a)+3*width)*(x/Lx); // the total region size is then about 26 to about 35 cm 
  //x = b - width + 3*width*(x/Lx);
  y = y*(2.0*M_PI/Ly);
  // double xx = x_new/a


  //q = q0*pow(x/a,2.0)/(1.0-pow(1.0-x/a,nu+1.0)*double());  


  double x_new;
  double y_new;
  double x_new2;

  double C = ((2*m*l*pow(a,2.0))/(R*q0*pow(b,2.0)))*eps;
  //output<<x<<" "<<y<<endl;
  x_new = x;
  y_new = y;
  qmax = q0*pow(b/a,2.0); 
  while(count < max_orbit and not hit_divert){
    // x_new = x;
    // y_new = y;
    x_new = x_new/(1-aa*sin(y_new));
    //q = q0*pow((x_new/a),2.0);

    q = q0* pow(x_new/a,2.0)/(1.0-double(x_new<a)*pow(1.0-x_new/a,nu+1.0)); 
    if (q >qmax) 
      q = qmax;
    // q = q+ double(x_new>b)*q0*pow(b/a,2.0)/(1.0-pow(1.0-b/a,nu+1.0)); 
    
    C = ((2*m*l*pow(a,2.0))/(R*q*pow(b,2.0)))*eps;
    y_new =  (y_new+ 2*M_PI/q + aa*cos(y_new));
    y_new = fmod(y_new,2*M_PI);
    
    x_new2 = Newton_root(x_new,y_new,b,C,m);
    
    //output<< (-1.0*x_new+ x_new2 +(m*b*C)/(m-1.0) * pow(x_new2/b, m-1.0) * sin(m*y_new))<<endl;

    //output<< "old: " <<(-1.0*x_new+ x_new +(m*b*C)/(m-1.0) * pow(x_new/b, m-1.0) * sin(m*y_new))<<endl;

    //chi = (-x_new + x_out +(m*b*C)/(m-1)*(x_out/b)**(m-1) *np.sin(m*y_new))**2
    //x_new2 = (newton_krylov(func,x_new));
    
    //q = q0*pow(x_new2/a,2.0)/(1.0-w(1.0-x_new2/a,nu+1.0));  
    //C = ((2*m*l*pow(a,2.0))/(R*q*pow(b,2.0)))*eps;
    
    y_new = (y_new - C*pow(x_new2/b , m-2) * cos(m*y_new));
    y_new = fmod(y_new,2*M_PI);
    x_new = x_new2;
    //output <<x_new<<endl;
    hit_divert = (x_new > b or x>1.2*b or x_new <0);// or (x_new <  and x < b);
    count++;

    if (!hit_divert) {
      L = L + 2.0*M_PI*q*R;
    }
    else{
      L = L + 2.0*M_PI*qmax*R;
    }

    if(count == max_orbit){
      L = -1;
    }
    
  }
  //output<<L<<endl;
  //if (L == max_orbit) {L = 10.0;}
  return L;
}

BoutReal alphamap(double x, double Lx,double y,double Ly,
		  double k,double q0 ,double R,int max_orbit,double period,
		  bool count_turns){ 
  

  
  x = x*(2.0*M_PI/Lx)-2.0*M_PI;
  y = y*(2.0*M_PI/Ly);
  //output << "[" << x  << "], "<<endl;
  int count = 0;
  bool hit_divert = false;
  bool inSOL;
  double q;

  double L = 0.0;
  
  //count_turns = true;
  while(count < max_orbit and not hit_divert){
    inSOL = x>0.0; //are we in SOL now?
    
    //update the field line
    x = x+k*sin(period * y);
    y = fmod((x+y),(2.0*M_PI)); //one can argue that x = 2*M_PI*fmod(q(R),1)

    q =q0; //+x/(2*M_PI);
   
    if (count_turns)
      L = L +1.0;
    else
      L = L + q *(R/3.0)*2.0*M_PI; //just assume that the minor radius is 1/3 the major radius, ~ DIII-D
   
    hit_divert = x>0.0;
   
    count++;
    }
  // output << L << endl;
  return L;

}


BoutReal Newton_root(double x_in,double y_in,double b, double C, double m){
  double x_out = x_in;
  double atol = 1.0;
  int iter = 0;
  int max_iter = 300;
  double f;
  double J;

  while (atol > .0001 && iter < max_iter)
    {
      f = (-1.0*x_in + x_out +(m*b*C)/(m-1.0) * pow(x_out/b, m-1.0) * sin(m*y_in));
      J = 1.0 + (m*C) * pow(x_out/b,m-2.0) * sin(m*y_in);
      atol = fabs(f);
      x_out = x_out - 1.0*f/J;
      iter++;
      
    }
  //output<< (x_out - x_in)/fabs(x_in)<<endl;
  return x_out;
}
