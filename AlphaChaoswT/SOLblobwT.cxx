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
Field3D u, n,n_prev, u_prev,Te, Te_prev,Te_,sqrt_Te; //vorticity, density


//for the jacobian
Field3D uO,phiO,nO,TeO;
//derived variables
Field3D phi,brkt,bias_phi;
int phi_flags;

//other fields
Field3D test1, test2, ReyN, ReyU,B0,R,beta;
Field2D n0;
Field3D Lambda, Lambda_eff;

//Constrained 
Field3D C_phi;

//other params
BoutReal nu, mu,gam,alpha_c, eps,fmei,kpar,AA,ZZ,a_dw;
BoutReal TIMESTEP,DT,t0;
BoutReal bias_phi_0;

BoutReal x0;

Field3D alpha, temp,edgefld,alpha_s, alpha_j,source,sink_out,sink_in;
Field3D nave,uave,uDC,nDC,n_rms,u_rms;

Field3D alpha_mask,div_jpar;
BoutReal Te0;
//solver options
bool use_jacobian, use_precon;
bool evolve_te;
//experimental

bool withsource,wave_bc,diff_bc,withsink;
bool use_constraint;
string chaosalpha;
bool inc_jpar, alpha_inp;
bool log_n,log_Te;
BoutReal n_sol,rho_s;
int max_orbit,counter;

int m;

int MZ,NX;

FieldGroup comms; // Group of variables for communications

const Field3D mybracket(const Field3D &phi, const Field3D &A);
int jacobian(BoutReal t); // Jacobian-vector multiply
int precon(BoutReal t, BoutReal cj, BoutReal delta); // Preconditioner


//BoutReal alphamap(BoutReal x,BoutReal z);
BoutReal alphamap(double x, double Lx, double y,double Ly,
		  double k=1.20,double q0=3.0,double R=100,
		  int max_orbit =4000,double period=1.0,
		  bool count_turn = 0);

BoutReal Ullmann(double x, double Lx, double y,double Ly, BoutReal eps, double m);
BoutReal Newton_root(double x_in,double y_in,double b = 50.0, double C=.01, double m = 3.0);

const Field3D smooth_xz(const Field3D &f); 

BoutReal ramp(double t,double t_window);
//const Field3D GT(const Field3D &f,Field3D &g);

//int alphamapPy();
//int precon_phi(BoutReal t, BoutReal cj, BoutReal delta);
//int jacobian_constrain(BoutReal t); // Jacobian-vector multiply

int physics_init(bool restarting)
{
  
  

  Options *globaloptions = Options::getRoot();
  Options *options = globaloptions->getSection("sol");
  Options *solveropts = globaloptions->getSection("solver");

  OPTION(options, phi_flags, 0);
  OPTION(options, alpha_c, 3e-5);
  OPTION(options, nu, 2e-3);
  OPTION(options, mu, 2e-3);
  OPTION(options, gam, 1e1);
  //OPTION(options, beta, 6e-4);
  OPTION(options, inc_jpar,false);
  OPTION(options, alpha_inp,false);
  OPTION(options,chaosalpha,"chaos");
  
  OPTION(options,x0,100.0);
  
  OPTION(options, m, 3);
  OPTION(options,eps,.1);
  OPTION(options, max_orbit, 100);
  OPTION(options, bias_phi_0, 0.0);
  OPTION(options,rho_s,1.0);
  
  OPTION(globaloptions,MZ,33);
  //OPTION(globaloptions,nx,100);
  OPTION(globaloptions,TIMESTEP,1);
  DT = TIMESTEP;
  (globaloptions->getSection("mesh"))->get("nx", NX,  100);
  (globaloptions->getSection("te"))->get("evolve", evolve_te,   false);

  OPTION(solveropts,use_precon,false);
  OPTION(solveropts,use_jacobian,false);
  OPTION(solveropts,use_constraint,false);

  OPTION(options,withsource,false);
  OPTION(options,withsink,true);
  OPTION(options,wave_bc,true);
  OPTION(options,diff_bc,false);
  OPTION(options,log_n,true);
  OPTION(options,log_Te,false);

  output.write("sol settings: %g: \n",withsource);

  OPTION(options,Te0,1e0);
  //OPTION(options,n0,1e0);
  OPTION(options, AA, 2.0); //deutrium ?
  OPTION(options, ZZ, 1.0); //singly ionized
  
  
  //sweet way to create stuff from inp
  FieldFactory f(mesh);
  string s;
  (globaloptions->getSection("B0"))->get("function",s,"");
  B0 = f.create3D(s);
  //output.write("B0: %g\n",B0[0][0][0]);
  dump.add(B0,"B0",0);


  (globaloptions->getSection("beta"))->get("function",s,"");
  beta = rho_s * f.create3D(s);

  if (alpha_inp) {
    (globaloptions->getSection("alpha"))->get("function",s,"");
    alpha = rho_s * f.create3D(s);
    //alpha_mask =
  } else {
    alpha.allocate();
    alpha_j.allocate();
    alpha_mask.allocate();
    

    BoutReal Lx = mesh->dx[0][0]*NX;
   
    BoutReal ***a = alpha.getData();
    BoutReal ***a_j = alpha_j.getData();
    BoutReal ***a_m = alpha_mask.getData();
    double x_sol = 60.0;
    BoutReal Lxz = 0.0;
    for(int jz=0;jz<mesh->ngz;jz++) 
      for(int jx=0;jx<mesh->ngx;jx++){
	//Lxz is in cm  - not normalized . . .!
	Lxz = Ullmann(x0+ Lx*mesh->GlobalX(jx),1.0,mesh->dz*jz,mesh->zlength,eps,m);
	//output.write("x %g \n",x0+ Lx*mesh->GlobalX(jx));
	for(int jy=0;jy<mesh->ngy;jy++){
	  if  ("jump" == chaosalpha) {
	    a[jx][jy][jz]=alpha_c*double((x0 + Lx*mesh->GlobalX(jx)) > x_sol);
	    a_m[jx][jy][jz]= double(a[jx][jy][jz] == 0.); //double(mesh->GlobalX(jx) <= x_sol);
	  }
	  else {
	    a[jx][jy][jz]=(Lxz>0)*(2.0*rho_s/Lxz);
	    a_m[jx][jy][jz]=(Lxz<0);
	  }
	
	  if (mesh->firstX()){
	    a_m[0][jy][jz]=0;
	    a[0][jy][jz]= 0;
	  }
	}
      }
  } // 

  if ("smooth" == chaosalpha){
    alpha = lowPass(alpha,0);
    alpha_mask = lowPass(alpha_mask,0);
  }

  bout_solve(u, "u");
  comms.add(u);
  
  
  bout_solve(n, "n");
  comms.add(n);
 
  
  
  if(withsource){
    source = f.create3D("gauss(x-0.05,0.1)");
    dump.add(source,"source",0);
    
  }

  

  if(withsink){
    //sink = 1.0 - f.create3D("h(x-.9)");
    sink_out = f.create3D("gauss(x-1.0,.05)");
    //sink_in = f.create3D("gauss(x,.04)");

    sink_out = sink_out/(max(sink_out,1));
    // sink_in = sink_in/(max(sink_in,1));			   

    //dump.add(sink_in,"sink_in",0);
    dump.add(sink_out,"sink_out",0);
    
  }

  if(inc_jpar){
    dump.add(div_jpar,"div_jpar",1);
    div_jpar.setBoundary("phi");
    
    // comms.add(div_jpar);
  }
  
  if(evolve_te) {
    bout_solve(Te, "Te");
    comms.add(Te);
    //rhscomms.add(ddt(Te));
    output.write("te\n");
    Te_prev = Te;
    output.write("in setup TE: %g\n",min(Te,true));



  }else
      Te = Te0;
  
  
  //phi = invert_laplace(u, phi_flags);
  static Field2D A = 0.0;
  static Field2D C = 1e-12;
  static Field2D D = 1.0;
  
  Lambda = 4.7;
  phi = Te*Lambda;
  bias_phi = bias_phi_0 * Lambda*(f.create3D("h(x-.25)") + f.create3D("h(.5-x)")-1.0) ;
  int smooth_n=5;
  while (smooth_n>0) {
    bias_phi = smooth_x(bias_phi);
    --smooth_n;
  }


  //bias_phi = smooth_x(bias_phi_0 * Lambda*(f.create3D("h(x-.25)") + f.create3D("h(.5-x)")-1.0)) ;
  dump.add(bias_phi,"bias_phi",0);
 

  
  //alpha = alpha_c;
  fmei  = 1./1836.2/AA;
  kpar = alpha_c;  //simplest possible 
  a_dw = pow(kpar,2.0)/(fmei*.51*.1);
  a_dw = .25;


  // BoutReal edge[mesh->ngz];
  // BoutReal zoomfactor = 3.0;
  // BoutReal lowR = .55;
  // BoutReal Lxz = 0;
  // BoutReal x_sol = .3;
  // BoutReal rho_s = .5;
  
  dump.add(alpha,"alpha",0);
  dump.add(alpha_mask,"alpha_mask",0);
  dump.add(alpha_s,"alpha_smooth",0);
  // dump.add(eps,"eps",0);
  dump.add(m,"m",0);
  dump.add(a_dw,"a_dw",0);
  


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
  output.write("with source %i \n",withsource);
  output.write("use phi_flags %i \n",phi_flags);

  output.write("DONE WITH PHYSICS_INIT\n");

  n0 = 1e-4;https://github.com/boutproject/BOUT
  n_prev =n;
  u_prev = u;
  nave = n;
  uave = u;


  if (log_n){
    n = log(n+n0);
    n0 = log(n0);
    n_prev =n;
  }

 if (log_Te){
    Te = log(Te+1e-4);
    Te0 = log(Te0);
    dump.add(Te_,"Te_",1);
  }



  ddt(n).setBoundary ("ddt[n]") ;
  // phi.setBoundary("phi");
  ddt(u).setBoundary ("ddt[u]") ;
 
  mesh->communicate(comms);

  if (restarting){
    t0 = solver->getCurrentTimestep();
    output.write("t0: %g  \n" ,t0);
  }

  return 0;
}

#define bracket3D(f, g) ( b0xGrad_dot_Grad(f, g) )
#define LapXZ(f)(mesh->g11*D2DX2(f) + mesh->g33*D2DZ2(f))
#define LapZ(f)(mesh->g33*D2DZ2(f))
//#define LapXZ(f)(D2DX2(f));// + mesh->g33*D2DZ2(f))

int physics_run(BoutReal t)
{
  
  if (counter == 0){
    t0 = t;
    counter++;
    
  
  }
  //output.write("t: %g\n",t);
  // Run communications
  mesh->communicate(comms);
  //phi = invert_laplace(u, phi_flags);
  
  //output.write("in run start Te: %g\n",min(Te,true));
  if (log_Te){
    Te_ = exp(Te);
    sqrt_Te = sqrt(Te_);
  }
  else{
    sqrt_Te = sqrt(Te);
    Te_ = Te;
  }
  //output.write("in run start 2 Te: %g\n",min(Te,true));
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
	    // n[i][j][k] =(.05*n[i+1][j][k] + n_prev[i][j][k])/(1.0+.05);
	    n[i][j][k] = n_prev[i][j][k];
	    //n[i][j][k] =(.5*(nDC[i+1][j][k]) + (nprevDC)[i][j][k])/(1.0+.5);
	    //u[i][j][k] =(.5*u[i+1][j][k] + u[i][j][k])/(1.0+.5);
	    //n[i][j][k] = (n.DC)[i][j][k];
	    //u[i][j][k] = (u.DC)[i][j][k];
	  }
    

    
    
  } else{
    n.applyBoundary();
    //output.write("Te: %g\n",Te[2][1][0]);

    
    // if(catch_neg)
    //   Te = Te
  }
  if (evolve_te)
    Te.applyBoundary(); //need this
  static Field2D A = 0.0;
  static Field2D C = 1e-10;
  static Field2D D = 1.0;
  FieldFactory f(mesh);

  u.applyBoundary(); //BIG speed up
  // ostringstream convert; 
  // convert << Lambda[0][0];
  u_prev =u;

  //awkward (but ok) way to set bc on phi
  if (mesh->firstX())
    for(int i=1;i>=0;i--)
      for(int j =0;j< mesh->ngy;j++)
	for(int k=0;k < mesh->ngz; k++)
	  u_prev[i][j][k] = (Lambda*Te)[i][j][k];

  if (mesh->lastX())
    for(int i=1;i>=0;i--)
      for(int j =0;j< mesh->ngy;j++)
	for(int k=0;k < mesh->ngz; k++)
	  u_prev[mesh->ngx-1-i][j][k] = (Lambda*Te)[mesh->ngx-1-i][j][k];
	    


  phi = invert_laplace(u_prev, 49152,&A,&C,&D);
  //fast way to incluce bias_phi in the physics without recoding
  Lambda_eff  = Lambda;// - ramp(t-t0,20)*bias_phi/Te_;
  //output.write("ramp  %g \n",ramp(t-t0,20*DT));
  
  mesh->communicate(comms);
  //mesh->communicate(phi);
  ddt(u)=0;
  ddt(n)=0;
  
  if (evolve_te)
    ddt(Te) = 0;
 

  //output.write("n in run: %g\n",n[5][1][0]);
 
  // B0 = 1.0;


 
  ddt(u) -= (1.0/B0)*bracket3D(phi,u);
  //ddt(u) += alpha * sqrt_Te*(1 - exp(Lambda -phi/Te_ - ramp(t-t0,20*DT)*bias_phi/Te_)); //
  ddt(u) += alpha * sqrt_Te*(1.0 - exp(Lambda_eff - phi/Te_));
  ddt(u) += nu * LapXZ(u);



  if (log_n){
    //ddt(u) += 2*B0^2*beta*(Te*DDZ(n) + DDZ(Te)); //slooow, B^2 != B*B, go figure, 

    ddt(u) += 2.0*B0*B0*beta*(Te_*DDZ(n) + DDZ(Te_));
    ddt(n) -= (1.0/B0)*bracket3D(phi,n);
       
    ddt(n) += mu * (LapXZ(n) + Grad(n)*Grad(n)) ; //boundary issues?
  
    ddt(n) -= alpha* sqrt_Te*exp(Lambda_eff -phi/Te_);
    ddt(n) += 2.0*(beta/B0)*(DDZ(Te_-phi)+Te_*DDZ(n)); //slow
   
  } else {
    // ddt(u) += beta* DDZ(n+n0)/(n+n0);
   
    // //output.write ("n log_n \n");
    // // ReyN = bracket3D(phi,n)/(mu * LapXZ(n)+1e-5);
    
    // ddt(n) -= bracket3D(phi,n);
    // ddt(n) += mu * (LapXZ(n)) ;
    // ddt(n) -= alpha* n;
  }
  
  //ddt(n) = n*0;
  //ddt(u) = u*0;
  
  
 
  if(withsource ){
    if (log_n)
      ddt(n) += (1.0e0 * max(alpha,1) * source)/exp(n);

      //ddt(n) += (5.0e-1 * alpha_c * source)/exp(n);
    else
      ddt(n) += (1.0e0 * alpha_c * source);
    // if (evolve_te)
    //   ddt(Te) += (1.0e0 * max(alpha,1)*source);
  }

  if(withsink){
    BoutReal target_sol,target_core;
   

    Field3D region_select = f.create3D("h(.98-x) +h(x-.7)")-1.0;
    target_sol = min((n* region_select).DC(),true)-.02;

    region_select = f.create3D("h(.3-x) +h(x-.02)")-1.0;
    target_core = min((n* region_select).DC(),true)-.02;

    //output.write("tarval_val  %g \n",target_val);
    //target_val = -6.0;
    //ddt(n) -= (1e0*alpha*sink_in)*(1.0-exp(target_core)/exp(n));
    ddt(n) -= (1e0*alpha*sink_out)*(1.0-exp(target_sol)/exp(n));
    //ddt(u) = lowPass(ddt(u)*sink_out,4) +  ddt(u)*(1.0-sink_out);
    //ddt(u) = lowPass(ddt(u)*(sink_sol+sink_core),int(MZ/8.0)) +  ddt(u)*(1.0-sink_sol-sink_core);
    // if (evolve_te){
    //   // region_select = f.create3D("h(.95-x) +h(x-.7)")-1.0;
    //   // target_sol = min((Te* region_select).DC(),true)*.98;

    //   // region_select = f.create3D("h(.3-x) +h(x-.05)")-1.0;
    //   // target_core = min((Te* region_select).DC(),true)*.98;
      
    //   // ddt(Te) -= (1e1*alpha*sink_core)*(Te - target_core);
    //   // ddt(Te) -= (1e1*alpha*sink_sol)*(Te - target_sol);
    //   //try low-pass filter on boundary
    
    //   ddt(u) = lowPass(ddt(u)*(sink_sol+sink_core),5) +  ddt(u)*(1.0-sink_sol-sink_core);
      
      

    // }
    
   
    
    //ddt(u) -= (3e0*alpha_c *sink)*(u - u.DC());


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
      div_jpar = -a_dw*(n*Te0 - phi);//*(log(n)*
    else {
      //phi = phi - lazy_log(n[0][0][0])*Te0;
      div_jpar = -a_dw*(log(abs(n+n0))*Te0 - phi);//*(log(n)

    }
    div_jpar = div_jpar - div_jpar.DC();
    div_jpar.applyBoundary();
    // div_jpar.applyBoundary();
    // //for values where alpha  = min

    ddt(u) += smooth_xz(smooth_xz(smooth_xz(alpha_mask)))*div_jpar;
    ddt(n) += smooth_xz(smooth_xz(smooth_xz(alpha_mask)))*div_jpar;
 
  }

  ddt(Te) = 0.0;
  if(evolve_te) {
    //Te.applyBoundary();
    ddt(Te)  -= (1.0/B0)*bracket3D(phi,Te);
    ddt(Te) += mu * (LapXZ(Te));
    if (log_n)
      ddt(Te) += (2./3.)* Te * (beta/B0)*(DDZ(Te - phi)+Te *DDZ(n)); //basically ddt(n) 

    //output.write("in run Te: %g\n",min(Te,true));

    ddt(Te) += 5./3. * Te*(beta/B0) * DDZ(Te);
    //the following term seems to cause negative Te .  not sure why
    ddt(Te) -= (2./3. * alpha*Te*sqrt_Te*(1.71*exp(Lambda_eff -phi/Te) - .71));

    //output.write("max alpha  %g\n",max(alpha,1));


    if (withsource)
      ddt(Te) += (1.0e0 * max(alpha,1)*source);
    if (withsink){
      //ddt(Te) = lowPass(ddt(Te)*sink_out,5.0) +  ddt(Te)*(1.0-sink_out);
      //ddt(Te) += 1e2*sink_out*mu * (LapXZ(Te));
      // Field3D region_select = f.create3D("h(.98-x) +h(x-.7)")-1.0;
      // test1 = 100.*(1.0-region_select)+(Te* region_select);
      // BoutReal target_sol;
      // target_sol = min(test1.DC(),true);
      
      // //target_Te = smooth_x(smooth_x(Te.DC()));
      // ddt(Te) -= (3e0*alpha*sink_out)*(Te-target_sol);
      //ddt(Te) += 5e0*sink_out*mu * (LapXZ(Te));
      ddt(Te) += 5e0*sink_out*mu * (LapZ(Te));
      
      //ddt(Te) = sink_out*smooth_x(ddt(Te).DC()) + (1.0 - sink_out)*ddt(Te);
      
      //output.write("in run ddt(Te:) %g\n",DT * min(ddt(Te),true));

    }
    // if (log_Te){
    //   if (log_n)
    // 	ddt(Te) += (2./3.)* (beta/B0)*(DDZ(Te_-phi)+Te_*DDZ(n)); //basically ddt(n) 

    //   ddt(Te) += 5./3. * (beta/B0) * DDZ(Te_);
    //   ddt(Te) -= 2./3. * alpha*sqrt(Te_)*(1.71*exp(Lambda -phi/Te_) - .71);
    // }
    // else {
    //   if (log_n)
    // 	ddt(Te) += (2./3.)* Te * (beta/B0)*(DDZ(Te - phi)+Te *DDZ(n)); //basically ddt(n) 

    //   ddt(Te) += 5./3. * Te*(beta/B0) * DDZ(Te);
    //   ddt(Te) -= 2./3. * alpha*Te*sqrt(Te)*(1.71*exp(Lambda -phi/Te) - .71);
    // }

    

    // if (log_Te)
    //   ddt(Te) += (mu/10.) * (LapXZ(Te) + Grad(Te)*Grad(Te));
    // else
    //   ddt(Te) += (mu/10.) * (LapXZ(Te));

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
  mesh->communicate(ddt(u),ddt(n),ddt(Te));
  

  
  static Field2D A = 0.0;
  static Field2D C = 1e-12;
  static Field2D D = 1.0;
  
  ddt(phi) = invert_laplace(ddt(u), phi_flags,&A,&C,&D);
  //ddt(phi) = invert_laplace(ddt(u), phi_flags); 
  // ddt(u).applyBoundary("dirichlet");
  // ddt(n).applyBoundary("dirichlet");
  // ddt(phi).applyBoundary("dirichlet");
  mesh->communicate(ddt(phi),ddt(u),ddt(n),ddt(Te));
  mesh->communicate(comms); //need this !
  // u=0;
  // n=0;
  // Te=0
    
  uO =u;
  nO = n;
  TeO = Te;
  phiO = phi;
 
  u -= (1.0/B0)*bracket3D(phiO,ddt(u));
  u -= (1.0/B0)*bracket3D(ddt(phi),uO+ddt(u));

  //u += alpha * sqrt(TeO)*(-Lambda_eff + phiO/TeO + ddt(phi)/TeO - ddt(Te)*phiO/TeO^2); //seems to speed up
  u += alpha * sqrt(TeO)*(ddt(phi)/TeO - ddt(Te)*phiO/TeO^2);
  u += alpha * .5/sqrt(TeO)*ddt(Te)*(1 - exp(Lambda_eff - phiO/TeO));
  
  

  if(evolve_te) {
    
    Te  -= (1.0/B0)*bracket3D(ddt(phi),TeO + ddt(Te));
    Te  -= (1.0/B0)*bracket3D(phiO,ddt(Te) );


    Te += 5./3. * TeO *(beta/B0) * DDZ(ddt(Te));
    Te += 5./3. * ddt(Te) *(beta/B0) * DDZ(TeO);

    Te -= 2./3. * alpha*3./2.*sqrt(TeO)*ddt(Te)*((1.71*exp(Lambda_eff -phiO/TeO)) - .71);

    Te -= 2./3. * alpha*TeO*sqrt(TeO)*(1.71*(1.0 + (0 - ddt(phi)/TeO + phiO*ddt(Te)/TeO^2)) - .71);

    Te += 5./3. * ddt(Te)*(beta/B0) * DDZ(TeO + ddt(Te)); 
    //Te += mu * (LapXZ(Te));

  }
  if (log_n){
    
    n -= bracket3D(ddt(phi),nO+ddt(n));
    n -= bracket3D(phiO,ddt(n));
    

    n += 2.0*(beta/B0)*(DDZ(ddt(Te)-ddt(phi)));
    
    
    
    u += 2.0*B0*B0*beta*((TeO+ddt(Te))*DDZ(ddt(n)) + DDZ(ddt(Te)));
    u += 2.0*B0*B0*beta*ddt(Te)*DDZ(nO);

    //n += mu * (LapXZ(ddt(n))+ 2.0*Grad(nO)*Grad(ddt(n))) ; //slows things down
   

    n -= alpha* .5/sqrt(TeO)*ddt(Te)*exp(Lambda_eff - phiO/TeO);  
    n -= alpha* sqrt(TeO)*(1 + 0*(Lambda_eff -phiO/TeO) - ddt(phi)/TeO + phiO*ddt(Te)/TeO); 

   
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
  Te.applyBoundary("dirichlet");
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

BoutReal ramp(double t,double t_window){
  if (t>t_window)
    return 1.0;
  else
    return t/t_window;
}

BoutReal Ullmann(double x, double Lx, double y,double Ly,double eps,
		 double m){
  
  
  int count = 0;
  bool hit_divert = false;
  bool inSOL;
  double q,qmax,C;

  //int max_orbit = 100;
  
  double L = 0.0;
  //double eps = .5;
  double aa = -.03;
  //double m = 3.0;
  double l = 10.0;
  double R = 85; // make this vary ? no keep it fixed 
  double q0 = 3.0;
  double b = 68.0;
  double a= 60.0;
  
  double nu = 2.0;
 
 
  y = y*(2.0*M_PI/Ly);
 

  double x_new;
  double y_new;
  double x_new2;


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
   
    y_new = (y_new - C*pow(x_new2/b , m-2) * cos(m*y_new));
    y_new = fmod(y_new,2*M_PI);
    x_new = x_new2;
    //output <<x_new<<endl;
    hit_divert = (x_new > b or x>b or x_new <0);// or (x_new <  and x < b);
    count++;

    if (!hit_divert) {
      //L = L + 2.0*M_PI*q*R;
      L = L + q*R;
    }
    else{
      //L = L + 2.0*M_PI*qmax*R;
      L = L + qmax*R;
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
