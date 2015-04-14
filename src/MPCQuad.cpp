/*
 *    This file is part of ACADO Toolkit.
 *
 *    ACADO Toolkit -- A Toolkit for Automatic Control and Dynamic Optimization.
 *    Copyright (C) 2008-2014 by Boris Houska, Hans Joachim Ferreau,
 *    Milan Vukov, Rien Quirynen, KU Leuven.
 *    Developed within the Optimization in Engineering Center (OPTEC)
 *    under supervision of Moritz Diehl. All rights reserved.
 *
 *    ACADO Toolkit is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    ACADO Toolkit is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with ACADO Toolkit; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


 /**
  *    \file   examples/ocp/rocket.cpp
  *    \author Boris Houska, Hans Joachim Ferreau
  *    \date   2009
  */

#include <acado_optimal_control.hpp>
#include <acado_gnuplot.hpp>
#include <spline.h>
#include <robot.h>
#include <quadconf.h>
#include <locale>
#include <trajectory.h>
#include <acado_code_generation.hpp>
#include <acado/process/process.hpp>
#include <acado_toolkit.hpp>

#include <time.h>

/* >>> start tutorial code >>> */
clock_t t;


#define ALGO 0 //MPC, OCP or GENE
#define INIT 1



int main( ){

    USING_NAMESPACE_ACADO

    // INTRODUCE THE VARIABLES:
    // -------------------------
    DifferentialState     x,y,z,vx,vy,vz,phi,theta,psi,p,q,r,u1,u2,u3,u4;
    Control               vu1,vu2,vu3,vu4;
    Parameter		  T;
    DifferentialEquation  f;


    // Quad constants
    const double c = 0.00001;
    const double Cf = 0.00065;
    const double d = 0.250;
    const double Jx = 0.018;
    const double Jy = 0.018;
    const double Jz = 0.026;
    const double m = 0.9;
    const double g = 9.81;
    const double Cx = 0.1;

    // Minimization Weights
    double coeffU = 0.00000000000000000000001;
    double coeffX = 0.00001;
    double coeffX2 = 10.;

    // final values
    double xf = 1, yf = 0., zf = 8.;

    // Temps max
    double tps = 6.;
    int nb_it = 20;
    double tmpc = 0.2;
    int mpc_nb_it = 2;



    // DEFINE A DIFFERENTIAL EQUATION:
    // -------------------------------
    OCP ocp( 0.0, tps, nb_it );
//    ocp.minimizeMayerTerm( 0 );
//    ocp.minimizeLagrangeTerm(coeffX*(sqrt(1+(x-xf)*(x-xf)) + sqrt(1+(y-yf)*(y-yf)) + sqrt(1+(z-zf)*(z-zf)) + r*r*50 + p*p*50 + q*q*50 +exp((16 - (x*x+(z-5)*(z-5)))*0.01)) + coeffU*( sqrt(1+(vu1)*(vu1))+sqrt(1+(vu2)*(vu2))+sqrt(1+(vu3)*(vu3))+sqrt(1+(vu4)*(vu4))));
//    ocp.minimizeLagrangeTerm(((x-xf)*(x-xf) + (y-yf)*(y-yf) + (z-zf)*(z-zf))*coeffX);

    //  ------------------------- LSQ MINIMISATION ------------------------ //
    Function h, hf;
    h << x;
    h << y;
    h << z;
    h << p;
    h << q;
    h << r;
    h << vu1;
    h << vu2;
    h << vu3;
    h << vu4;
//    h << exp((16 - (x*x+(z-5)*(z-5)))*0.01);

    hf << x;
    hf << y;
    hf << z;

    DMatrix Q(10,10), Qf(3,3);
    Q(0,0) = coeffX;
    Q(1,1) = coeffX;
    Q(2,2) = coeffX;
    Q(3,3) = coeffX*coeffX2;
    Q(4,4) = coeffX*coeffX2;
    Q(5,5) = coeffX*coeffX2;
    Q(6,6) = coeffU;
    Q(7,7) = coeffU;
    Q(8,8) = coeffU;
    Q(9,9) = coeffU;
//    Q(10,10) = coeffX*0.03;

    Qf(0,0) = coeffX*100;
    Qf(1,1) = coeffX*100;
    Qf(2,2) = coeffX*100;

    DVector reff(3), ref(10);
    ref(0) = xf;
    ref(1) = yf;
    ref(2) = zf;
    ref(3) = 0.;
    ref(4) = 0.;
    ref(5) = 0.;
    ref(6) = 0.;
    ref(7) = 0.;
    ref(8) = 0.;
    ref(9) = 0.;
//    ref(10) = 0.;

    reff(0) = xf;
    reff(1) = yf;
    reff(2) = zf;


    ocp.minimizeLSQ ( Q, h, ref);
//    ocp.minimizeLSQEndTerm(Qf, hf, reff);
    // ------------------------------------------------------------------ //


    f << dot(x) == vx;
    f << dot(y) == vy;
    f << dot(z) == vz;
    f << dot(vx) == -Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*sin(theta)/m;
    f << dot(vy) == Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*sin(phi)*cos(theta)/m;
    f << dot(vz) == Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*cos(phi)*cos(theta)/m - g;
    f << dot(phi) == p + sin(theta)*r;
    f << dot(theta) == cos(phi)*q-sin(phi)*cos(theta)*r;
    f << dot(psi) == sin(phi)*q+cos(phi)*cos(theta)*r;
    f << dot(p) == (d*Cf*(u1*u1-u2*u2)+(Jy-Jz)*q*r)/Jx;
    f << dot(q) == (d*Cf*(u4*u4-u3*u3)+(Jz-Jx)*p*r)/Jy;
    f << dot(r) == (c*(u1*u1+u2*u2-u3*u3-u4*u4)+(Jx-Jy)*p*q)/Jz;
    f << dot(u1) == vu1;
    f << dot(u2) == vu2;
    f << dot(u3) == vu3;
    f << dot(u4) == vu4;

    // DEFINE AN OPTIMAL CONTROL PROBLEM:
    // ----------------------------------
    
    //Dynamic
    ocp.subjectTo( f );

// ------------------------------------------------------------------------------------ //



    // ----------------------------- Set initial guess -------------------------- //



#if INIT==0 || INIT==1
    Grid timeGrid(0.0,tps,nb_it+1);
    VariablesGrid x_init(16, timeGrid);
#endif

#if INIT==0     // init with libquadspline
    double vmax=100;
    double amax=20;
    double jmax=40;
    double smax=320;

    Robot robot=Robot(1.0);

    robot.addDof(vmax,amax,jmax,smax);
    robot.addDof(vmax,amax,jmax,smax);
    robot.addDof(vmax,amax,jmax,smax);
    robot.addDof(1,2,1,1,true);

    LocalPath path(QuadConf(&robot, 0,0,0,0,0,0,0,0,0), QuadConf(&robot, xf,yf,zf,0,0,0,0,0,0));

    VariablesGrid acc(3,timeGrid);
    for (int i = 0 ; i<nb_it+1 ; i++ ) {
        QuadConf conf = path.getQuadConfAt(double(i)*tps/double(nb_it));
        x_init(i,0) = conf.getX();
        x_init(i,1) = conf.getY();
        x_init(i,2) = conf.getZ();
        x_init(i,3) = conf.getVx();
        x_init(i,4) = conf.getVy();
        x_init(i,5) = conf.getVz();
        x_init(i,6) = conf.getRoll();
        x_init(i,7) = -conf.getPitch();
        x_init(i,8) = conf.getYaw();
        x_init(i,9) = conf.getDroll();
        x_init(i,10) = -conf.getDpitch();
        x_init(i,11) = conf.getDyaw();
        x_init(i,12) = 58;
        x_init(i,13) = 58;
        x_init(i,14) = 58;
        x_init(i,15) = 58;
        acc(i,0) = conf.getAx();
        acc(i,1) = conf.getAy();
        acc(i,2) = conf.getAz();
    }
#endif

#if INIT==1      // init with quasistatic solution
    for (int i = 0 ; i<nb_it+1 ; i++ ) {
        x_init(i,0) = xf*double(i)/double(nb_it);
        x_init(i,1) = yf*double(i)/double(nb_it);
        x_init(i,2) = zf*double(i)/double(nb_it);
        x_init(i,3) = 0.;
        x_init(i,4) = 0.;
        x_init(i,5) = 0.;
        x_init(i,6) = 0.;
        x_init(i,7) = 0.;
        x_init(i,8) = 0.;
        x_init(i,9) = 0.;
        x_init(i,10) = 0.;
        x_init(i,11) = 0.;
        x_init(i,12) = 58;
        x_init(i,13) = 58;
        x_init(i,14) = 58;
        x_init(i,15) = 58;
    }
#endif
//    x_init.print();
 // ----------------------------------------------------------------------------------- //

    // ---------------------------- DEFINE CONTRAINTES --------------------------------- //


    // Min/Max time
//    ocp.subjectTo(1 <= T <= 10);

#if ALGO==1 //|| ALGO==0
    //Start conditions
    ocp.subjectTo( AT_START, x == 0 );
    ocp.subjectTo( AT_START, y == 0  );
    ocp.subjectTo( AT_START, z ==  0 );
    ocp.subjectTo( AT_START, vx == 0  );
    ocp.subjectTo( AT_START, vy ==  0 );
    ocp.subjectTo( AT_START, vz ==  0 );
//    ocp.subjectTo( AT_START, phi ==  0.0 );
//    ocp.subjectTo( AT_START, theta ==  -0. );
//    ocp.subjectTo( AT_START, psi ==  -0.);
    ocp.subjectTo( AT_START, p ==  0.0 );
    ocp.subjectTo( AT_START, q ==  0.0 );
    ocp.subjectTo( AT_START, r ==  0.0 );
    ocp.subjectTo( AT_START, (d*Cf*(u1*u1-u2*u2)+(Jy-Jz)*q*r)/Jx == 0 );
    ocp.subjectTo( AT_START, (d*Cf*(u4*u4-u3*u3)+(Jz-Jx)*p*r)/Jy == 0 );
    ocp.subjectTo( AT_START, (c*(u1*u1+u2*u2-u3*u3-u4*u4)+(Jx-Jy)*p*q)/Jz ==  0.0 );
    ocp.subjectTo( AT_START, -Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*sin(theta)/m ==  0);
    ocp.subjectTo( AT_START, Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*sin(phi)*cos(theta)/m ==   0 );
    ocp.subjectTo( AT_START, Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*cos(phi)*cos(theta)/m - g ==  0 );
#endif

    //end conditions
//    ocp.subjectTo( AT_END  , (x-xf)*(x-xf) <= 100);//(x_init(nb_it,0)-xf)*(x_init(nb_it,0)-xf) );
//    ocp.subjectTo( AT_END  , (y-yf)*(y-yf) <= 10);//(x_init(nb_it,1)-yf)*(x_init(nb_it,1)-yf) );
//    ocp.subjectTo( AT_END  , (z-zf)*(z-zf) <= 10);//(x_init(nb_it,2)-zf)*(x_init(nb_it,2)-zf) );

//    ocp.subjectTo( AT_END  , x ==  xf );
//    ocp.subjectTo( AT_END  , y ==  yf );
//    ocp.subjectTo( AT_END  , z ==  zf );
//    ocp.subjectTo( AT_END  , vx ==  0 );
//    ocp.subjectTo( AT_END  , vy ==   0);
//    ocp.subjectTo( AT_END  , vz ==  0 );
////    ocp.subjectTo( AT_END  , phi ==  0.0 );
////    ocp.subjectTo( AT_END  , theta ==  -0. );
////    ocp.subjectTo( AT_END  , psi == -0.0 );
//    ocp.subjectTo( AT_END  , p ==  0.0 );
//    ocp.subjectTo( AT_END  , q ==  0.0 );
//    ocp.subjectTo( AT_END  , r ==  0.0 );
//    ocp.subjectTo( AT_END, (d*Cf*(u1*u1-u2*u2)+(Jy-Jz)*q*r)/Jx ==  0.0 );
//    ocp.subjectTo( AT_END, (d*Cf*(u4*u4-u3*u3)+(Jz-Jx)*p*r)/Jy ==  0.0 );
//    ocp.subjectTo( AT_END, (c*(u1*u1+u2*u2-u3*u3-u4*u4)+(Jx-Jy)*p*q)/Jz ==  0.0);
//    ocp.subjectTo( AT_END, -Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*sin(theta)/m ==  0 );
//    ocp.subjectTo( AT_END, Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*sin(phi)*cos(theta)/m == 0.  );
//    ocp.subjectTo( AT_END, Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*cos(phi)*cos(theta)/m - g ==  0 );

    //Command constraints
    //on each motor
    ocp.subjectTo( 16 <= u1 <= 95 );
    ocp.subjectTo( 16 <= u2 <= 95 );
    ocp.subjectTo( 16 <= u3 <= 95 );
    ocp.subjectTo( 16 <= u4 <= 95 );
    ocp.subjectTo( -100 <= vu1 <= 100 );
    ocp.subjectTo( -100 <= vu2 <= 100 );
    ocp.subjectTo( -100 <= vu3 <= 100 );
    ocp.subjectTo( -100 <= vu4 <= 100 );

    //Obstacle constraints
#if ALGO==0 || ALGO==1
    ocp.subjectTo( 40000 <= (x*x+2*(z-3)*(z-3))*10000 );
    ocp.subjectTo( 40000 <= ((x-2)*(x-2)+2*(z-6)*(z-6))*10000 );
//    ocp.subjectTo( z<=5);
#endif
     // ----------------------------------------------------------------------------------- //




#if ALGO == 0
    // SETTING UP THE (SIMULATED) PROCESS:
    // -----------------------------------
    OutputFcn identity;
    DynamicSystem dynamicSystem( f,identity );
    Process process( dynamicSystem,INT_RK45 );

    // SETTING UP THE MPC CONTROLLER:
    // ------------------------------
    RealTimeAlgorithm alg( ocp,tmpc );
    alg.set( MAX_NUM_ITERATIONS, mpc_nb_it );
    StaticReferenceTrajectory zeroReference;
    Controller controller( alg,zeroReference );

    alg.initializeDifferentialStates(x_init);

    // SET OPTION AND PLOTS WINDOW
    // ---------------------------
    alg.set(DISCRETIZATION_TYPE, SINGLE_SHOOTING);
//        alg.set( HESSIAN_APPROXIMATION, EXACT_HESSIAN );
    GnuplotWindow window1(PLOT_AT_EACH_ITERATION);
    window1.addSubplot( x,"DifferentialState x" );
    window1.addSubplot( vx,"DifferentialState vx" );
    window1.addSubplot( -Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*sin(theta)/m,"acc x" );
    window1.addSubplot( z,"DifferentialState z" );
    window1.addSubplot( (x*x+(z-5)*(z-5)),"R^2" );
    alg<<window1;


    // SETTING UP THE SIMULATION ENVIRONMENT,  RUN THE EXAMPLE...
    // ----------------------------------------------------------
    SimulationEnvironment sim( 0.0,10.,process,controller );

    DVector x0(16);
    x0.setZero();
    x0(12) = 58.;
    x0(13) = 58.;
    x0(14) = 58.;
    x0(15) = 58.;

    t = clock();
    if (sim.init( x0 ) != SUCCESSFUL_RETURN)
        exit( EXIT_FAILURE );
    if (sim.run( ) != SUCCESSFUL_RETURN)
        exit( EXIT_FAILURE );
    t = clock() - t;
    std::cout << "total time : " << (((float)t)/CLOCKS_PER_SEC)<<std::endl;

    // ...AND PLOT THE RESULTS
    // ----------------------------------------------------------
    VariablesGrid sampledProcessOutput;
    sim.getSampledProcessOutput( sampledProcessOutput );
    std::ofstream file;
    file.open("/tmp/log.txt",std::ios::out);
    sampledProcessOutput.print(file);

    VariablesGrid feedbackControl;
    sim.getFeedbackControl( feedbackControl );

    GnuplotWindow window;
    window.addSubplot( sampledProcessOutput(0), "x " );
    window.addSubplot( sampledProcessOutput(1), "y " );
    window.addSubplot( sampledProcessOutput(2), "z " );
    window.addSubplot( sampledProcessOutput(6),"phi" );
        window.addSubplot( sampledProcessOutput(7),"theta" );
            window.addSubplot( sampledProcessOutput(8),"psi" );

    //               window.addSubplot( sampledProcessOutput(3), "Wheel Velocity [m/s]" );
    //               window.addSubplot( feedbackControl(1),      "Damping Force [N]" );
    //               window.addSubplot( feedbackControl(0),      "Road Excitation [m]" );
    window.plot( );

#endif

#if ALGO == 1

    // DEFINE A PLOT WINDOW:
    // ---------------------
    GnuplotWindow window2(PLOT_AT_EACH_ITERATION);
        window2.addSubplot( z,"DifferentialState z" );
        window2.addSubplot( x,"DifferentialState x" );
//        window2.addSubplot( x_init(0),"DifferentialState x" );

        window2.addSubplot( y,"DifferentialState y" );
        window2.addSubplot( vx,"DifferentialState vx" );
//        window2.addSubplot( vz,"DifferentialState vz" );
//        window2.addSubplot( Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*cos(phi)*cos(theta)/m - g,"acc z" );

        //window2.addSubplot( vx,"DifferentialState vx" );
        window2.addSubplot( -Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*sin(theta)/m,"acc x" );
//        window2.addSubplot( acc(0),"acc x libquadspline" );

        //window2.addSubplot( vy,"DifferentialState vy" );
        //window2.addSubplot( Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*sin(phi)*cos(theta)/m,"acc y" );
//        window2.addSubplot( u1,"u1" );
//        window2.addSubplot( u2,"u2" );
//        window2.addSubplot( u3,"u3" );
//        window2.addSubplot( u4,"u4" );
        //window2.addSubplot( u4+u1+u2+u3,"f tot" );
        //window2.addSubplot( vu1,"vu1" );
       //window2.addSubplot( vu2,"vu2" );
        //window2.addSubplot( vu3,"vu3" );
        //window2.addSubplot( vu4,"vu4" );
        //window2.addSubplot( phi,"phi" );
//        window2.addSubplot( q,"theta" );
//        window2.addSubplot( x_init(10),"theta" );
//        window2.addSubplot(x_init(7), "pith libquad");
//        window2.addSubplot( sqrt((y-yf)*(y-yf)),"|y-yf|" );
        //window2.addSubplot( y,"y" );

//        window2.addSubplot( PLOT_KKT_TOLERANCE,"KKT Tolerance" );


    // ------------------------- INITIALIZE OCP PROBLEM ------------------------------- //
       OptimizationAlgorithm algorithm(ocp);
#if INIT==0 || INIT==1
       algorithm.initializeDifferentialStates(x_init);
#endif

    algorithm.set(DISCRETIZATION_TYPE, SINGLE_SHOOTING);
//    algorithm.set( DYNAMIC_SENSITIVITY,  FORWARD_SENSITIVITY );

//        algorithm.set( HESSIAN_APPROXIMATION, EXACT_HESSIAN );
    algorithm.set( MAX_NUM_ITERATIONS, 100 );
    algorithm.set( KKT_TOLERANCE, 1e-20 );
// 	algorithm.set( MAX_NUM_INTEGRATOR_STEPS, 4 );

//    algorithm << window2;
        LogRecord logdiff(LOG_AT_END , PS_DEFAULT);
    logdiff << LOG_DIFFERENTIAL_STATES;
//logdiff << LOG_KKT_TOLERANCE;
        algorithm << logdiff;

        t = clock();

    algorithm.solve();

    t = clock() - t;
    std::cout << "total time : " << (((float)t)/CLOCKS_PER_SEC)<<std::endl;

    algorithm.getLogRecord(logdiff);
    std::ofstream file;
    file.open("/tmp/log.txt",std::ios::out);
    logdiff.print(file);

    // -------------------- Export Results ------------------ //
//    VariablesGrid states, controls;
//    algorithm.getDifferentialStates(x_init);
//    x_init.print();

//    algorithm.getDifferentialStates(states);
//    algorithm.getControls("/tmp/controls.txt");
//    algorithm.getControls(controls);
#endif

#if ALGO==2
    // ------------------------------- CODe GENERATION TOOL -------------- //
//    OCPexport mpc( ocp );

//    mpc.set( HESSIAN_APPROXIMATION,       GAUSS_NEWTON    );
//    mpc.set( DISCRETIZATION_TYPE,         SINGLE_SHOOTING );
//    mpc.set( INTEGRATOR_TYPE,             INT_RK4         );
//    mpc.set( NUM_INTEGRATOR_STEPS,        30              );

//    mpc.set( QP_SOLVER,                   QP_QPOASES      );
//    mpc.set( HOTSTART_QP,                 YES             );
//// 	mpc.set( LEVENBERG_MARQUARDT,         1.0e-4          );
//    mpc.set( GENERATE_TEST_FILE,          YES             );
//    mpc.set( GENERATE_MAKE_FILE,          YES             );
////    mpc.set( GENERATE_MATLAB_INTERFACE,   YES             );
////    mpc.set( GENERATE_SIMULINK_INTERFACE, YES             );
//// 	mpc.set( USE_SINGLE_PRECISION,        YES             );

//    if (mpc.exportCode( "getting_started_export" ) != SUCCESSFUL_RETURN)
//        exit( EXIT_FAILURE );

//    mpc.printDimensionsQP( );
#endif

    return 0;
}
/* <<< end tutorial code <<< */



//  algorithm.set( DISCRETIZATION_TYPE, MULTIPLE_SHOOTING );
//  algorithm.set( DISCRETIZATION_TYPE, SINGLE_SHOOTING   );
//
//  algorithm.set( DYNAMIC_SENSITIVITY, FORWARD_SENSITIVITY );
//  algorithm.set( DYNAMIC_SENSITIVITY, BACKWARD_SENSITIVITY );

//  algorithm.set( INTEGRATOR_TYPE, INT_RK45 );
//  algorithm.set( INTEGRATOR_TYPE, INT_RK78 );
//  algorithm.set( INTEGRATOR_TYPE, INT_BDF );
//
//  algorithm.set( KKT_TOLERANCE, 1e-4 );
//  algorithm.set( MAX_NUM_ITERATIONS, 20 );
