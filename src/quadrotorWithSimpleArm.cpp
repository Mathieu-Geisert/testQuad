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

//#include <acado_optimal_control.hpp>
//#include <acado_gnuplot.hpp>
#include <acado_toolkit.hpp>
#include <include/acado_gnuplot/gnuplot_window.hpp>

#include <time.h>
#include <math.h>

/* >>> start tutorial code >>> */
int main( ){

    USING_NAMESPACE_ACADO

    // INTRODUCE THE VARIABLES:
    // -------------------------
    DifferentialState     x,y,z,vx,vy,vz,phi,theta,psi,p,q,r,u1,u2,u3,u4,arm1x,arm1v,arm2x,arm2v;
    Control               vu1,vu2,vu3,vu4;
    //Parameter		  T;
    DifferentialEquation  f(0.0,15);

    const double c = 0.00001;
    const double Cf = 0.00065;
    const double d = 0.250;
    const double Jx = 0.018;
    const double Jy = 0.018;
    const double Jz = 0.026;
    const double m = 0.9;
    const double g = 9.81;
    const double Cx = 0.1;

    const double Iarm1 = 0 ;

    const double xf = 20.;
    const double yf = 0.;
    const double zf = 0.;

    double coeffU = 0.001;
    double coeffX = .01;

    // DEFINE A DIFFERENTIAL EQUATION:
    // -------------------------------
    OCP ocp( 0.0, 15, 60 );
    //ocp.minimizeMayerTerm( 1000*coeffX*((x-xf)*(x-xf) + (y-yf)*(y-yf) + (z-zf)*(z-zf)) );
    ocp.minimizeLagrangeTerm(coeffX*(sqrt(1+(x-xf)*(x-xf)) + sqrt(1+(y-yf)*(y-yf)) + sqrt(1+(z-zf)*(z-zf))) + coeffU*( sqrt(1+(u1-58.27)*(u1-58.27))+sqrt(1+(u2-58.27)*(u2-58.27))+sqrt(1+(u3-58.27)*(u3-58.27))+sqrt(1+(u4-58.27)*(u4-58.27))));

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

    //Start conditions
    ocp.subjectTo( AT_START, x ==  0 );
    ocp.subjectTo( AT_START, y ==  0 );
    ocp.subjectTo( AT_START, z ==  0 );
    ocp.subjectTo( AT_START, vx ==  0 );
    ocp.subjectTo( AT_START, vy ==  0 );
    ocp.subjectTo( AT_START, vz ==  0 );
    ocp.subjectTo( AT_START, phi ==  0.0 );
    ocp.subjectTo( AT_START, theta ==  -0. );
    ocp.subjectTo( AT_START, psi ==  -0.);
    ocp.subjectTo( AT_START, p ==  0.0 );
    ocp.subjectTo( AT_START, q ==  0.0 );
    ocp.subjectTo( AT_START, r ==  0.0 );
    ocp.subjectTo( AT_START, (d*Cf*(u1*u1-u2*u2)+(Jy-Jz)*q*r)/Jx ==  0.0 );
    ocp.subjectTo( AT_START, (d*Cf*(u4*u4-u3*u3)+(Jz-Jx)*p*r)/Jy ==  0.0 );
    ocp.subjectTo( AT_START, (c*(u1*u1+u2*u2-u3*u3-u4*u4)+(Jx-Jy)*p*q)/Jz ==  0.0 );
    ocp.subjectTo( AT_START, -Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*sin(theta)/m ==  0 );
    ocp.subjectTo( AT_START, Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*sin(phi)*cos(theta)/m ==  0 );
    ocp.subjectTo( AT_START, Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*cos(phi)*cos(theta)/m - g ==  0 );

    //Command constraints
    //on each motor
    //ocp.subjectTo( (u1*u1+u2*u2-u3*u3-u4*u4) == 0);
    ocp.subjectTo( 16 <= u1 <= 95 );
    ocp.subjectTo( 16 <= u2 <= 95 );
    ocp.subjectTo( 16 <= u3 <= 95 );
    ocp.subjectTo( 16 <= u4 <= 95 );
    ocp.subjectTo( -100 <= vu1 <= 100 );
    ocp.subjectTo( -100 <= vu2 <= 100 );
    ocp.subjectTo( -100 <= vu3 <= 100 );
    ocp.subjectTo( -100 <= vu4 <= 100 );
    //ocp.subjectTo( -15 <= vx <= 15);
    //ocp.subjectTo( -0.1 <= z <= 0.1);
    //ocp.subjectTo( -0.1 <= y <= 0.1);

    //on total power
    //ocp.subjectTo( u1*u1+u2*u2+u3*u3+u4*u4 <= 800 );


    // DEFINE A PLOT WINDOW:
    // ---------------------
    GnuplotWindow window2(PLOT_AT_EACH_ITERATION);
        window2.addSubplot( z,"DifferentialState z" );
        window2.addSubplot( x,"DifferentialState x" );
        //window2.addSubplot( y,"DifferentialState y" );
        window2.addSubplot( vx,"DifferentialState vx" );
        window2.addSubplot( vz,"DifferentialState vz" );
        window2.addSubplot( Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*cos(phi)*cos(theta)/m - g,"acc z" );
        //window2.addSubplot( vx,"DifferentialState vx" );
        window2.addSubplot( -Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*sin(theta)/m,"acc x" );
        //window2.addSubplot( vy,"DifferentialState vy" );
        //window2.addSubplot( Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*sin(phi)*cos(theta)/m,"acc y" );
        window2.addSubplot( u1,"u1" );
        window2.addSubplot( u2,"u2" );
        window2.addSubplot( u3,"u3" );
        window2.addSubplot( u4,"u4" );
        //window2.addSubplot( u4+u1+u2+u3,"f tot" );
        //window2.addSubplot( vu1,"vu1" );
       //window2.addSubplot( vu2,"vu2" );
        //window2.addSubplot( vu3,"vu3" );
        //window2.addSubplot( vu4,"vu4" );
        //window2.addSubplot( phi,"phi" );
        //window2.addSubplot( theta,"theta" );
        //window2.addSubplot( z,"z" );
        //window2.addSubplot( y,"y" );


        //window.addSubplot( PLOT_KKT_TOLERANCE,"KKT Tolerance" );
//         window.addSubplot( 0.5 * m * v*v,"Kinetic Energy" );


    // DEFINE AN OPTIMIZATION ALGORITHM AND SOLVE THE OCP:
    // ---------------------------------------------------
    OptimizationAlgorithm algorithm(ocp);
    //algorithm.set(DISCRETIZATION_TYPE, SINGLE_SHOOTING);
    //VariableGrid T_init(1
   // algorithm.initializeParameters(4.0);

// 	algorithm.set( INTEGRATOR_TYPE, INT_BDF );

//     algorithm.set( INTEGRATOR_TOLERANCE, 1e-6 );
//     algorithm.set( KKT_TOLERANCE, 1e-3 );

  //algorithm.set( DYNAMIC_SENSITIVITY,  FORWARD_SENSITIVITY );
    Grid timeGrid(0.0,1.0,41);
    VariablesGrid u_init(16, timeGrid);
    for (int i = 0 ; i<41 ; i++ ) {
      if(i>24 && i<33) {
      //srand (time(NULL));
      u_init(i,12) = 50;
      u_init(i,15) = 35;
      u_init(i,14) = 50;
      u_init(i,13) = 50;
      }/*
      if(i>15 && i<25) {u_init(i,12) = 15.0;
      u_init(i,15) = 90.0;
      u_init(i,14) = 15.0;
      u_init(i,13) = 15.0;}
      if(i>24 && i<34) {u_init(i,12) = 90.0;
      u_init(i,15) = 90.0;
      u_init(i,14) = 90.0;
      u_init(i,13) = 90.0;}
      if(i>34 && i<41) {u_init(i,12) = 90.0;
      u_init(i,15) = 15.0;
      u_init(i,14) = 15.0;
      u_init(i,13) = 15.0;}*/
      //u_init(i,8) = i*1.6/40;
    }
    //u_init.print();
    //algorithm.initializeDifferentialStates(u_init);
    //algorithm.set( HESSIAN_APPROXIMATION, EXACT_HESSIAN );
    algorithm.set( MAX_NUM_ITERATIONS, 2000 );
    algorithm.set( KKT_TOLERANCE, 1e-18 );
// 	algorithm.set( MAX_NUM_INTEGRATOR_STEPS, 4 );

    //algorithm << window3;
    algorithm << window2;
    algorithm.solve();


    VariablesGrid states, controls;
    algorithm.getDifferentialStates("/tmp/states.txt");
    algorithm.getDifferentialStates(states);
    algorithm.getControls("/tmp/controls.txt");
    algorithm.getControls(controls);

    /*VariablesGrid output(5,timeGrid);
    for ( int i =0 ; i<41 ; i++ ) {
      output(i,0) = states(i,0);
      output(i,1) = states(i,3);
     if(i>0) output(i,2) = (states(i,3)-states(i-1,3))/0.025;
     if(i>1) output(i,3) = (output(i,2)-output(i-1,2));///0.025;
     if(i>2) output(i,4) = (output(i,3)-output(i-1,3));///0.025;
    }*/
    //output.print();
    //states.print();

// 	BlockMatrix sens;
// 	algorithm.getSensitivitiesX( sens );
// 	sens.print();

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

//  algorithm.set( PRINT_SCP_METHOD_PROFILE, YES );

