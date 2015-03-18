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
//#include <acado_controller.hpp>
#include <acado_toolkit.hpp>
#include <acado_gnuplot.hpp>


/* >>> start tutorial code >>> */
int main( ){

    USING_NAMESPACE_ACADO

    // INTRODUCE THE VARIABLES:
    // -------------------------
    DifferentialState     x,y,z,vx,vy,vz,phi,theta,psi,p,q,r;//,u1,u2,u3,u4;
    Control               vu1,vu2,vu3,vu4;

    //Disturbance R;
    //Parameter		  T;
    DifferentialEquation  f;
    OutputFcn out;

    const double c = 8.004e-4;
    const double d = 0.315;
    const double Jx = 0.0820;
    const double Jy = 0.0845;
    const double Jz = 0.1377;
    const double m = 4.34;
    const double g = 9.81;
    const double Cx = 0.1;


    // DEFINE A DIFFERENTIAL EQUATION:
    // -------------------------------

    f << dot(x) == vx;
    f << dot(y) == vy;
    f << dot(z) == vz;
    f << dot(vx) == -(vu1+vu2+vu3+vu4)*sin(theta)/m;
    f << dot(vy) == (vu1+vu2+vu3+vu4)*sin(phi)*cos(theta)/m;
    f << dot(vz) == (vu1+vu2+vu3+vu4)*cos(phi)*cos(theta)/m - g;
    f << dot(phi) == p + sin(theta)*r;
    f << dot(theta) == cos(phi)*q-sin(phi)*cos(theta)*r;
    f << dot(psi) == sin(phi)*q+cos(phi)*cos(theta)*r;
    f << dot(p) == (d*(vu1-vu2)+(Jy-Jz)*q*r)/Jx;
    f << dot(q) == (d*(vu4-vu3)+(Jz-Jx)*p*r)/Jy;
    f << dot(r) == (c*(vu1+vu2-vu3-vu4)+(Jx-Jy)*p*q)/Jz;
    //f << dot(u1) == vu1;
    //f << dot(u2) == vu2;
    //f << dot(u3) == vu3;
    //f << dot(u4) == vu4;

    DynamicSystem dynSys(f,out);
 
    Process myProcess;
    myProcess.setDynamicSystem(dynSys, INT_RK45);
    myProcess.set(ABSOLUTE_TOLERANCE,1.0e-8);
 
    DVector x0(12);
    x0.setZero();
    x0(0) = 10.0;

    //myProcess.initializeStartValues(x0);
    //myProcess.set(PLOT_RESOLUTION,HIGH);


    // DEFINE MPC CONTROL
    OCP ocp(0.0,1.0,20);
     
    //define cost function
    Function h,i;
    h << x;
    h << y;
    h << z;
    h << vu1;
    h << vu2;
    h << vu3;
    h << vu4;
    i << x;
    i << y;
    i << z;

    DMatrix Q(7,7),P(3,3);
    Q(0,0) = 10;
    Q(1,1) = 10;
    Q(2,2) = 10;
    Q(3,3) = 1;
    Q(4,4) = 1;
    Q(5,5) = 1;
    Q(6,6) = 1;
    P(0,0) = 25;
    P(1,1) = 25;
    P(2,2) = 25;

    DVector v(7),w(3);
    v.setAll(0.0);
    v(0) = 0.0;
    w.setAll(0.0);
    w(0) = 0.0;

    ocp.minimizeLSQ(Q,h,v);
    //ocp.minimizeLSQEndTerm(P,i,w);
    ocp.subjectTo(f);
    /*ocp.subjectTo( 5 <= u1 <= 40 );
    ocp.subjectTo( 5 <= u2 <= 40 );
    ocp.subjectTo( 5 <= u3 <= 40 );
    ocp.subjectTo( 5 <= u4 <= 40 );
    ocp.subjectTo( -4000 <= vu1 <= 4000 );    
    ocp.subjectTo( -4000 <= vu2 <= 4000 );
    ocp.subjectTo( -4000 <= vu3 <= 4000 );
    ocp.subjectTo( -4000 <= vu4 <= 4000 );*/

    //SETTING THE REAL-TIME ALGORITHM:
    RealTimeAlgorithm alg(ocp,0.025);
    alg.set(MAX_NUM_ITERATIONS,500);
    //alg.set(DISCRETIZATION_TYPE, SINGLE_SHOOTING );
    alg.set(KKT_TOLERANCE, 10e-4 );
    //alg.set(PLOT_RESOLUTION,HIGH);
    StaticReferenceTrajectory zeroReference;
    Controller controller(alg);//, zeroReference);

    //controller.init(0.0, x0);
    //controller.step(0.0,x0);
    //SETTING UP THE SIMULATION ENVIRONMENT
    //SimulationEnvironment sim(0.0,15,myProcess,controller);
    //sim.init(x0);

    // DEFINE A PLOT WINDOW:
    // ---------------------
    GnuplotWindow window2(PLOT_AT_START);
        window2.addSubplot( z,"DifferentialState z" );
        window2.addSubplot( x,"DifferentialState x" );
        window2.addSubplot( y,"DifferentialState y" );
        window2.addSubplot( vx,"DifferentialState vx" );
        //window2.addSubplot( vz,"DifferentialState vz" );
        //window2.addSubplot( (u1+u2+u3+u4)*cos(phi)*cos(theta)/m - g,"acc z" );
        //window2.addSubplot( -(u1+u2+u3+u4)*sin(theta)/m,"acc x" );
        //window2.addSubplot( vy,"DifferentialState vy" );
        //window2.addSubplot( (u1+u2+u3+u4)*sin(phi)*cos(theta)/m,"acc y" );
        //window2.addSubplot( u1,"u1" );
        //window2.addSubplot( u2,"u2" );
        //window2.addSubplot( u3,"u3" );
        //window2.addSubplot( u4,"u4" );
        //window2.addSubplot( u4+u1+u2+u3,"f tot" );
        //window2.addSubplot( vu1,"vu1" );
       //window2.addSubplot( vu2,"vu2" );
        //window2.addSubplot( vu3,"vu3" );
        //window2.addSubplot( vu4,"vu4" );
        //window2.addSubplot( phi,"phi" );
        //window2.addSubplot( theta,"theta" );
        //window2.addSubplot( psi,"psi" );


    //alg << window2;
    controller.init(0.0, x0);
    controller.step(0.0,x0);
    //sim.run();  
//    VariablesGrid u(4,0.0,10.0,10);
//    u.setAll(10.70);
//    u.print();

//    myProcess.init(0.0);
//    myProcess.run(u);

    /*VariablesGrid xSim;
    myProcess.getLast( LOG_SIMULATED_DIFFERENTIAL_STATES,xSim );
    window2.addSubplot(xSim(0),"x");
    window2.addSubplot(xSim(1),"y");
    window2.addSubplot(xSim(2),"z");
    window2.plot();*/
    return 0;
}
/* <<< end tutorial code <<< */

