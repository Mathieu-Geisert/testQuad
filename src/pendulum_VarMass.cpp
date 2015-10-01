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


#define ALGO 0 // OCP; MPC; OCP then MPC; GENE
#define INIT 1 //libquadspline, quasistatic, quasistatic with avoidance or static
#define OBS 0
#define PLOT 1 //plots, file
#define AVOID_SINGULARITIES 1 //No, yes

// PENDULE AVEC POIDS, COORDONNEES Euler(3,2,1)

int main( ){

    std::ofstream file;
    file.open("/tmp/log.txt",std::ios::out);

    USING_NAMESPACE_ACADO

    // INTRODUCE THE VARIABLES:
    // -------------------------
    DifferentialState     x,y,z,vx,vy,vz,phi,theta,psi,p,q,r,u1,u2,u3,u4, pendRoll, pendPitch, Wpendx, Wpendy, mtot;
    Control               vu1,vu2,vu3,vu4;
//    Parameter		  T1,T2;
    DifferentialEquation  f;
    Disturbance D;


    // Quad constants
    const double c = 0.00001;
    const double Cf = 0.00065;
    const double d = 0.250;
    const double Jx = 0.018;
    const double Jy = 0.018;
    const double Jz = 0.026;
    const double mq = 0.9;
    const double mp = 0.45;
    const double m = mq + mp;
    const double Lpend =4.;
    const double Xqg = mp/m*Lpend;
    const double Xpg = mq/m*Lpend;
    const double Ipend = Xqg*Xqg*mq + Xpg*Xpg*mp;
    const double g = 9.81;
    const double Cx = 0.1;

    // Minimization Weights
    double coeffU = 0.00000000000000000000001;
    double coeffX = .00001;
    double coeffX2 = coeffX * 1.;
    double coeffX3 = coeffX * 0.00001;
    double coeffO = -coeffX * 0.1;

    // final values
    double xf = 13., yf = 0., zf = 2.;

    // Temps max
    double tps = 8.;
    int nb_it = 60.;
    double tmpc = 0.2;
    int mpc_nb_it = 1;
    int ocp_nb_it = 250;

    DVector x0(16);
    x0.setZero();
    x0(0) = 0.;
    x0(12) = 58.;
    x0(13) = 58.;
    x0(14) = 58.;
    x0(15) = 58.;

    // DEFINE A DIFFERENTIAL EQUATION:
    // -------------------------------
    OCP ocp( 0.0, tps, nb_it );

    //  ------------------------- LSQ MINIMISATION ------------------------ //
    Function h, hf;
    h << x;
//    h << y;
//    h << z;

//    h << vu1;
//    h << vu2;
//    h << vu3;
//    h << vu4;

//    h << p;
//    h << q;
//    h << r;

//    h << phi;
//    h << theta;

//    h << exp((16 - (x*x+(z-5)*(z-5)))*0.01);

    hf << x;
    hf << y;
    hf << z;

    DMatrix Q(1,1), Qf(3,3);
    Q(0,0) = coeffX;
//    Q(1,1) = coeffX;
//    Q(2,2) = coeffX;

//    Q(3,3) = coeffU;
//    Q(4,4) = coeffU;
//    Q(5,5) = coeffU;
//    Q(6,6) = coeffU;

//    Q(7,7) = coeffX2;
//    Q(8,8) = coeffX2;
//    Q(9,9) = coeffX2;

//    Q(10,10) = coeffX3;
//    Q(11,11) = coeffX3;

//    Q(12,12) = coeffX*0.03;

    Qf(0,0) = coeffX;//*100;
    Qf(1,1) = coeffX;//*100;
    Qf(2,2) = coeffX;//*100;

    DVector reff(3), ref(1);
    ref(0) = xf;
//    ref(1) = yf;
//    ref(2) = zf;
//    ref(3) = 0.;
//    ref(4) = 0.;
//    ref(5) = 0.;
//    ref(6) = 0.;
//    ref(7) = 0.;
//    ref(8) = 0.;
//    ref(9) = 0.;
//    ref(10) = 0.;
//    ref(11) = 0.;
//    ref(10) = 0.;

    reff(0) = xf;
    reff(1) = yf;
    reff(2) = zf;


//    ocp.minimizeLSQ ( Q, h, ref);
//    ocp.minimizeLSQEndTerm(Qf, hf, reff);

//    ocp.minimizeMayerTerm( (x-xf)*(x-xf) + (y-yf)*(y-yf) + (z-zf)*(z-zf));
        Grid timeGrid2(0.0,tps,nb_it+1);
        VariablesGrid Weight(1,timeGrid2);
        for (int i = 0 ; i<nb_it+1 ; i++ ) {
            if (i > 20 && i <50 ) { Weight(i,0) = 1.;}
            else { Weight(i,0) = 0.;}
        }

        TIME ti;
//        ocp.subjectTo( 0 <= T1 <= 8);
//        ocp.subjectTo( 0 <= T2 <= 8);
//        ocp.subjectTo( 0 <= T2-T1 );

        ocp.minimizeLagrangeTerm( exp(-(ti-3.)*(ti-3.)/0.25)*(((x-sin(pendRoll)*cos(pendPitch)*4.)-5.)*((x-sin(pendRoll)*cos(pendPitch)*4.)-5.) + (y+sin(pendPitch)*4.)*(y+sin(pendPitch)*4.) + (z-cos(pendPitch)*cos(pendRoll)*4.)*(z-cos(pendPitch)*cos(pendRoll)*4.) + 1.*((vx - Lpend*mq/mtot*(Wpendx*cos(pendRoll)-Wpendy*sin(pendRoll)*sin(pendPitch)))*(vx - Lpend*mq/mtot*(Wpendx*cos(pendRoll)-Wpendy*sin(pendRoll)*sin(pendPitch))) + (vy + Lpend*mq/mtot*(Wpendy*cos(pendPitch)))*(vy + Lpend*mq/mtot*(Wpendy*cos(pendPitch))) + (vz + Lpend*mq/mtot*(Wpendx*sin(pendRoll)+Wpendy*sin(pendPitch)*cos(pendRoll)))*(vz + Lpend*mq/mtot*(Wpendx*sin(pendRoll)+Wpendy*sin(pendPitch)*cos(pendRoll))))) + exp(-(ti-6.)*(ti-6.)/0.25)*(((x-sin(pendRoll)*cos(pendPitch)*4.)-10.)*((x-sin(pendRoll)*cos(pendPitch)*4.)-10.) + (y+sin(pendPitch)*4.)*(y+sin(pendPitch)*4.) + (z-cos(pendPitch)*cos(pendRoll)*4.)*(z-cos(pendPitch)*cos(pendRoll)*4.) + 1.*((vx - Lpend*mq/mtot*(Wpendx*cos(pendRoll)-Wpendy*sin(pendRoll)*sin(pendPitch)))*(vx - Lpend*mq/mtot*(Wpendx*cos(pendRoll)-Wpendy*sin(pendRoll)*sin(pendPitch))) + (vy + Lpend*mq/mtot*(Wpendy*cos(pendPitch)))*(vy + Lpend*mq/mtot*(Wpendy*cos(pendPitch))) + (vz + Lpend*mq/mtot*(Wpendx*sin(pendRoll)+Wpendy*sin(pendPitch)*cos(pendRoll)))*(vz + Lpend*mq/mtot*(Wpendx*sin(pendRoll)+Wpendy*sin(pendPitch)*cos(pendRoll))))) );
//    ocp.minimizeLagrangeTerm( (log(1+((x-sin(pendRoll)*cos(pendPitch)*Xpg)-5)*((x-sin(pendRoll)*cos(pendPitch)*Xpg)-5)) + log(1+(y+sin(pendPitch)*Xpg)*(y+sin(pendPitch)*Xpg)) + log(1+(z-cos(pendPitch)*cos(pendRoll)*Xpg)*(z-cos(pendPitch)*cos(pendRoll)*Xpg))) );
//      ocp.minimizeLagrangeTerm( x*x+y*y+z*z );

//    ocp.minimizeLagrangeTerm(coeffX*(sqrt(1+(x-xf)*(x-xf)) + sqrt(1+(y-yf)*(y-yf)) + sqrt(1+(z-zf)*(z-zf))));
//    ocp.minimizeLagrangeTerm(coeffX*(sqrt(1+(x-ref(0))*(x-ref(0))) + sqrt(1+(y-ref(1))*(y-ref(1))) + sqrt(1+(z-ref(2))*(z-ref(2))) + (r-ref(7))*(r-ref(7)) + (p-ref(8))*(p-ref(8)) + (q-ref(9))*(q-ref(9)) ) + coeffU*( sqrt(1+(vu1-ref(3))*(vu1-ref(3)))+sqrt(1+(vu2-ref(4))*(vu2-ref(4)))+sqrt(1+(vu3-ref(5))*(vu3-ref(5)))+sqrt(1+(vu4-ref(6))*(vu4-ref(6)))));
//    ocp.minimizeLagrangeTerm(((x-xf)*(x-xf) + (y-yf)*(y-yf) + (z-zf)*(z-zf))*coeffX);

    // ------------------------------------------------------------------ //


        f << dot(x) == vx + Lpend/mtot*(mtot-mq)*(Wpendx*cos(pendRoll)-Wpendy*sin(pendRoll)*sin(pendPitch) /*+ 2*4*(exp(-(ti-T1)*(ti-T1)/0.025) - exp(-(ti-T2)*(ti-T2)/0.025) )/m*sin(pendRoll)*cos(pendRoll)*/);
        f << dot(y) == vy - Lpend/mtot*(mtot-mq)*(Wpendy*cos(pendPitch) /*+ 2*4*(exp(-(ti-T1)*(ti-T1)/0.025) - exp(-(ti-T2)*(ti-T2)/0.025) )/m*sin(pendPitch)*/) ;
        f << dot(z) == vz - Lpend/mtot*(mtot-mq)*(Wpendx*sin(pendRoll)+Wpendy*sin(pendPitch)*cos(pendRoll) /*- 2*4*(exp(-(ti-T1)*(ti-T1)/0.025) - exp(-(ti-T2)*(ti-T2)/0.025) )/m*cos(pendPitch)*cos(pendRoll)*/);
        f << dot(vx) == Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*sin(theta)/mtot + mq*Lpend*D/(mtot*(mtot+D))*(sin(pendPitch)*sin(pendRoll)*Wpendy-cos(pendRoll)*Wpendx);
        f << dot(vy) == - Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*sin(psi)*cos(theta)/mtot + mq*Lpend*D/(mtot*(mtot+D))*(cos(pendPitch)*Wpendy);
        f << dot(vz) == Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*cos(psi)*cos(theta)/mtot - g + mq*Lpend*D/(mtot*(mtot+D))*(sin(pendPitch)*cos(pendRoll)*Wpendy + sin(pendRoll)*Wpendx);
    f << dot(phi) == -cos(phi)*tan(theta)*p+sin(phi)*tan(theta)*q+r;
    f << dot(theta) == sin(phi)*p+cos(phi)*q;
    f << dot(psi) == cos(phi)/cos(theta)*p-sin(phi)/cos(theta)*q;
    f << dot(p) == (d*Cf*(u1*u1-u2*u2)+(Jy-Jz)*q*r)/Jx;
    f << dot(q) == (d*Cf*(u4*u4-u3*u3)+(Jz-Jx)*p*r)/Jy;
    f << dot(r) == (c*(u1*u1+u2*u2-u3*u3-u4*u4)+(Jx-Jy)*p*q)/Jz;
    f << dot(u1) == vu1;
    f << dot(u2) == vu2;
    f << dot(u3) == vu3;
    f << dot(u4) == vu4;
    f << dot(pendRoll) == Wpendx/cos(pendPitch);
    f << dot(pendPitch) == Wpendy;
    f << dot(Wpendx) == (cos(pendRoll)*Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*sin(theta) - sin(pendRoll)*Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*cos(psi)*cos(theta))/(mq*Lpend);
    f << dot(Wpendy) ==  - ( - cos(pendPitch)*Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*sin(psi)*cos(theta) + sin(pendPitch)*sin(pendRoll)*Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*sin(theta) + sin(pendPitch)*cos(pendRoll)*Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*cos(psi)*cos(theta) ) / (mq*Lpend);
    f << dot(mtot) == D;


    // DEFINE AN OPTIMAL CONTROL PROBLEM:
    // ----------------------------------

    //Dynamic
    ocp.subjectTo( f );

// ------------------------------------------------------------------------------------ //



    // ----------------------------- Set initial guess -------------------------- //



#if INIT==0 || INIT==1 || INIT==2 || INIT==3
    Grid timeGrid(0.0,tps,nb_it+1);
    VariablesGrid x_init(21, timeGrid);
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
        x_init(i,16) = 0.;
        x_init(i,17) = 0.;
        x_init(i,18) = 0.;
        x_init(i,19) = 0.;
        acc(i,0) = conf.getAx();
        acc(i,1) = conf.getAy();
        acc(i,2) = conf.getAz();
    }
#endif

#if INIT==1      // init with quasistatic solution
    for (int i = 0 ; i<nb_it+1 ; i++ ) {
//        if (i<<50)
            x_init(i,0) = xf*double(i)/double(nb_it);;
//        else
//            x_init(i,0) = 10 - (double(i)-50);
        x_init(i,1) = yf*double(i)/double(nb_it);
        x_init(i,2) = 2.;//zf*double(i)/double(nb_it);
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
        x_init(i,16) = 0.;
        x_init(i,17) = 0.;
        x_init(i,18) = 0.;
        x_init(i,19) = 0.;
        x_init(i,20) = .9;
    }
#endif

#if INIT==2      // init with quasistatic solution and obtacle avoidance
    for (int i = 0 ; i<nb_it+1 ; i++ ) {
        if (i <= nb_it/4) {
        x_init(i,0) = 3.*double(i)/double(nb_it/4);
        x_init(i,1) = 0.*double(i)/double(nb_it/4);
        x_init(i,2) = 3.*double(i)/double(nb_it/4);
        }
        if (i > nb_it/4 && i <= nb_it/2) {
        x_init(i,0) = 3. + -3.*(double(i)-nb_it/4)/double(nb_it/4);
        x_init(i,1) = 0.*(double(i)-nb_it/4)/double(nb_it/4);
        x_init(i,2) = 3. + 3.*(double(i)-nb_it/4)/double(nb_it/4);
        }
        if (i > nb_it/2) {
        x_init(i,0) = 0. + 1.*(double(i)-2*nb_it/4)/double(nb_it/2);
        x_init(i,1) = 0.*(double(i)-2*nb_it/4)/double(nb_it/2);
        x_init(i,2) = 6. + 9.*(double(i)-2*nb_it/4)/double(nb_it/2);
        }
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
        x_init(i,16) = 0.;
        x_init(i,17) = 0.;
        x_init(i,18) = 0.;
        x_init(i,19) = 0.;
    }
#endif

#if INIT==3      // init with static
    for (int i = 0 ; i<nb_it+1 ; i++ ) {
        x_init(i,0) = x0(0);
        x_init(i,1) = x0(1);
        x_init(i,2) = 2.;//x0(2);
        x_init(i,3) = 0.;
        x_init(i,4) = 0.;
        x_init(i,5) = 0.;
        x_init(i,6) = 0.;
        x_init(i,7) = 0.;
        x_init(i,8) = 0.;
        x_init(i,9) = 0.;
        x_init(i,10) = 0.;
        x_init(i,11) = 0.;
        x_init(i,12) = 58.;
        x_init(i,13) = 58.;
        x_init(i,14) = 58.;
        x_init(i,15) = 58.;
        x_init(i,16) = 0.;
        x_init(i,17) = 0.;
        x_init(i,18) = 0.;
        x_init(i,19) = 0.;
                x_init(i,20) = 1.35;
    }
#endif

//    x_init.print();
//    GnuplotWindow windowXinit;
//    windowXinit.addSubplot( x_init(0), "x " );
//    windowXinit.addSubplot( x_init(1), "y " );
//    windowXinit.addSubplot( x_init(2), "z " );
//    windowXinit.plot();
 // ----------------------------------------------------------------------------------- //

    // ---------------------------- DEFINE CONTRAINTES --------------------------------- //


    // Min/Max time
//    ocp.subjectTo(1 <= T <= 10);

#if ALGO==0 //|| ALGO==1
    //Start conditions
    ocp.subjectTo( AT_START, x == 0 );
    ocp.subjectTo( AT_START, y == 0  );
    ocp.subjectTo( AT_START, z ==  2. );
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
    ocp.subjectTo( AT_START, Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*sin(psi)*cos(theta)/m ==   0 );
    ocp.subjectTo( AT_START, Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*cos(psi)*cos(theta)/m - g ==  0 );
    ocp.subjectTo( AT_START, Wpendx == 0. );
    ocp.subjectTo( AT_START, Wpendy == 0. );
    ocp.subjectTo( AT_START, pendPitch == 0.0 );
    ocp.subjectTo( AT_START, pendRoll == 0. );
        ocp.subjectTo( AT_START, mtot == 0.9 );
#endif

    //end conditions
//    ocp.subjectTo( AT_END  , (x-xf)*(x-xf) <= 100);//(x_init(nb_it,0)-xf)*(x_init(nb_it,0)-xf) );
//    ocp.subjectTo( AT_END  , (y-yf)*(y-yf) <= 10);//(x_init(nb_it,1)-yf)*(x_init(nb_it,1)-yf) );
//    ocp.subjectTo( AT_END  , (z-zf)*(z-zf) <= 10);//(x_init(nb_it,2)-zf)*(x_init(nb_it,2)-zf) );

    ocp.subjectTo( AT_END  , x ==  xf );
    ocp.subjectTo( AT_END  , y ==  yf );
    ocp.subjectTo( AT_END  , z ==  zf );
    ocp.subjectTo( AT_END  , vx ==  0 );
    ocp.subjectTo( AT_END  , vy ==   0);
    ocp.subjectTo( AT_END  , vz ==  0 );
////    ocp.subjectTo( AT_END  , phi ==  0.0 );
////    ocp.subjectTo( AT_END  , theta ==  -0. );
////    ocp.subjectTo( AT_END  , psi == -0.0 );
    ocp.subjectTo( AT_END  , p ==  0.0 );
    ocp.subjectTo( AT_END  , q ==  0.0 );
    ocp.subjectTo( AT_END  , r ==  0.0 );
    ocp.subjectTo( AT_END, (d*Cf*(u1*u1-u2*u2)+(Jy-Jz)*q*r)/Jx ==  0.0 );
    ocp.subjectTo( AT_END, (d*Cf*(u4*u4-u3*u3)+(Jz-Jx)*p*r)/Jy ==  0.0 );
    ocp.subjectTo( AT_END, (c*(u1*u1+u2*u2-u3*u3-u4*u4)+(Jx-Jy)*p*q)/Jz ==  0.0);
    ocp.subjectTo( AT_END, -Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*sin(theta)/m ==  0 );
    ocp.subjectTo( AT_END, Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*sin(psi)*cos(theta)/m == 0.  );
    ocp.subjectTo( AT_END, Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*cos(psi)*cos(theta)/m - g ==  0 );
//    ocp.subjectTo( AT_END, Wpendx == 0. );
//    ocp.subjectTo( AT_END, Wpendy == 0. );
//    ocp.subjectTo( AT_END, pendPitch == 0.0 );
//    ocp.subjectTo( AT_END, pendRoll == 0. );

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


#if AVOID_SINGULARITIES == 1
    //Constraint to avoid singularity
    ocp.subjectTo( -1.45 <= theta <= 1.45);
    ocp.subjectTo(-1.45 <= pendPitch <= 1.45);
#endif

    //Obstacle constraints
#if OBS==1
    ocp.subjectTo( 40000 <= (x*x+2*(z-3)*(z-3))*10000 );
    ocp.subjectTo( 40000 <= ((x-2)*(x-2)+2*(z-6)*(z-6))*10000 );
//    ocp.subjectTo( z<=5);
#endif

#if OBS==2
    ocp.subjectTo( 36 <= (36*(x-5)*(x-5)+(z+4)*(z+4)) );
    ocp.subjectTo( 36 <= (36*(x-5)*(x-5)+(z-10)*(z-10)) );
    ocp.subjectTo( 36 <= (36*(x-5)*(x-5)+(y+4)*(y+4)) );
    ocp.subjectTo( 36 <= (36*(x-5)*(x-5)+(y-10)*(y-10)) );

    ocp.subjectTo( 36 <= (36*((x+sin(pendRoll)*cos(pendPitch)*Xqg)-5)*((x+sin(pendRoll)*cos(pendPitch)*Xqg)-5)+((z+cos(pendPitch)*cos(pendRoll)*Xqg)+4)*((z+cos(pendPitch)*cos(pendRoll)*Xqg)+4)) );
    ocp.subjectTo( 36 <= (36*((x+sin(pendRoll)*cos(pendPitch)*Xqg)-5)*((x+sin(pendRoll)*cos(pendPitch)*Xqg)-5)+((z+cos(pendPitch)*cos(pendRoll)*Xqg)-10)*((z+cos(pendPitch)*cos(pendRoll)*Xqg)-10)) );
    ocp.subjectTo( 36 <= (36*((x+sin(pendRoll)*cos(pendPitch)*Xqg)-5)*((x+sin(pendRoll)*cos(pendPitch)*Xqg)-5)+((y-sin(pendPitch)*Xqg)+4)*((y-sin(pendPitch)*Xqg)+4)) );
    ocp.subjectTo( 36 <= (36*((x+sin(pendRoll)*cos(pendPitch)*Xqg)-5)*((x+sin(pendRoll)*cos(pendPitch)*Xqg)-5)+((y-sin(pendPitch)*Xqg)-10)*((y-sin(pendPitch)*Xqg)-10)) );

    ocp.subjectTo( 36 <= (36*((x-sin(pendRoll)*cos(pendPitch)*Xpg)-5)*((x-sin(pendRoll)*cos(pendPitch)*Xpg)-5)+((z-cos(pendPitch)*cos(pendRoll)*Xpg)+4)*((z-cos(pendPitch)*cos(pendRoll)*Xpg)+4)) );
    ocp.subjectTo( 36 <= (36*((x-sin(pendRoll)*cos(pendPitch)*Xpg)-5)*((x-sin(pendRoll)*cos(pendPitch)*Xpg)-5)+((z-cos(pendPitch)*cos(pendRoll)*Xpg)-10)*((z-cos(pendPitch)*cos(pendRoll)*Xpg)-10)) );
    ocp.subjectTo( 36 <= (36*((x-sin(pendRoll)*cos(pendPitch)*Xpg)-5)*((x-sin(pendRoll)*cos(pendPitch)*Xpg)-5)+((y+sin(pendPitch)*Xpg)+4)*((y+sin(pendPitch)*Xpg)+4)) );
    ocp.subjectTo( 36 <= (36*((x-sin(pendRoll)*cos(pendPitch)*Xpg)-5)*((x-sin(pendRoll)*cos(pendPitch)*Xpg)-5)+((y+sin(pendPitch)*Xpg)-10)*((y+sin(pendPitch)*Xpg)-10)) );

//    ocp.subjectTo( 1000 <= (1/9*(x-5)*(x-5)+1/36*(z+4)*(z+4))*1000 );
//    ocp.subjectTo( 1000 <= (1/9*(x-5)*(x-5)+1/36*(z-10)*(z-10))*1000 );
#endif

    VariablesGrid dist(1,timeGrid);
    dist.setZero();
    dist(23,0) = 4;
    dist(45,0) = -4;
//    dist.print();

    ocp.subjectTo(D == dist);
     // ----------------------------------------------------------------------------------- //






#if ALGO == 0 || ALGO==2

    // DEFINE A PLOT WINDOW:
    // ---------------------
    GnuplotWindow window2(PLOT_AT_EACH_ITERATION);
        window2.addSubplot( z,"DifferentialState z" );
        window2.addSubplot( x,"DifferentialState x" );
        window2.addSubplot( y,"DifferentialState y" );

        window2.addSubplot( pendPitch,"DifferentialState pendPitch" );
        window2.addSubplot( pendRoll,"DifferentialState pendRoll" );
                window2.addSubplot( mtot,"DifferentialState mtot" );
//        window2.addSubplot( vz,"DifferentialState vz" );
//        window2.addSubplot( Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*cos(phi)*cos(theta)/m - g,"acc z" );

        //window2.addSubplot( vx,"DifferentialState vx" );
//        window2.addSubplot( -Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*sin(theta)/m,"acc x" );
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

//        VariablesGrid p_init(2,timeGrid);
//        for (int i=0 ; i<nb_it+1; i++) {
//            p_init(i,0) = 3.;
//            p_init(i,1) = 6.;
//        }




    // ------------------------- INITIALIZE OCP PROBLEM ------------------------------- //
       OptimizationAlgorithm algorithm(ocp);
#if INIT==0 || INIT==1 || INIT==2 || INIT==3
       algorithm.initializeDifferentialStates(x_init);
#endif

       algorithm.initializeDisturbances(dist);
//       algorithm.getDisturbances(dist);
//       dist.print();
//       algorithm.initializeParameters(p_init);
//    algorithm.set(DISCRETIZATION_TYPE, SINGLE_SHOOTING);
//    algorithm.set( DYNAMIC_SENSITIVITY,  FORWARD_SENSITIVITY );

//        algorithm.set( HESSIAN_APPROXIMATION, EXACT_HESSIAN );
       algorithm.set(INTEGRATOR_TOLERANCE, 10e-1);
       algorithm.set(ABSOLUTE_TOLERANCE, 10e-1);
    algorithm.set( MAX_NUM_ITERATIONS, ocp_nb_it );
    algorithm.set( KKT_TOLERANCE, 1e-20 );
// 	algorithm.set( MAX_NUM_INTEGRATOR_STEPS, 4 );

#if PLOT == 0
    algorithm << window2;
#endif

        LogRecord logdiff(LOG_AT_END , PS_DEFAULT);
    logdiff << LOG_DIFFERENTIAL_STATES;
//logdiff << LOG_KKT_TOLERANCE;

#if PLOT == 1
        algorithm << logdiff;
#endif
        t = clock();

    algorithm.solve();

    t = clock() - t;
    std::cout << "total time : " << (((float)t)/CLOCKS_PER_SEC)<<std::endl;

    algorithm.getLogRecord(logdiff);

#if ALGO==0
    logdiff.print(file);
#endif

    VariablesGrid control(4,timeGrid);
    VariablesGrid state(16,timeGrid);
    DVector parameters(2);
    algorithm.getDifferentialStates(state);
    algorithm.getControls(control);
    algorithm.getParameters(parameters);
    parameters.print();

//    state.print();
    // -------------------- Export Results ------------------ //
//    VariablesGrid states, controls;
//    algorithm.getDifferentialStates(x_init);
//    x_init.print();

//    algorithm.getDifferentialStates(states);
//    algorithm.getControls("/tmp/controls.txt");
//    algorithm.getControls(controls);

#endif

#if ALGO == 1 || ALGO == 2
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

#if ALGO==2
    alg.initializeDifferentialStates(state);
    alg.initializeControls(control);
#endif

#if ALGO==1 && (INIT==0 || INIT==1 || INIT==2 || INIT==3 )
    alg.initializeDifferentialStates(x_init);
#endif

    // SET OPTION AND PLOTS WINDOW
    // ---------------------------
//    alg.set(DISCRETIZATION_TYPE, SINGLE_SHOOTING);
//        alg.set( HESSIAN_APPROXIMATION, BLOCK_BFGS_UPDATE );
//        alg.set( HESSIAN_APPROXIMATION, GAUSS_NEWTON_WITH_BLOCK_BFGS );
//            alg.set( HESSIAN_APPROXIMATION, GAUSS_NEWTON );
//    alg.set( HESSIAN_APPROXIMATION, FULL_BFGS_UPDATE );
//    alg.set( HESSIAN_APPROXIMATION, CONSTANT_HESSIAN );
//        alg.set( HESSIAN_APPROXIMATION, EXACT_HESSIAN );
        alg.set( GLOBALIZATION_STRATEGY, GS_LINESEARCH );

    GnuplotWindow window1(PLOT_AT_EACH_ITERATION);
    window1.addSubplot( x,"DifferentialState x" );
    window1.addSubplot( vx,"DifferentialState vx" );
    window1.addSubplot( -Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*sin(theta)/m,"acc x" );
    window1.addSubplot( z,"DifferentialState z" );
    window1.addSubplot( ((x-2)*(x-2)+2*(z-6)*(z-6)),"R2^2" );
    window1.addSubplot( (x*x+2*(z-3)*(z-3)),"R1^2" );
//    alg<<window1;


    // SETTING UP THE SIMULATION ENVIRONMENT,  RUN THE EXAMPLE...
    // ----------------------------------------------------------
    SimulationEnvironment sim( 0.0,10.,process,controller );

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

#if ALGO==3
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


