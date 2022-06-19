#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstring>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


#define M_HBARC 0.197

double alphaEM=1.0/137.0; double mllSqr=0.0; double qFSqrSum=1.0/9.0+4.0/9.0+1.0/9.0;

int QUARK_SUPPRESSION=1;

double Nc=3.0;
double Nf=3.0;
double nuG=2.0*(Nc*Nc-1.0);
double nuQ=4.0*Nc*Nf;

// ENERGY DENSITIES OF QUARKS AND GLUONS //
double eGlEq(double T){
    return (M_PI*M_PI/30.0)*nuG*T*T*T*T;
}

double eQuEq(double T){
    return (7.0*M_PI*M_PI/240.0)*nuQ*T*T*T*T;
}

double eEq(double T){
    return eGlEq(T)+eQuEq(T);
}

double eNEq(double T,double qSupp){
    return eGlEq(T)+qSupp*eQuEq(T);
}

// STANDARD RANDOM NUMBER GENERATOR //
double rng(){
    return drand48();
}

#include "PhaseSpaceDistribution.cpp"


// COMMANDLINE OPTIONS //
#include "IO/cfile.c"


    // GSL INTERPOLATION OBJECTS //
    gsl_interp_accel *EAcc,*QAcc;
    gsl_spline *EInt,*QInt;

    double wTMin; double wTMax;

    void Setup(){

        // SET DATA //
        double wTildeValues[161];
        double EValues[161];
        double QValues[161];

        std::ifstream EInStream;    EInStream.open("DATA/ECurveQCD.txt");
        std::ifstream QInStream;    QInStream.open("DATA/QCurveQCD.txt");


        int i=0;

        while(EInStream.good() && QInStream.good()){

            double wTE; double EVal;
            double wTQ; double QVal;

            EInStream >> wTE; EInStream >> EVal;
            QInStream >> wTQ; QInStream >> QVal;


            if(wTE==wTQ){
                wTildeValues[i]=wTE; EValues[i]=EVal; QValues[i]=QVal;
                i++;
            }
            else{
                std::cerr << "#ERROR -- COULD NOT PARSE HYDROATTRACTOR INPUT FILES CORRECTLY" << std::endl;
                exit(0);
            }

        }

        // SETUP SPLINE //
        EAcc = gsl_interp_accel_alloc ();
        EInt=gsl_spline_alloc(gsl_interp_cspline,161);
        gsl_spline_init(EInt,wTildeValues,EValues,161);

	cout << "Là ?" <<endl;
        QAcc = gsl_interp_accel_alloc ();
        QInt=gsl_spline_alloc(gsl_interp_cspline,161);
        gsl_spline_init(QInt,wTildeValues,QValues,161);


        wTMin=wTildeValues[0];
        wTMax=wTildeValues[160];

	cout << "ReLà ?" <<endl;
    }

    // ENERGY ATTRACTOR CURVE //
    double E(double wT){

        if(wT<wTMin){
            return 1.14480658493*std::pow(wT,4.0/9.0);
        }
        else if(wT>wTMax){
            return 1.0-2.0/(3.0*M_PI*wT);
        }
        else{
            return gsl_spline_eval(EInt,wT,EAcc);
        }

    }

    // DERIVATIVE OF ENERGY ATTRACTOR CURVE //
    double EPrime(double wT){

        if(wT<wTMin){
            return 1.14480658493*4.0/9.0*std::pow(wT,-5.0/9.0);
        }
        else if(wT>wTMax){
            return 2.0/(3.0*M_PI*wT*wT);
        }
        else{
            return gsl_spline_eval_deriv(EInt,wT,EAcc);
        }
    }

    // pL/E ATTRACTOR CURVE //
    double P(double wT){

        if(wT<wTMin){
            return (wT/wTMin)*P(wTMin);
        }
        else{
            return (-4.0*E(wT)+9.0*wT*EPrime(wT))/(-12.0*E(wT)+3.0*wT*EPrime(wT));
        }
    }

    // eQ/eG ATTRACTOR CURVE //
    double Q(double wT){

        if(wT<wTMin){
            return 0.0;
        }
        else if(wT>wTMax){
            return eQuEq(1.0)/eGlEq(1.0);
        }
        else{
            return gsl_spline_eval(QInt,wT,QAcc);
        }

    }

    // GET VALUES OF T,wT,e,pL,eQ/eG //
    void GetValues(double dNchdEta,double Area,double etaOverS,double Tau,double &T,double &wTilde,double &e,double &pL,double &eQOvereG){

        // DETERMINE (tau T^3)_hydro //
        double SByNCh=6.7; double nuEff=32.0;
        double tauTCubeHydro=8.26*(dNchdEta/1900.0)/(Area/110.0)*(SByNCh/6.7)/(nuEff/32.0); // fm^{-2}

        // DETERMINE (e(tau) tau^{4/3})_{infty} //
        double eTau43Infty=eEq(1.0)*std::pow(tauTCubeHydro,4.0/3.0);

        //////////////////////////////////////////////////////////
        // DETERMINE TEMPERATURE SELF-CONSISTENTLY ACCORDING TO //
        // e(T)tau^{4/3} = E(wTilde) (e(tau) tau^{4/3})_{infty} //
        //          wTilde= (T tau)/(4pi eta/s)                 //
        //////////////////////////////////////////////////////////

        double TLow=0.0; double THigh=std::pow(eTau43Infty/(eEq(1.0)*std::pow(Tau,4.0/3.0)),1.0/4.0);

        double TMid=(THigh+TLow)/2.0;
        double wTildeMid=(TMid*Tau)/(4.0*M_PI*etaOverS);


        while(THigh-TLow>1E-4*TMid){

            if(E(wTildeMid)/std::pow(TMid,4)>eEq(1.0)*std::pow(Tau,4.0/3.0)/eTau43Infty){
                TLow=TMid;
            }
            else{
                THigh=TMid;
            }

            TMid=(THigh+TLow)/2.0;
            wTildeMid=(TMid*Tau)/(4.0*M_PI*etaOverS);

        }

        wTilde=(TMid*Tau)/(4.0*M_PI*etaOverS);

        T=TMid*M_HBARC;
        e=eEq(T);
        pL=e*P(wTilde);
        eQOvereG=Q(wTilde);

//	cout << "T=" << T << ", e=" << e << ", pL=" << pL << endl;


    }

namespace DileptonRates{

    // NON-EQUILIBIRIUM CURRENT-CURRENT CORRELATION FUNCTION -- yQ IS PSEUDORPAIDITY OF DILEPTON PAIR //
    double SampleTracePi(double q0,double qT,double PhiQ, double yQ,double EtaX,double Xi,double Teff,double qSupp){

        // SET JACOBIAN TO UNITY //
        double Jacobian=1.0;

        // SET KINETMATIC VARIABLES DERIVED FROM q //
        double qZ=qT*sinh(yQ);
        double qAbs=qT*cosh(yQ);
//        double CosThetaQ=qZ/qAbs;
        double CosThetaQ=tanh(yQ);
        double SinThetaQ=sqrt(1.0-CosThetaQ*CosThetaQ);

        // GET COORDINATE SYSTEM FOR p INTEGRATION //
        double eq[3];

        eq[0]=cos(PhiQ)/sqrt(1.0+sinh(yQ)*sinh(yQ));
        eq[1]=sin(PhiQ)/sqrt(1.0+sinh(yQ)*sinh(yQ));
        eq[2]=sinh(yQ)/sqrt(1.0+sinh(yQ)*sinh(yQ));


        double el[3];

        el[0]=(0.0-CosThetaQ*eq[0])/SinThetaQ;
        el[1]=(0.0-CosThetaQ*eq[1])/SinThetaQ;
        el[2]=(1.0-CosThetaQ*eq[2])/SinThetaQ;


        double et[3];

        et[0]=(eq[1]*el[2]-eq[2]*el[1]);
        et[1]=(eq[2]*el[0]-eq[0]*el[2]);
        et[2]=(eq[0]*el[1]-eq[1]*el[0]);

        // SAMPLE p VECTOR //
        double pMin=(q0-qAbs)/2.0;
        double pMax=(q0+qAbs)/2.0;

        double pAbs=pMin+(pMax-pMin)*rng();
        double PhiP=2.0*M_PI*drand48();

        Jacobian*=2.0*M_PI*(pMax-pMin);

        double CosThetaPQ=(qAbs*qAbs+2.0*pAbs*q0-q0*q0)/(2.0*pAbs*qAbs);
        double SinThetaPQ=sqrt(1.0-CosThetaPQ*CosThetaPQ);

        // GET RELEVANT VECTORS //
        double qVec[3]; double pVec[3];  double qMpVec[3];


        qVec[0]=qAbs*eq[0];
        qVec[1]=qAbs*eq[1];
        qVec[2]=qAbs*eq[2];

        pVec[0]=pAbs*(SinThetaPQ*cos(PhiP)*et[0]+SinThetaPQ*sin(PhiP)*el[0]+CosThetaPQ*eq[0]);
        pVec[1]=pAbs*(SinThetaPQ*cos(PhiP)*et[1]+SinThetaPQ*sin(PhiP)*el[1]+CosThetaPQ*eq[1]);
        pVec[2]=pAbs*(SinThetaPQ*cos(PhiP)*et[2]+SinThetaPQ*sin(PhiP)*el[2]+CosThetaPQ*eq[2]);

        qMpVec[0]=qVec[0]-pVec[0];
        qMpVec[1]=qVec[1]-pVec[1];
        qMpVec[2]=qVec[2]-pVec[2];

        double pT=sqrt(pVec[0]*pVec[0]+pVec[1]*pVec[1]);
        double yP=atanh(pVec[2]/pAbs);

        double qMpT=sqrt(qMpVec[0]*qMpVec[0]+qMpVec[1]*qMpVec[1]);
        double qMpAbs=sqrt(qMpVec[0]*qMpVec[0]+qMpVec[1]*qMpVec[1]+qMpVec[2]*qMpVec[2]);
        double yqMP=atanh(qMpVec[2]/qMpAbs);

        // GET PHASE SPACE DISTRIBUTION //
        double fp=PhaseSpaceDistribution::fQ(pT,yP,EtaX,Xi,Teff,qSupp);
        double fqMp=PhaseSpaceDistribution::fQ(qMpT,yqMP,EtaX,Xi,Teff,qSupp);

        // GET POLARIZATION TENSOR //
        return Nc*(q0*q0-qAbs*qAbs)/(4.0*M_PI*M_PI*qAbs)*Jacobian*fp*fqMp;

    }


    // SAMPLE DILPETON PRODUCTION -- QSqr INVARIANT MASS SQUARE, qT TRANSVERSE MOMENTUM , EtaQ RAPIDITY OF DILEPTON PAIR //
    void SampledNdQdy(double QSqr,double qTMin,double qTMax,double TauMin,double TauMax,double EtaQ,double dNchdEta,double Area,double etas,double &dN,double &dNPreEq,double &dNHydro){

        // SAMPLE INTEGRATION POINT //
        double Jacobian=1.0;

        double EtaMin=-8+EtaQ;
        double EtaMax=8+EtaQ;
        double EtaX=EtaMin+(EtaMax-EtaMin)*rng();
        double qT=qTMin+(qTMax-qTMin)*rng();
        double PhiQ=2.0*M_PI*rng();

        double Q=sqrt(QSqr);

        // ENERGY AND MOMENTUM OF DILEPTON PAIR //
        double qZ=std::sqrt(QSqr+qT*qT)*sinh(EtaQ);
        double qAbs=std::sqrt(qT*qT+(QSqr+qT*qT)*sinh(EtaQ)*sinh(EtaQ));
        double q0=std::sqrt(QSqr+qAbs*qAbs);


        // PSEUDO-RAPIDITY OF DILEPTON PAIR //
        double yQ=atanh(qZ/qAbs);

        // JACOBIAN  -- d^4Q=QdQ dy d^2qT //
        Jacobian*=2.0*M_PI*(qTMax-qTMin)*(EtaMax-EtaMin)*qT*Q;

        // EVOLUTION TIME //
        double Tau=TauMin+(TauMax-TauMin)*rng();

        Jacobian*=Tau*(TauMax-TauMin)*Area/(M_HBARC*M_HBARC*M_HBARC*M_HBARC);

        // GET EVOLUTION OF MACROSCOPIC FIELDS //
        double T,wTilde,e,pL,eQOvereG;
        GetValues(dNchdEta,Area,etas,Tau,T,wTilde,e,pL,eQOvereG);

        // GET PARAMETERS OF PHASE-SPACE DISTRIBUTION //
        double Xi,Teff,qSupp;
        PhaseSpaceDistribution::GetPhaseSpaceDistributionParameters(e,pL,eQOvereG,Xi,Teff,qSupp);

        // SAMPLE DILEPTON PRODUCTION //
        double PreFactor=alphaEM*alphaEM/(6.0*M_PI*M_PI*M_PI*QSqr)*(1.0+mllSqr/QSqr)*sqrt(1.0-4.0*mllSqr/QSqr)*qFSqrSum;
        double dNlld4xd4Q=PreFactor*SampleTracePi(q0,qT,PhiQ,yQ,EtaX,Xi,Teff,qSupp);

        // GET PRODUCTION YIELD //
        dN=Jacobian*dNlld4xd4Q;

        // SEPARATE INTO PRE-EQ AND HYDRO //
        if(wTilde<1.0){
            dNPreEq=dN; dNHydro=0.0;
        }
        else{
            dNPreEq=0.0; dNHydro=dN;
        }


    }
    // SAMPLE DILPETON PRODUCTION -- QSqr INVARIANT MASS SQUARE, qT TRANSVERSE MOMENTUM , EtaQ RAPIDITY OF DILEPTON PAIR //
    void SampledNd2qTdQdy(double QMin,double QMax,int NQ,double qTMin,double qTMax,int NqT,double TauMin,double TauMax,double EtaQ,double dNchdEta,double Area,double *dN,double *dNPreEq,double *dNHydro,double eta_over_s){

        // GET EVOLUTION OF MACROSCOPIC FIELDS //

        for(int ib=0;ib<NQ*NqT;ib++){
            dN[ib]=0.0; dNHydro[ib]=0.0;dNPreEq[ib]=0.0;
        }

        // SET JACOBIAN FOR MONTE-CARLO INTEGRATION //
        double Jacobian=1.0;

        // SAMPLE EVOLUTION TIME //
        double Tau=TauMin+(TauMax-TauMin)*rng();
        Jacobian*=Tau*(TauMax-TauMin)*Area/(M_HBARC*M_HBARC*M_HBARC*M_HBARC);

        double T,wTilde,e,pL,eQOvereG;
        GetValues(dNchdEta,Area,eta_over_s,Tau,T,wTilde,e,pL,eQOvereG);

            // CALCULATE DILEPTON PRODUCTION dN/dydQd2qT FOR ALL BINS //
            for(int iQ=0;iQ<NQ;iQ++){

                double Q=QMin+(QMax-QMin)*iQ/(NQ-1);

                for(int iqT=0;iqT<NqT;iqT++){

                    double qT=qTMin+(qTMax-qTMin)*iqT/(NqT-1);

                    // SAMPLE INTEGRATION POINT //
                    double EtaMin=-8+EtaQ;
                    double EtaMax=8+EtaQ;
                    double EtaX=EtaMin+(EtaMax-EtaMin)*rng();
                    double PhiQ=2.0*M_PI*rng();

                    double QSqr=Q*Q;

                    // ENERGY AND MOMENTUM OF DILEPTON PAIR //
                    double qZ=std::sqrt(QSqr+qT*qT)*sinh(EtaQ);
                    double qAbs=std::sqrt(qT*qT+(QSqr+qT*qT)*sinh(EtaQ)*sinh(EtaQ));
                    double q0=std::sqrt(QSqr+qAbs*qAbs);

                    // PSEUDO-RAPIDITY OF DILEPTON PAIR //
                    double yQ=atanh(qZ/qAbs);


                    // JACOBIAN  -- d^4Q=QdQ dy d^2qT //
                    Jacobian*=(EtaMax-EtaMin)*Q;


                    // CALCULATE DILEPTON PRODUCTION //
                    double Xi,Teff,qSupp;
                    PhaseSpaceDistribution::GetPhaseSpaceDistributionParameters(e,pL,eQOvereG,Xi,Teff,qSupp);
                    double PreFactor=alphaEM*alphaEM/(6.0*M_PI*M_PI*M_PI*QSqr)*(1.0+2.0*mllSqr/QSqr)*sqrt(1.0-4.0*mllSqr/QSqr)*qFSqrSum;
                    double dNlld4xd4Q=PreFactor*DileptonRates::SampleTracePi(q0,qT,PhiQ,yQ,EtaX,Xi,Teff,qSupp);

                    double dNSamp=Jacobian*dNlld4xd4Q;

		if(dNSamp>1){
		   std::cout << "Q = " << Q << ", qT = " << qT << ", Rate = " << dNSamp << std::endl;
		}

                    // DISTINGUSIH BETWEEN PRE-EQUILIBRIUM AND HYDRO //
                    dN[iqT+NqT*iQ]+=dNSamp;

                    if(wTilde<1.0){
                        dNPreEq[iqT+NqT*iQ]=dNSamp;
                        dNHydro[iqT+NqT*iQ]=0;
                    }
                    else{
                        dNHydro[iqT+NqT*iQ]=dNSamp;
                        dNPreEq[iqT+NqT*iQ]=0;
                    }

                }
	    }
	}

    void SampledNd2qTdQdy_v2(double Q,double qT,double TauMin,double TauMax,double EtaQ,double dNchdEta,double Area,double &dN,double &dNPreEq,double &dNHydro,double eta_over_s){


        // SET JACOBIAN FOR MONTE-CARLO INTEGRATION //
        double Jacobian=1.0;

        // SAMPLE EVOLUTION TIME //
        double Tau=TauMin+(TauMax-TauMin)*rng();
        Jacobian*=Tau*(TauMax-TauMin)*Area/(M_HBARC*M_HBARC*M_HBARC*M_HBARC);

        double T,wTilde,e,pL,eQOvereG;
        GetValues(dNchdEta,Area,eta_over_s,Tau,T,wTilde,e,pL,eQOvereG);

            // CALCULATE DILEPTON PRODUCTION dN/dydQd2qT FOR ALL BINS //

                    // SAMPLE INTEGRATION POINT //
                    double EtaMin=-8+EtaQ;
                    double EtaMax=8+EtaQ;
                    double EtaX=EtaMin+(EtaMax-EtaMin)*rng();
                    double PhiQ=2.0*M_PI*rng();

                    double QSqr=Q*Q;

                    // ENERGY AND MOMENTUM OF DILEPTON PAIR //
                    double qZ=std::sqrt(QSqr+qT*qT)*sinh(EtaQ);
                    double qAbs=std::sqrt(qT*qT+(QSqr+qT*qT)*sinh(EtaQ)*sinh(EtaQ));
                    double q0=std::sqrt(QSqr+qAbs*qAbs);

                    // PSEUDO-RAPIDITY OF DILEPTON PAIR //
                    double yQ=atanh(qZ/qAbs);


                    // JACOBIAN  -- d^4Q=QdQ dy d^2qT //
                    Jacobian*=(EtaMax-EtaMin)*Q*qT*2*M_PI;


                    // CALCULATE DILEPTON PRODUCTION //
                    double Xi,Teff,qSupp;
                    PhaseSpaceDistribution::GetPhaseSpaceDistributionParameters(e,pL,eQOvereG,Xi,Teff,qSupp);
                    double PreFactor=alphaEM*alphaEM/(6.0*M_PI*M_PI*M_PI*QSqr)*(1.0+2.0*mllSqr/QSqr)*sqrt(1.0-4.0*mllSqr/QSqr)*qFSqrSum;
                    double dNlld4xd4Q=PreFactor*DileptonRates::SampleTracePi(q0,qT,PhiQ,yQ,EtaX,Xi,Teff,qSupp);

                    double dNSamp=Jacobian*dNlld4xd4Q;

//		   std::cout << "Q = " << Q << ", qT = " << qT << ", Jacob = " << Jacobian << ", dNd4x = "<< dNlld4xd4Q << std::endl;

                    // DISTINGUSIH BETWEEN PRE-EQUILIBRIUM AND HYDRO //
                    dN=dNSamp;

                    if(wTilde<1.0){
                        dNPreEq=dNSamp;
                        dNHydro=0;
                    }
                    else{
                        dNHydro=dNSamp;
                        dNPreEq=0;
                    }

    }

    void SampledNd2qTdQdeta(double Q,double qT,double TauMin,double TauMax,double EtaQ,double dNchdEta,double Area,double &dN,double &dNPreEq,double &dNHydro,double eta_over_s){


        // SET JACOBIAN FOR MONTE-CARLO INTEGRATION //
        double Jacobian=1.0;

        // SAMPLE EVOLUTION TIME //
        double Tau=TauMin+(TauMax-TauMin)*rng();
        Jacobian*=Tau*(TauMax-TauMin)*Area/(M_HBARC*M_HBARC*M_HBARC*M_HBARC);

        double T,wTilde,e,pL,eQOvereG;
        GetValues(dNchdEta,Area,eta_over_s,Tau,T,wTilde,e,pL,eQOvereG);

            // CALCULATE DILEPTON PRODUCTION dN/dydQd2qT FOR ALL BINS //

                    // SAMPLE INTEGRATION POINT //
                    double EtaMin=-8+EtaQ;
                    double EtaMax=8+EtaQ;
                    double EtaX=EtaMin+(EtaMax-EtaMin)*rng();
                    double PhiQ=2.0*M_PI*rng();

                    double QSqr=Q*Q;

                    // ENERGY AND MOMENTUM OF DILEPTON PAIR //
                    double qZ=std::sqrt(QSqr+qT*qT)*sinh(EtaQ);
                    double qAbs=std::sqrt(qT*qT+(QSqr+qT*qT)*sinh(EtaQ)*sinh(EtaQ));
                    double q0=std::sqrt(QSqr+qAbs*qAbs);

                    // PSEUDO-RAPIDITY OF DILEPTON PAIR //
                    double yQ=atanh(qZ/qAbs);


                    // JACOBIAN  -- d^4Q=QdQ dy d^2qT //
                    Jacobian*=(EtaMax-EtaMin)*Q*qT*2*M_PI;


                    // CALCULATE DILEPTON PRODUCTION //
                    double Xi,Teff,qSupp;
                    PhaseSpaceDistribution::GetPhaseSpaceDistributionParameters(e,pL,eQOvereG,Xi,Teff,qSupp);
                    double PreFactor=alphaEM*alphaEM/(6.0*M_PI*M_PI*M_PI*QSqr)*(1.0+2.0*mllSqr/QSqr)*sqrt(1.0-4.0*mllSqr/QSqr)*qFSqrSum;
                    double dNlld4xd4Q=PreFactor*DileptonRates::SampleTracePi(q0,qT,PhiQ,yQ,EtaX,Xi,Teff,qSupp);

                    double dNSamp=Jacobian*dNlld4xd4Q*(qAbs/q0);

    //		   std::cout << "Q = " << Q << ", qT = " << qT << ", Jacob = " << Jacobian << ", dNd4x = "<< dNlld4xd4Q << std::endl;

                    // DISTINGUSIH BETWEEN PRE-EQUILIBRIUM AND HYDRO //
                    dN=dNSamp;

                    if(wTilde<1.0){
                        dNPreEq=dNSamp;
                        dNHydro=0;
                    }
                    else{
                        dNHydro=dNSamp;
                        dNPreEq=0;
                    }

                }



}

        // PSEUDO-RAPIDITY OF DILEPTON PAIR //
int main(int argc, char **argv) {

    std::ios_base::sync_with_stdio(false); // Faster output

    // SET COMMANDLINE ARGUMENTS //
    Konfig CommandlineArguments(argc,argv);

    // COLLISION PARAMETERS //
    double EtaOverS=0.32; double dNchdEta=1900; double Area=110;

    CommandlineArguments.Getval("etas",EtaOverS);
    CommandlineArguments.Getval("Nch",dNchdEta);
    CommandlineArguments.Getval("area",Area);
    CommandlineArguments.Getval("Q",QUARK_SUPPRESSION);

    std::cerr << "#CALCULATING DILEPTON PRODUCTION FOR dNchdEta=" << dNchdEta << " Area=" << Area << " fm^2 AND Eta/s=" << EtaOverS << " QUARK SUPPRESION " << QUARK_SUPPRESSION << std::endl;


    // DILEPTON PARAMTERS //
    double QMin=1.0; double QMax=5.0;

    CommandlineArguments.Getval("QMin",QMin);
    CommandlineArguments.Getval("QMax",QMax);

    double qTMin=1; double qTMax=5.0; double TauMin=0.0; double TauMax=40.0;

    CommandlineArguments.Getval("qTMin",qTMin);
    CommandlineArguments.Getval("qTMax",qTMax);
    CommandlineArguments.Getval("TauMin",TauMin);
    CommandlineArguments.Getval("TauMax",TauMax);

    std::cerr << "#KINEMATIC CUTS ARE qT=" << qTMin << " - " << qTMax << " AND  tau=" << TauMin << "-" << TauMax << " fm" << std::endl;


    // MONTE CARLO SAMPLING //
    int NSamples=512000;
    CommandlineArguments.Getval("NSamples",NSamples);

    // SEED RANDOM NUMBER GENERATOR //
    srand48(time(0));

    // SETUP INTERPOLATORS FOR HYDRO ATTRACTORS //
        double wTildeValues[161];
        double EValues[161];
        double QValues[161];

        std::ifstream EInStream;    EInStream.open("DATA/ECurveQCD.txt");
        std::ifstream QInStream;    QInStream.open("DATA/QCurveQCD.txt");


        int i=0;

        while(EInStream.good() && QInStream.good()){

            double wTE; double EVal;
            double wTQ; double QVal;

            EInStream >> wTE; EInStream >> EVal;
            QInStream >> wTQ; QInStream >> QVal;


            if(wTE==wTQ){
                wTildeValues[i]=wTE; EValues[i]=EVal; QValues[i]=QVal;
                i++;
            }
            else{
                std::cerr << "#ERROR -- COULD NOT PARSE HYDROATTRACTOR INPUT FILES CORRECTLY" << std::endl;
                exit(0);
            }

        }

        // SETUP SPLINE //
        EAcc = gsl_interp_accel_alloc ();
        EInt=gsl_spline_alloc(gsl_interp_cspline,161);
        gsl_spline_init(EInt,wTildeValues,EValues,161);

        QAcc = gsl_interp_accel_alloc ();
        QInt=gsl_spline_alloc(gsl_interp_cspline,161);
        gsl_spline_init(QInt,wTildeValues,QValues,161);


        wTMin=wTildeValues[0];
        wTMax=wTildeValues[160];

    // CALCULATE DILEPTON PRODUCTION -- dN/dQdyQ //
    int NQ=10; int NqT=10; double yQ=2.0;

    CommandlineArguments.Getval("NQ",NQ);
    CommandlineArguments.Getval("NqT",NqT);
    CommandlineArguments.Getval("yQ",yQ);

    std::cerr << "#CALCULATING FOR Q=" << QMin << " - " << QMax << " IN " << NQ << " BINS AT yQ=" <<  yQ  << " WITH " << NSamples << " SAMPLES PER BIN" << std::endl;

    // WRITE HEADER //
    std::cout << "#1-Q [GeV] 2-qT[GeV] 3--dN/dqTdQdY [GeV-2] 4--dN_{PreEq}/dqTdQdY [GeV-2] 5--dN_{Hydro}/dqTdQdY [GeV-2]" << std::endl;


        double dNlldQdY[NQ*NqT];
        double dNlldQdYPreEq[NQ*NqT];
        double dNlldQdYHydro[NQ*NqT];

//        double dN[NQ*NqT],dNPreEq[NQ*NqT],dNHydro[NQ*NqT];

        for(int ib=0;ib<NQ*NqT;ib++){
            dNlldQdY[ib]=0.0; dNlldQdYPreEq[ib]=0.0;dNlldQdYHydro[ib]=0.0;
        }

//        for(int i=0;i<NSamples;i++){
//	        for(int ib=0;ib<NQ*NqT;ib++){
//        	    dN[ib]=0.0; dNPreEq[ib]=0.0;dNHydro[ib]=0.0;
//        	}
//            DileptonRates::SampledNd2qTdQdy(QMin,QMax,NQ,qTMin,qTMax,NqT,TauMin,TauMax,yQ,dNchdEta,Area,dN,dNPreEq,dNHydro,EtaOverS);
//	    for(int ib=0;ib<NQ*NqT;ib++){
//            	dNlldQdY[ib]+=dN[ib]; dNlldQdYPreEq[ib]+=dNPreEq[ib]; dNlldQdYHydro[ib]+=dNHydro[ib];
//	    }
//        }
        for(int i=0;i<NSamples;i++){
            std::cerr << i << "/" << NSamples << '\r';
            for(int iQ=0;iQ<NQ;iQ++){
                double Q=QMin+(QMax-QMin)*iQ/(NQ-1);
                for(int iqT=0;iqT<NqT;iqT++){
                    double qT=qTMin+(qTMax-qTMin)*iqT/(NqT-1);
        	    double dN=0.0; double dNPreEq=0.0; double dNHydro=0.0;
            	    DileptonRates::SampledNd2qTdQdeta(Q,qT,TauMin,TauMax,yQ,dNchdEta,Area,dN,dNPreEq,dNHydro,EtaOverS);
		    dNlldQdY[iqT+NqT*iQ]+=dN;
		    dNlldQdYPreEq[iqT+NqT*iQ]+=dNPreEq;
		    dNlldQdYHydro[iqT+NqT*iQ]+=dNHydro;
		}
	    }
	}

    std::cerr << endl;

	for(int ib=0;ib<NQ*NqT;ib++){
	        dNlldQdY[ib]/=double(NSamples);
        	dNlldQdYPreEq[ib]/=double(NSamples);
	        dNlldQdYHydro[ib]/=double(NSamples);
	}


	double sum[NQ];
	for(int iQ=0;iQ<NQ;iQ++){
		double Q=QMin+(QMax-QMin)*iQ/(NQ-1);
		sum[iQ]=0;

		for(int iqT=0;iqT<NqT;iqT++){
                    double qT=qTMin+(qTMax-qTMin)*iqT/(NqT-1);
      	            std::cout << Q << " " << qT << " " << dNlldQdY[iqT+NqT*iQ] << " " << dNlldQdYPreEq[iqT+NqT*iQ] << " " << dNlldQdYHydro[iqT+NqT*iQ] << std::endl;
		    sum[iQ]+=dNlldQdY[iqT+NqT*iQ]*(qTMax-qTMin)/NqT;
		}
	}
	for(int iQ=0;iQ<NQ;iQ++){
		double Q=QMin+(QMax-QMin)*iQ/(NQ-1);
      	            std::cout << "\n" << Q << " " << sum[iQ] << std::endl;

	}

//	double N=0.0; double Npre=0.0; double Nhy=0.0;
//	DileptonRates::SampledNdQdy(1*1,1,1,10,10,yQ,dNchdEta,Area,EtaOverS,N,Npre,Nhy);
  //    	            std::cout << "Old calculation : " << N << " " << Npre << " " << Nhy << std::endl;


//    // CALCULATE DILEPTON PRODUCTION -- dN/dy //
//    double yMin=-2.0; double yMax=+2.0; int NY=15;
//
//    CommandlineArguments.Getval("yMin",QMin);
//    CommandlineArguments.Getval("yMax",QMax);
//    CommandlineArguments.Getval("NY",NY);
//
//    std::cerr << "#CALCULATING FOR y=" << yMin << " - " << yMax << " IN " << NY << " BINS WITH Q=" <<  QMin  << " - " << QMax << " WITH " << NSamples << " SAMPLES PER BIN" << std::endl;
//
//    for(int iY=0;iY<NY;iY++){
//
//        double yQ=yMin+iY*(yMax-yMin)/double(NY-1);
//
//        double dNlldY=0.0;
//
//        for(int i=0;i<NSamples;i++){
//            dNlldY+=DileptonRates::SampledNdy(QMin,QMax,qTMin,qTMax,TauMin,TauMax,yQ,dNchdEta,Area,EtaOverS);
//        }
//        dNlldY/=double(NSamples);
//
//        std::cout << yQ << " " << dNlldY << std::endl;
//
//    }

    // EXIT //
    exit(1);


}
