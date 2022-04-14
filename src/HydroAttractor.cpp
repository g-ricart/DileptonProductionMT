namespace HydroAttractor {
    
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
        
        QAcc = gsl_interp_accel_alloc ();
        QInt=gsl_spline_alloc(gsl_interp_cspline,161);
        gsl_spline_init(QInt,wTildeValues,QValues,161);
        
        
        wTMin=wTildeValues[0];
        wTMax=wTildeValues[160];
        
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
        
        
    }
    
}
