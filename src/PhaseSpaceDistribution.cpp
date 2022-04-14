namespace PhaseSpaceDistribution {
    
    // BOOST INVARIANT PARAMETRIZATION OF THE PHASE-SPACE DISTRIBUTION OF GLUONS //
    double fG(double pT,double yP,double EtaX,double Xi,double Teff,double qSupp){
        
        double Arg=sqrt(1.0+Xi*Xi*sinh(yP-EtaX)*sinh(yP-EtaX))*pT/Teff;
        
        return 1.0/(exp(Arg)-1.0);
    }
    
    
    // BOOST INVARIANT PARAMETRIZATION OF THE PHASE-SPACE DISTRIBUTION OF QUARKS //
    double fQ(double pT,double yP,double EtaX,double Xi,double Teff,double qSupp){
        
        double Arg=sqrt(1.0+Xi*Xi*sinh(yP-EtaX)*sinh(yP-EtaX))*pT/Teff;
        
        return qSupp*1.0/(exp(Arg)+1.0);
    }
    
    double fQBolt(double pT,double yP,double EtaX,double Xi,double Teff,double qSupp){
        
        double Arg=sqrt(1.0+Xi*Xi*sinh(yP-EtaX)*sinh(yP-EtaX))*pT/Teff;
        
        return qSupp*1.0/(exp(Arg));
    }
    
    
    // ANISOTROPY FUNCTIONS C(Xi) AND S(Xi) //
    double CXi(double xi){
        
        if(xi<1){
            std::cerr << "#ERROR -- BEHAVIOR OF C(Xi) NOT DEFINED FOR Xi=" << xi << std::endl;
            exit(0);
        }
        else if(xi<1.000001){
            return 1.0-4.0/3.0*(xi-1.0);
        }
        else{
            return 0.5*(1.0/(xi*xi)+atan(sqrt(xi*xi-1.0))/sqrt(xi*xi-1.0));
        }
    }
    
    double SXi(double xi){
        
        if(xi<1){
            std::cerr << "#ERROR -- BEHAVIOR OF S(Xi) NOT DEFINED FOR Xi=" << xi << std::endl;
            exit(0);
        }
        else if(xi<1.000001){
            return 1.0-4.0/5.0*(xi-1.0);
        }
        else{
            return 0.5*(1.0/(xi*xi-xi*xi*xi*xi)+atan(sqrt(xi*xi-1.0))/std::pow(xi*xi-1.0,3.0/2.0));
        }
    }
    
    
    // LONG. PRESSURE TO ENERGY DENSITY RATIO //
    double pLOverE(double xi){
        return SXi(xi)/CXi(xi);
    }
    
    // GET PARAMETERS OF THE PHASE-SPACE DISTRIBUTION //
    void GetPhaseSpaceDistributionParameters(double e,double pL,double eQOvereG,double &Xi,double &TEff,double &qSupp){
        
        // DETERMINE ANSIOTROPY PARAMETER Xi //
        double XiHigh=sqrt(e/pL); double XiLow=1.0; double XiMid=(XiHigh+XiLow)/2.0;
        
        while(XiHigh-XiLow>1e-6*XiMid){
            
            if(pL/e>pLOverE(XiMid)){
                XiHigh=XiMid;
            }
            else{
                XiLow=XiMid;
            }
            
            XiMid=(XiHigh+XiLow)/2.0;
            
        }
        
        Xi=XiMid;
        
        // DETERMINE QUARK SUPPRESSION FACTOR //
        if(QUARK_SUPPRESSION==1){
            qSupp=eQOvereG*(eGlEq(1.0)/eQuEq(1.0));
        }
        else{
            qSupp=1.0;
        }
        
        // DETERMINE EFFECTIVE TEMPERATURE TEff //
        TEff=std::pow(e/(eNEq(1.0,qSupp)*CXi(Xi)),1.0/4.0);
        
        
    }
    
}
