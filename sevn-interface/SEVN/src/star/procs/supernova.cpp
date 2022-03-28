//
// Created by spera on 30/12/18.
//



#include <supernova.h>
#include <star.h>
#include <remnant.h>

supernova::supernova(Star *s){

    fallback_frac   = -1.0; //it means "not-set, yet"
    Average_remnant = -1.0; //it means "not-set, yet"
    Average_ejected = -1.0; //it means "not-set, yet"
    remnant_type    = Lookup::Remnants::NotARemnant; //not a remnant

    if(s != nullptr) { //This means that I am calling the SN constructor from the specific instance of the SN type, and not just for registring

        std::string kicksmodel = s->get_svpar_str("sn_kicks");
        kick = Kicks::Instance(kicksmodel);

        if (kick == nullptr)
            svlog.critical("Unknown SN kick model: [" + kicksmodel + "]", __FILE__, __LINE__);
    }


}

void supernova::initialise_remnant(Star *s, double Mass_remnant, Lookup::Remnants Remnant_type){


    ///Check input
    if (Mass_remnant>0 and Remnant_type==Lookup::Remnants::Empty)
        svlog.critical("It is not possible to use a Mass (" + utilities::n2s(Mass_remnant,__FILE__,__LINE__) +
        ") larger than 0 with Remnant_type=Empty",__FILE__,__LINE__,sevnstd::sevnio_error(""));
    else if (Mass_remnant<=0 and Remnant_type!=Lookup::Remnants::Empty)
        svlog.critical("It is not possible to use a Mass<=0 (" + utilities::n2s(Mass_remnant,__FILE__,__LINE__) +
        ") if the Remnant_type is not Empty",__FILE__,__LINE__,sevnstd::sevnio_error(""));
    else if (Mass_remnant>s->get_svpar_num("sn_Mchandra") and (Remnant_type==Lookup::Remnants::HeWD or Remnant_type==Lookup::Remnants::COWD or Remnant_type==Lookup::Remnants::ONeWD))
        svlog.critical("It is not possible to use a Mass (" + utilities::n2s(Mass_remnant,__FILE__,__LINE__) +
        ") larger than sn_Mchandra with Remnant_type=WD",__FILE__,__LINE__,sevnstd::sevnio_error(""));
    else if (Mass_remnant>s->get_svpar_num("sn_max_ns_mass") and (Remnant_type==Lookup::Remnants::NS_CCSN or Remnant_type==Lookup::Remnants::NS_ECSN) )
        svlog.critical("It is not possible to use a Mass (" + utilities::n2s(Mass_remnant,__FILE__,__LINE__) +
        ") larger than sn_max_ns_mass with Remnant_type=NS",__FILE__,__LINE__,sevnstd::sevnio_error(""));


    ///
    if (Mass_remnant>0.0 ) {
        Mremnant = Mass_remnant;
        remnant_type =  Remnant_type;
        set_staremnant(s);
        s->set_remnant();
    }
    else {
        Mremnant = 0.0;
        remnant_type = Lookup::Remnants::Empty;
        s->set_empty();
    }

    Mejected = 0;


    return;

}

supernova* supernova::Instance(std::string const &name, Star *s){


    auto it = GetStaticMap().find(name);
    if(it!= GetStaticMap().end()) GetUsed()[std::distance(GetStaticMap().begin(), it)] = 1;
    return it == GetStaticMap().end() ? nullptr : (it->second)->instance(s); //when I call the instance create a new object with the non trivial constructor
}

/**
 * Apply the PPSIN correction to the remnant mass
 * @param mass  mass of the remnant
 * @param s Pointer to the star that generated the remnant
 * @return The new mass of the remnant corrected for the PISN  (it can also be 0).
 * Notice we limit the mass of the remnant to be larger than the maximum Neutron star mass, e.g. we cannot produce NS from PISN correction
 */
double supernova::pisn(const double mass, Star *s) {

    //default: no correction, i.e. pisn_correction factor = 1.0;
    pisn_correction = 1.0;

    //Helium fraction
    double he_fin = s->getp(MHE::ID);
    double m_fin = s->getp(Mass::ID);
    double he_frac = he_fin/m_fin;

    //k_param (pisn_correction dependence upon helium fraction)
    double k_param = 0.67*he_frac + 0.1;

    //pulsation pair-instability supernova: pulses
    if(he_fin > 32.0 && he_fin < 64.0){
        //use Table 2 fit from Woosley 2006
        if(he_frac < 0.9){

            if(he_fin <= 37.0)
                pisn_correction = (k_param - 1.0)/5.0*he_fin + (37.0 - 32.0*k_param)/5.0;
            else if (he_fin > 37.0 && he_fin <= 60.0)
                pisn_correction = k_param;
            else
                pisn_correction = -(k_param/4.0)*he_fin + 16.0*k_param;
        }

            //use WR table 1 fit from Woosley 2016
        else{
            if(he_fin <= 37.0)
                pisn_correction = (0.5226*he_frac-0.52974)*(he_fin-32.0) + 1.0;
            else if (he_fin > 37.0 && he_fin <= 56.0){
                double val = (0.5226*he_frac-0.52974)*5.0 + 1.0;
                if(val < 0.82916)
                    pisn_correction = val;
                else
                    pisn_correction = (-0.1381*he_frac+0.1309)*(he_fin-56.0) + 0.82916;
            }
            else
                pisn_correction = -0.103645*he_fin + 6.63328;
        }

    }

    //pair-instability supernova: disintegrated
    else if (he_fin >= 64.0 && he_fin < 135.0) {
        pisn_correction = 0.0;
    }
        //end of PISN... again standard direct collapse
    else
        pisn_correction = 1.0;

    double corrected_mass=pisn_correction*mass;


    //Disable the creation of NS due to PISN correction
    if (corrected_mass<=s->get_svpar_num("sn_max_ns_mass") and pisn_correction!=1)
        return 0.0;

    return corrected_mass;
}

//TODO we can implement several versions of the following functions...
// we may want to add some input options for this
double supernova::neutrino_mass_loss(const double mass, _UNUSED Star *s) {

    double mremnant = 6.6667*(sqrt(1.0+0.3*mass)-1.0); //Lattimer & Yahil mass loss, valid for NSs...
    M_neutrinos = mass - mremnant;
    return (M_neutrinos < 0.5 ? mremnant : mass-0.5);

    //alternative (with 10% mass loss from the CO core, instead of )
    // return s->getp(MCO::ID) < 5.0 ? mass - 0.1*s->getp(MCO::ID) : mass - 0.5;

}

void supernova::main(Star *s) {

    double Mtot_before=s->getp(Mass::ID);
    double MHE_before=s->getp(MHE::ID);
    double MCO_before=s->getp(MCO::ID);

    remnant_properties(s);
    s->set_remnant();
    std::string w="";

    if (remnant_type== Lookup::Remnants::HeWD or remnant_type== Lookup::Remnants::COWD or remnant_type== Lookup::Remnants::ONeWD){
        w = utilities::log_print("WD",s,Mtot_before,MHE_before,MCO_before,Mremnant,remnant_type);
    }
    else {
        kick->apply(s);
        w = utilities::log_print("SN",s,Mtot_before,MHE_before,MCO_before,Mremnant,remnant_type,
                kick->get_random_kick(),s->vkick[3],s->vkick[0],s->vkick[1],s->vkick[2]);
    }



    s->print_to_log(w);

}

void supernova::WDformation(Star *s) {

    //TODO the following should be made consistent with the tracks of the look-up tables (we KNOW which are the stars that do not ignite He!!)
    double _Z = std::log10(s->get_Z()/0.02);
    double MHE_flash = 1.995 + 0.25*_Z + 0.087*_Z*_Z;  //Eq. 2 Hurley https://academic.oup.com/mnras/article/315/3/543/972062


    if( (s->get_zams() < MHE_flash or s->getp(MCO::ID)<1e-3) and s->getp(MHE::ID)<s->get_svpar_num("sn_Mchandra") ) {
        Mremnant = s->getp(MHE::ID); //WD mass taken as the mass of the HE core
        remnant_type = Lookup::Remnants::HeWD;
    }
    else {
        Mremnant = s->getp(MCO::ID); //WD mass taken as the mass of the CO core
        if (s->getp(MHE::ID) < 1.6) //TODO 1.6 limit is taken from Hurley (M_{c,BAGB}), Section 6 of https://academic.oup.com/mnras/article/315/3/543/972062
            remnant_type = Lookup::Remnants::COWD;
        else
            remnant_type = Lookup::Remnants::ONeWD;
    }

}

void supernova::ECSN(Star *s) {
    Mremnant = neutrino_mass_loss(s->getp(MCO::ID), s); //TODO s->getp(MCO::ID) or CO_LOWER_ECSN ?????
    remnant_type = Lookup::Remnants::NS_ECSN;
    fallback_frac = 0.; //No fallback for ECSN
}

void supernova::set_staremnant(Star *s) {

    if(remnant_type == Lookup::Remnants::NS_CCSN) {
        s->set_staremnant(new NSCCrem(s,Mremnant));
    }
    else if (remnant_type == Lookup::Remnants::NS_ECSN){
        s->set_staremnant(new NSECrem(s,Mremnant));
    }
    else if (remnant_type == Lookup::Remnants::HeWD){
        s->set_staremnant(new HeWDrem(s,Mremnant));
    }
    else if(remnant_type == Lookup::Remnants::COWD ){

        s->set_staremnant(new COWDrem(s,Mremnant));
    }
    else if(remnant_type == Lookup::Remnants::ONeWD){
        s->set_staremnant(new ONeWDrem(s,Mremnant));
    }
    else if(remnant_type == Lookup::Remnants::BH){
        s->set_staremnant(new BHrem(s,Mremnant));
    }
}

void supernova::remnant_properties(Star *s) {

    double ecsn_tshold = s->aminakedhelium() and s->get_svpar_num("sn_co_lower_ecsn_pureHe")>0 ?
            s->get_svpar_num("sn_co_lower_ecsn_pureHe")  : s->get_svpar_num("sn_co_lower_ecsn");
    double sn_tshold = s->aminakedhelium() and s->get_svpar_num("sn_co_lower_sn_pureHe")>0 ?
                         s->get_svpar_num("sn_co_lower_sn_pureHe")  : s->get_svpar_num("sn_co_lower_sn");

    ///1-Set Mremnant e Mejected
    if(s->getp(MCO::ID) < ecsn_tshold) //WD formation.. no SN
        WDformation(s);
    else if(s->getp(MCO::ID) < sn_tshold) //ECSN.. NSs
        ECSN(s);
    else //SN explosion (NSs or BHs)
        explosion(s); //it must be calculated with the values of the real star!!

    Mejected = s->amiempty() ? 0 :  s->getp(Mass::ID) - Mremnant; //Total ejected mass including neutrinos

    if (Mejected<0){
        svlog.critical("Mejected after a SN in negative: total Mass preSN="+
        utilities::n2s(s->getp(Mass::ID) ,__FILE__,__LINE__)+", Remnant mass="+
        utilities::n2s(Mremnant ,__FILE__,__LINE__),__FILE__,__LINE__);
    }

    ///2-Set the staremant object
    set_staremnant(s);

}

void supernova::explosion_SNI(Star *s) {
    s->set_empty();
    remnant_type = Lookup::Remnants::Empty;
}

void supernova::explosion_SNI(Star *s, Binstar *b) {
    s->set_empty_in_bse(b);
    remnant_type = Lookup::Remnants::Empty;
}

void supernova::set_Average_for_Unified(Star *s, double default_Average_Mremnant, double default_Average_Mejected) {

    Average_remnant = s->get_svpar_num("sn_Mremnant_average")!=-1.? s->get_svpar_num("sn_Mremnant_average") : default_Average_Mremnant;
    Average_ejected = s->get_svpar_num("sn_Mejected_average")!=-1.? s->get_svpar_num("sn_Mejected_average") : default_Average_Mejected;

}


//NSfromGau
double NSfromGau::get_NS_mass(Star *s){
   double Mrem;
   Mrem=generate_random_gau(s->get_svpar_num("sn_Mremnant_average_NS"),s->get_svpar_num("sn_Mremnant_std_NS"));
   //Limit the mass to the maximum allowed NS mass
   Mrem=std::max(1.1,std::min(Mrem,s->get_svpar_num("sn_max_ns_mass")));
   return Mrem;
}

void NSfromGau::ECSN(Star *s) {
    Mremnant = get_NS_mass(s);
    //Limit to the total mass of the star
    Mremnant = std::min(Mremnant,s->getp(Mass::ID));
    remnant_type = Lookup::Remnants::NS_ECSN;
    fallback_frac = 0.; //No fallback for ECSN
}


///Delayed
delayed::delayed(Star *s) : supernova(s) {
    if(s == nullptr) {
        Register(this, name());
    }
    else {
        //Default Mremnant, Mejected for this SN model
        double _default_Mremnant_average=1.36;
        double _default_Mejected_average=10.45;
        set_Average_for_Unified(s,_default_Mremnant_average,_default_Mejected_average);
    }
}
void delayed::explosion(Star *s) {

    double final_CO = s->getp(MCO::ID);
    double final_mass = s->getp(Mass::ID);


    double mproto, alpha_D, beta_D;

    if(final_CO < 2.5){
        mproto = 1.15;
        fallback_frac = 0.2/(final_mass - mproto);
    }
    else if(final_CO >= 2.5 && final_CO < 3.5){
        mproto = 1.15;
        fallback_frac = (0.5*final_CO - 1.05)/(final_mass - mproto);
    }
    else if(final_CO >= 3.5 && final_CO < 6){
        mproto = 1.2;
        alpha_D = 0.133 - 0.084/(final_mass - mproto);
        beta_D = 1. - 11.*alpha_D;
        fallback_frac = alpha_D*final_CO + beta_D;
    }
    else if(final_CO >= 6. && final_CO < 11.){
        mproto = 1.3;
        alpha_D = 0.133 - 0.084/(final_mass - mproto);
        beta_D = 1. - 11.*alpha_D;
        fallback_frac = alpha_D*final_CO + beta_D;
    }
    else if(final_CO >= 11.){
        mproto = 1.5;
        fallback_frac = 1.;
    }
    else{
        svlog.critical("Unexpected final_CO value in delayed::explosion",__FILE__,__LINE__,sevnstd::sn_error());
    }

    double mrem_bar_nopisn = mproto + fallback_frac*(final_mass-mproto);
    Mremnant = corrections(mrem_bar_nopisn, s); //PISN + neutrino mass loss

    ///Set Remnant type
    if (Mremnant<1e-10){
        //If Mremants is 0 after correction the SN disentegrated the stars (here we check a small number rather directly 0)
        remnant_type = Lookup::Remnants::Empty;
        s->set_empty();
    }
    else if(Mremnant >= s->get_svpar_num("sn_max_ns_mass")) //If the mass is larger than the max for a NS, this is a BH
        remnant_type = Lookup::Remnants::BH;
    else //Otherwise this is a Core Collapse NS
        remnant_type = Lookup::Remnants::NS_CCSN;


    return;


}

///Delayed with NS mass drawn from a Gaussian
delayed_gauNS::delayed_gauNS(Star *s) : supernova(s), delayed(s), NSfromGau(s){
    if(s == nullptr) {
        Register(this, name());
    }
    else {
        //Default Mremnant, Mejected for this SN model
        double _default_Mremnant_average=s->get_svpar_num("sn_Mremnant_average_NS"); //Input parameter
        double _default_Mejected_average=Average_ejected; //Same of the parent delayed class that has been already initialised
        set_Average_for_Unified(s,_default_Mremnant_average,_default_Mejected_average);
    }
}

void delayed_gauNS::explosion(Star *s){
    delayed::explosion(s);
    if (remnant_type==Lookup::Remnants::NS_CCSN or remnant_type==Lookup::Remnants::NS_ECSN){
        Mremnant = NSfromGau::get_NS_mass(s);
        Mremnant = std::min(Mremnant,s->getp(Mass::ID));
    }
}



///Rapid
rapid::rapid(Star *s) : supernova(s) {

    if ( s == nullptr) {
        Register(this, name());
    }
    else {
        //Default Mremnant, Mejected for this SN model
        double _default_Mremnant_average=1.27;
        double _default_Mejected_average=10.9;
        set_Average_for_Unified(s,_default_Mremnant_average,_default_Mejected_average);
    }

}
void rapid::explosion(Star *s) {

    double final_CO = s->getp(MCO::ID);
    double final_mass = s->getp(Mass::ID);

    double mproto = 1.1;
    double alpha_R = 0.25 - 1.275/(final_mass - mproto);
    double beta_R = 1. - 11.*alpha_R;

    if(final_CO < 2.5){
        fallback_frac = 0.2/(final_mass-mproto);
    }
    else if  (final_CO >= 2.5 && final_CO < 6.) {
        fallback_frac = (0.286*final_CO - 0.514)/(final_mass - mproto);
    }
    else if (final_CO >= 6. and final_CO < 7.){
        fallback_frac = 1;
    }
    else if ( final_CO >= 7. and final_CO < 11. ){
        fallback_frac = alpha_R*final_CO + beta_R;
    }
    else if (final_CO >= 11.){
        fallback_frac = 1;
    }
    else{
        svlog.critical("Unexpected final_CO value in delayed::explosion",__FILE__,__LINE__,sevnstd::sn_error());
    }

    double mrem_bar_nopisn = mproto + fallback_frac*(final_mass-mproto);
    Mremnant = corrections(mrem_bar_nopisn, s); //PISN + neutrino mass loss

    ///Set Remnant type
    if (Mremnant<1e-10){
        //If Mremants is 0 after correction the SN disentegrated the stars (here we check a small number rather directly 0)
        Mremnant = 0;
        remnant_type = Lookup::Remnants::Empty;
        s->set_empty();
    }
    else if(Mremnant >= s->get_svpar_num("sn_max_ns_mass")) //If the mass is larger than the max for a NS, this is a BH
        remnant_type = Lookup::Remnants::BH;
    else //Otherwise this is a Core Collapse NS
        remnant_type = Lookup::Remnants::NS_CCSN;

}

///Rapid with NS mass drawn from a Gaussian
rapid_gauNS::rapid_gauNS(Star *s) : supernova(s), rapid(s),  NSfromGau(s){
    if(s == nullptr) {
        Register(this, name());
    }
    else {
        //Default Mremnant, Mejected for this SN model
        double _default_Mremnant_average=s->get_svpar_num("sn_Mremnant_average_NS"); //Input parameter
        double _default_Mejected_average=Average_ejected; //Same of the parent delayed class that has been already initialised
        set_Average_for_Unified(s,_default_Mremnant_average,_default_Mejected_average);
    }
}

void rapid_gauNS::explosion(Star *s){
    rapid::explosion(s);
    if (remnant_type==Lookup::Remnants::NS_CCSN or remnant_type==Lookup::Remnants::NS_ECSN){
        Mremnant = NSfromGau::get_NS_mass(s);
        Mremnant = std::min(Mremnant,s->getp(Mass::ID));
    }
}




//Compactness
compactness::compactness(Star *s) : supernova(s) {
    if ( s == nullptr) {
        Register(this, name());
    }
    else{
        csi25_explosion_tshold = s->get_svpar_num("sn_compact_csi25_tshold");
        Average_Mremnant_NS    = s->get_svpar_num("sn_Mremnant_average_NS");
        Std_Mremnant_NS        = s->get_svpar_num("sn_Mremnant_std_NS");
        fallback_frac          = s->get_svpar_num("sn_compact_fallback");  //Fallback fraction in Eq. 3 Mapelli+20 (https://arxiv.org/pdf/1909.01371.pdf)
        auxiliary_table_name   = "xi25_explprobability_PS20.dat";

        //Use Average_Mremnant_NS as _default_Mremnant_average
        double _default_Mejected_average=10.45; //Same as delayed
        set_Average_for_Unified(s,Average_Mremnant_NS,_default_Mejected_average);

        //Load auxiliary table only if really needed
        if (csi25_explosion_tshold==-1 and csi25_vs_explosion_probability.empty())
            load_table(s);

    }

}

void compactness::explosion(Star *s) {

    double csi25 = csi25_mapelli20(s);
    double mproto = s->getp(MHE::ID);
    double final_mass = s->getp(Mass::ID);

    ///Following  Mapelli+20 (https://arxiv.org/pdf/1909.01371.pdf)
    //EXPLOSION
    if(triggering_explosion(csi25)){
        remnant_type = Lookup::Remnants::NS_CCSN;
        double Mgau=generate_random_gau(s->get_svpar_num("sn_Mremnant_average_NS"),s->get_svpar_num("sn_Mremnant_std_NS"));
        Mremnant = std::max(1.1,std::min(Mgau,s->get_svpar_num("sn_max_ns_mass")));
        fallback_frac = 0.; //Fallback for NS assumed 0.
    }
    //IMPLOSION
    else{
        remnant_type = Lookup::Remnants::BH;
        double mrem_bar_nopisn = mproto + fallback_frac*(final_mass-mproto);
        Mremnant = corrections(mrem_bar_nopisn, s); //PISN + neutrino mass loss
        //utilities::hardwait("Implosion",csi25,mproto,final_mass,Mremnant,__FILE__,__LINE__);
    }

    ///Check if empty
    if (Mremnant<1e-10){
        //If Mremants is 0 after correction the SN disintegrated the stars (here we check a small number rather directly 0)
        remnant_type = Lookup::Remnants::Empty;
        s->set_empty();
    }

}

double compactness::csi25_mapelli20(Star *s){
    double a=0.55, b=-1.1;
    return std::max(a+b/s->getp(MCO::ID),0.); //Mapelli+20, Eq. 2 (https://arxiv.org/pdf/1909.01371.pdf)
}

std::vector<std::vector<double>>  compactness::csi25_vs_explosion_probability={};

void compactness::load_table(Star *s){
    std::vector<std::vector<double>> _temp_Matrix=s->load_auxiliary_table(auxiliary_table_name);
    utilities::transpose(csi25_vs_explosion_probability,_temp_Matrix);
}

bool compactness::triggering_explosion(double csi25){

    if (csi25_explosion_tshold==-1){
        //Explosion probability and given xi25
        double prob_explosion = utilities::interpolate_1D(csi25,csi25_vs_explosion_probability[0],csi25_vs_explosion_probability[1],true);
        //Random number between 0 and 1
        double random_shoot   = rand_unif_0_1(utilities::mtrand);
        //utilities::hardwait("Prob",csi25,prob_explosion,random_shoot,__FILE__,__LINE__);
        //If random number is lower than the Explosion probability trigger explosion
        return random_shoot<=prob_explosion;
    }
    else if(csi25_explosion_tshold>0){
        //If lower than threshold trigger explosion
        return csi25<csi25_explosion_tshold;
    }
    else{
        svlog.critical("csi25 is negative and not -1, this is not allowes",__FILE__,__LINE__,sevnstd::params_error());
    }

    return false;

}

double compactness::generate_random_gau(double mean, double std) {
    return normal_dist(utilities::mtrand)*std + mean;
}

///Disabled
void disabled::explosion(Star *s) {


    Mremnant = s->getp(Mass::ID);

    if(Mremnant >= s->get_svpar_num("sn_max_ns_mass"))
        remnant_type = Lookup::Remnants::BH;
    else
        remnant_type = Lookup::Remnants::NS_CCSN;


    return;

}


//Direct collapse
directcollapse::directcollapse(Star *s) : supernova(s){
    if(s == nullptr) {
        Register(this, name());
    }
    //Fake values needed for Unified kick, but the kick will be always zero given that we have always a direct collapse
    Average_remnant = 1.;
    Average_ejected = 1.;

}
void directcollapse::explosion(Star *s) {

    fallback_frac = 1;
    Mremnant=s->getp(Mass::ID);

    ///Set Remnant type
    if (Mremnant<1e-10){
        //If Mremants is 0 after correction the SN disentegrated the stars (here we check a small number rather directly 0)
        Mremnant = 0;
        remnant_type = Lookup::Remnants::Empty;
        s->set_empty();
    }
    else if(Mremnant >= s->get_svpar_num("sn_max_ns_mass")) //If the mass is larger than the max for a NS, this is a BH
        remnant_type = Lookup::Remnants::BH;
    else //Otherwise this is a Core Collapse NS
        remnant_type = Lookup::Remnants::NS_CCSN;

}

/// DeathMatrix by Woosley+20 https://arxiv.org/pdf/2001.10492.pdf
std::vector<std::vector<double>> DeathMatrix::death_matrix;

DeathMatrix::DeathMatrix(Star *s) : supernova(s) {
    if(s == nullptr) {
        Register(this, name());
    } else{
        //Fake values needed for Unified kick
        Average_remnant = 1.;
        Average_ejected = 1.;
        //Check if we need to load the death matrix
        if(death_matrix.empty()){
            load_table(s);
        }
    }
}

void DeathMatrix::load_table(Star *s){
    std::vector<std::vector<double>> _temp_Matrix=s->load_auxiliary_table(auxiliary_table_name);
    utilities::transpose(death_matrix,_temp_Matrix);
}



void DeathMatrix::explosion(Star *s) {

    double MHE= s->getp(MHE::ID);

    //From DeathMatrix Tab.2 in Woosley+20
    if (MHE<=preSN_MHE_min_NS){
        Mremnant=NS_min_mass;
    }
    else if (MHE>preSN_MHE_max_PISN){
        Mremnant=0.;
    } else{
        Mremnant = utilities::interpolate_1D(MHE,death_matrix[0],death_matrix[1],false);
    }

    //TODO Can we do better?
    double fake_mproto=1.1; //Just to estimate a fallback
    ///Set Remnant type
    if (Mremnant<1e-10){
        //If Mremants is 0 after correction the SN disentegrated the stars (here we check a small number rather directly 0)
        Mremnant = 0;
        remnant_type = Lookup::Remnants::Empty;
        s->set_empty();
    }
    else if(Mremnant > NS_max_mass){   //If the mass is larger than the max for a NS, this is a BH
        remnant_type = Lookup::Remnants::BH;
        fallback_frac =(Mremnant-fake_mproto)/(s->getp(Mass::ID)-fake_mproto);
    }
    else{ //Otherwise this is a Core Collapse NS
        remnant_type = Lookup::Remnants::NS_CCSN;
        fallback_frac =(Mremnant-fake_mproto)/(s->getp(Mass::ID)-fake_mproto);
    }
}

void DeathMatrix::ECSN(Star *s) {
    explosion(s);
    remnant_type = Lookup::Remnants::NS_ECSN;
    fallback_frac = 0.0; //Assuming fallback 0 (explosion)
}
