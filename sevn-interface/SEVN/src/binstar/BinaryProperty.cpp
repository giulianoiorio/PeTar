//
// Created by mario on 13/02/19.
//

#include <BinaryProperty.h>
#include <Processes.h>
#include <binstar.h>
#include <star.h>

#define FOR_BPROC for(size_t _i = 0; _i < Process::all.size(); _i++)

int BinaryProperty::evolve(Binstar *binstar) {

    set_0(get());

    //double DeltaV = 0.0;
    double cumulative_V=get();

    FOR_BPROC {
        double _delta = binstar->getprocess(_i)->get((size_t)getID()); //variation caused by process "_i" on my property (ID)
        //svlog.debug("Inside BinaryProperty " + name() +" "  + utilities::n2s(get(), __FILE__, __LINE__) +"   "+utilities::n2s(_delta, __FILE__, __LINE__));
        //svlog.debug("Inside BinaryProperty " + name() +" "  + utilities::n2s(get()+_delta, __FILE__, __LINE__));
        //utilities::wait(__FILE__,__LINE__);
        //We have to check _delta+get() because is bad works on the the current value of the properties not on delta.
        //Moreover we are not checking this at then end, but for each single process so that we can tracks which brocess is causing trouble
        if(isbad(_delta+get()))
            svlog.critical("Got a bad value for [" + name() + "] from [" +binstar->getprocess(_i)->name() +  "]. Value is [" + std::to_string(_delta+get()) + "]: " +
                   "original value is [" + std::to_string(get()) +  "], delta is [" + std::to_string(_delta) + "]. \nWorldTime:" + std::to_string(binstar->getp(BWorldtime::ID)), __FILE__, __LINE__,sevnstd::bse_error());
        else
            cumulative_V += _delta;
    }


    check_boundaries(&cumulative_V, binstar);

    set(cumulative_V);


    return EXIT_SUCCESS;
}

void Semimajor::check_boundaries(double *val, Binstar *binstar) {

    //IF Semimajor is <0, force it to be equal to 0.99(R1+R2), then set force_tiny_dt to be make an extreme small step
    //where the only thing that happen is a collision and we can have a merger or a colasce.
    //GI 10/03: with the last update we disable the collision during RLO, therefore it can happen that we this method we
    //force a RLO instead of a collision
    if (*val<=0){
        svlog.pwarning("Binstar",binstar->get_name(),"ID(",binstar->get_ID(),"),""Semimajor axis after evolution is <0, it will be set to the sum of the stellar radii.",__FILE__,__LINE__);
        *val = 0.99*(binstar->getstar(0)->getp(Radius::ID) + binstar->getstar(1)->getp(Radius::ID));
        binstar->force_tiny_dt=true;
        //TODO Enable disable_a_check in this cases?
        //binstar->disable_a_check=true;
    }

}

int Eccentricity::evolve(_UNUSED Binstar *binstar){

    //If the circularised flag is true, put eccentricity to 0, then reset circularised
    //TODO I think this is not needed anymore
    if (binstar->circularised){
        set_0(get());
        double circular_ecc=0.0;
        check_boundaries(&circular_ecc,binstar); //Just to be safe
        set(circular_ecc);
        binstar->circularised=false;
    } else{
        BinaryProperty::evolve(binstar);
    }

    return EXIT_SUCCESS;
}

int dadt::evolve(Binstar *binstar) {
    set_0(get());

    //svlog.debug("A "+utilities::n2s(binstar->getp(Semimajor::ID),__FILE__,__LINE__,6)+
    //                    " A0 "+utilities::n2s(binstar->getp_0(Semimajor::ID),__FILE__,__LINE__,6)+
    //                    " T "+utilities::n2s(binstar->getp(BTimestep::ID),__FILE__,__LINE__,6));

    set((binstar->getp(Semimajor::ID) - binstar->getp_0(Semimajor::ID))/binstar->getp(BTimestep::ID));
    //svlog.debug("AR "+ utilities::n2s(get(),__FILE__,__LINE__,16),__FILE__,__LINE__);
    //utilities::wait();
    return 0;
}

int dedt::evolve(Binstar *binstar) {
    set_0(get());
    set((binstar->getp(Eccentricity::ID) - binstar->getp_0(Eccentricity::ID)) / binstar->getp(BTimestep::ID));
    return 0;
}

int BTimestep::evolve(Binstar *binstar){

    /*************************************
    * Start, set V0 and initialise max_variation
    *************************************/
    set_0(get());
    double max_variation = binstar->get_svpar_num("ts_maximum_variation");



    /*************************************
    * Check and repeat
    *************************************/

    ///Reset repeatstep
    binstar->repeatstep = false;

    /// CHeck if the Binary evolution trigger a special step
    /// in that case call a repeatstep setting a tiny dt
    if (binstar->mix or binstar->comenv or binstar->is_swallowed[0] or binstar->is_swallowed[1]){
        binstar->repeatstep = true;
        double dt=tiny_dt(binstar);
        //It could be possible that the last step done was smaller than tiny_dt so:
        dt=std::min(dt, 0.999*get());


        set(dt);
        set_0(get());
        svlog.pdebug(" Special  evolve REPEAT: ",__FILE__,__LINE__);

        return 1;
    }

    /// "Normal" check and repeat
    //TODO Check R over RLO
    if (get_0()!=binstar->get_svpar_num("ts_min_dt")){

        ///check if Semimajor axis has varied a lot, despiting the timestep control
        if (!binstar->disable_a_check){
            if(check_repeat(binstar, binstar->getp_0(Semimajor::ID), binstar->getp(Semimajor::ID), binstar->getp(dadt::ID) )){
                binstar->repeatstep = true;
                svlog.pdebug(" Semimajor axis REPEAT: ",__FILE__,__LINE__);

                //utilities::hardwait("A repeat",get(),__FILE__,__LINE__);
                return 1;

            }
        }
        else{
            binstar->disable_a_check=false; //Check skipped reset the flag
        }

        ///check if Eccentricity has varied a lot, despiting the timestep control
        if (!binstar->disable_e_check){
            if(check_repeat(binstar, binstar->getp_0(Eccentricity::ID), binstar->getp(Eccentricity::ID), binstar->getp(dedt::ID) )){
                binstar->repeatstep = true;
                svlog.pdebug(" Eccentricity REPEAT: ",__FILE__,__LINE__);
                //utilities::hardwait("Ecc repeat",get(),__FILE__,__LINE__);

                return 1;
            }
        }
        else{
            binstar->disable_e_check=false; //Check skipped reset the flag
        }

        ///Special check_repeat: Check the amount of mass transferred in the RLO
        if (!binstar->disable_DM_check){
            if(check_repeat_DM_RLO(binstar)){
                binstar->repeatstep = true;
                svlog.pdebug("DM RLO REPEAT: ",__FILE__,__LINE__);
                //utilities::hardwait("DM RLO REPEAT",get(),__FILE__,__LINE__);
                return 1;
            }
        } else {
            binstar->disable_DM_check=false; //Check skipped reset the flag
        }

        //TODO both the special check before should be considered also in proposing a new timestep
        //Otherwise we will spend a lot of time just repeating the evolution

        ///Special check repeat: OmegaSpin
        if (!binstar->disable_OmegaSPin_check and binstar->get_svpar_bool("ts_check_spin_bse")){
            if (check_repeat_OmegaSpin(binstar)){
                binstar->repeatstep = true;
                return 1;
            }
        } else{
            binstar->disable_OmegaSPin_check=false; //Check skipped reset the flag
        }


        ///Special check repeat: check the amount of spin changed for mass accreted in NS
        ///NOTICE WE DISABLE THIS BY DEFAULT BECAUSE IT CAN SLOW DOWN THE EVOLUTION
        ///AND IT IS NOT NEEDED IF WE ARE NOT INTERESTED TO THE NS SPIN
        //GI 11/03, Added parameter ts_check_Nspin to decide if use or not this check
        //THE REPOSITORY SEVN2pulsar use it
        if (!binstar->disable_OmegaRem_NS_check and binstar->get_svpar_bool("ts_check_NSspin")){
            //Activate this check only if at least of the two stars is a NS
            if (binstar->getstar(0)->amiNS() or binstar->getstar(1)->amiNS()){
                if(check_repeat_OmegaRem_NS(binstar)){
                    binstar->repeatstep = true;
                    //svlog.pdebug("dOmega  REPEAT: ",__FILE__,__LINE__);
                    //utilities::hardwait("dOmegaREPEAT",get(),__FILE__,__LINE__);
                    return 1;
                }
            }
        }
        else {
            binstar->disable_OmegaRem_NS_check=false; //Check skipped reset the flag
        }





    }



    /*************************************
    * Propose new dt
    *************************************/

    ///Initialise new_timestep as the maximum possible timestep
    double new_timestep = max_timestep(binstar);

    //If force_tiny_dt is active, set the new_dt, reset the flag and exit (force_tiny_dt is stronger than mix max dt).

    ///Check if force_tiny_dt is active.
    ///This is a special condition, it can be triggered when the cumulative binary evolution reduce a to 0
    ///and we want the next step is a step needed only for the CE evolution
    if (binstar->force_tiny_dt){
        new_timestep = tiny_dt(binstar);
        //Check
        if (new_timestep<0)
            svlog.critical("New timestep is negative or zero: "+utilities::n2s(new_timestep,__FILE__,__LINE__)
                    , __FILE__, __LINE__,sevnstd::bse_error());

        set(new_timestep);
        binstar->force_tiny_dt=false;
        return 0;
    }



    svlog.debug("Dt Max Proposed "+utilities::n2s(new_timestep,__FILE__,__LINE__),__FILE__,__LINE__);
    //utilities::wait();



    ///Condition from "expected" variation of DA
    if (binstar->getp(Semimajor::ID)==0.0)
        svlog.critical("When trying to estimate new dt. Semimajor axis is 0",__FILE__,__LINE__);
    else if(fabs(binstar->getp(dadt::ID)) != 0.0)
        new_timestep = std::min(new_timestep, max_variation*( binstar->getp(Semimajor::ID)/fabs(binstar->getp(dadt::ID)) ) );


    svlog.debug("Dt from Semimajor "+utilities::n2s(new_timestep,__FILE__,__LINE__),__FILE__,__LINE__);
    //utilities::wait();

    ///Condition from DE
    if(fabs(binstar->getp(dedt::ID)) != 0.0){ //should be consistent with the conditions in check repeat!!
        if(binstar->getp(Eccentricity::ID) > 0.1 && binstar->getp_0(Eccentricity::ID) > 0.1) {
            //Condition from expected variation of DE
            new_timestep = std::min(new_timestep, max_variation *
                                                  (binstar->getp(Eccentricity::ID) / fabs(binstar->getp(dedt::ID))));
        }
        else{
            new_timestep = std::min(new_timestep, 0.01 / fabs(binstar->getp(dedt::ID))); //limit eccentricity variation to 0.01 if e < 0.1
        }
        //Extra condition to refine when e is approching
        svlog.debug("Dt refine from eccentricity "+utilities::n2s(new_timestep,__FILE__,__LINE__),__FILE__,__LINE__);
        //utilities::wait();

    }


    svlog.debug("Final step "+utilities::n2s(new_timestep,__FILE__,__LINE__),__FILE__,__LINE__);

    for (size_t star_ID=0; star_ID<2; star_ID++){
        Star *s = binstar->getstar(star_ID);
        double DM_RL = std::abs(binstar->getprocess(RocheLobe::ID)->get_var(star_ID, Mass::ID));
        if (DM_RL>0){
            double DM_RL_dt = DM_RL/get();
            new_timestep = std::min(new_timestep, max_variation*s->getp(Mass::ID)/DM_RL_dt);
            svlog.pdebug("Dt refine from RL mass transfer ID: ",utilities::n2s(star_ID,__FILE__,__LINE__),"= "
                    ,utilities::n2s(new_timestep,__FILE__,__LINE__),DM_RL,DM_RL_dt,__FILE__,__LINE__);
        }
    }



    svlog.debug("Final step "+utilities::n2s(new_timestep,__FILE__,__LINE__,17),__FILE__,__LINE__);

    if (new_timestep<=0)
        svlog.critical("New timestep is negative or zero: "+utilities::n2s(new_timestep,__FILE__,__LINE__)
                , __FILE__, __LINE__,sevnstd::bse_error());
    else if (new_timestep<utilities::TINY)
        svlog.pwarning("Binstar",binstar->get_name(),"ID(",binstar->get_ID(),")","New dt estimated in Binary Timestep evolve  is extremely small ("+
                      utilities::n2s(new_timestep,__FILE__,__LINE__)+"). If this happens just few times, it can be ignored, but"
                                                                     "if this message it is frequent during the evoluton of a system, it "
                                                                     "can be an hint that something is broken in the bse.",__FILE__,__LINE__);


    //Check and modify if new_timestep>max_dt or new_timestep<min_dt
    check_dt_limits(new_timestep, binstar);




    set(new_timestep);
    /*********************************/


    return 0;

}

double BTimestep::max_timestep(Binstar *binstar) {

    if (binstar->get_svpar_num("ts_max_dt")>0)
        return binstar->get_svpar_num("ts_max_dt");

    return 1.0e30;
}

inline bool BTimestep::check_repeat(Binstar *binstar, const double &V0, const double &V, const double &derivative) {

    double new_dt        = max_timestep(binstar);
    double max_variation = binstar->get_svpar_num("ts_maximum_variation");
    bool shouldIrepeat = false;

    //std::cout<<" the two values = "<<V0<<"  "<<V<<std::endl;


    if(V <= 0.1 || V0 <= 0.1){
        if(fabs(V-V0) > 0.01){
            new_dt = (0.01/fabs(derivative)); //fix variation of eccentricity to 0.01 if the eccentricity is less than 0.1
            shouldIrepeat = true;
        }
    }
    else if(fabs(V-V0)/V > 2.0 * max_variation && V != 0.0){ //repeat tge step if we got twice the variation!
        new_dt = max_variation*(V/fabs(derivative));
        shouldIrepeat = true;
    }



    if(shouldIrepeat){
        //Check and modify if new_dt>max_dt or new_dt<min_dt
        check_dt_limits(new_dt, binstar);
        if(new_dt >= 0.8*get()) return false; // no need to change the timestep,
        // change it only if it is smaller than 80% of the current step

        if (new_dt<utilities::TINY)
            svlog.warning("New dt estimated in binary check and repeat   is extremely small ("+
                          utilities::n2s(new_dt,__FILE__,__LINE__)+"). If this happens just few times, it can be ignored, but"
                                                                         "if this message it is frequent during the evoluton of a system, it "
                                                                         "can be an hint that something is broken in the bse.",__FILE__,__LINE__);


        set(new_dt);
        set_0(get());
        return true;
    }
    else
        return false;

}

inline bool BTimestep::check_repeat_DM_RLO(Binstar *binstar){

    double eps_tollerance=1e-10; //Maximum Tollerance between fraction of Mass exchanged through RLO and max_variation
    double new_dt        = max_timestep(binstar);
    double max_variation = binstar->get_svpar_num("ts_maximum_variation");
    bool shouldIrepeat=false;

    for (size_t star_ID=0; star_ID<2; star_ID++) {
        Star *s = binstar->getstar(star_ID);
        double fDM_RL = std::abs(binstar->getprocess(RocheLobe::ID)->get_var(star_ID, Mass::ID)) / s->getp(Mass::ID);
        //NB We do not check fDM_RL>max_variation, because some time fDM_RL - max_variation is extremely small close
        //to the machine precision so that fDM_RL > max_variation is true, but max_variation/fDM_RL is equal to 1
        //and new_dt is the same of the old one and we have a never ending loop.
        if ((fDM_RL - max_variation) > eps_tollerance) {
            new_dt = get()*max_variation/fDM_RL;
            shouldIrepeat = true;
            break; //If true does not check the other star we already know we have to repeat
        }
    }


    if(shouldIrepeat){
        //Check and modify if new_dt>max_dt or new_dt<min_dt
        check_dt_limits(new_dt, binstar);
        if(new_dt >= 0.8*get()) return false; // no need to change the timestep,
        // change it only if it is smaller than 80% of the current step


        if (new_dt<utilities::TINY)
            svlog.warning("New dt estimated in binary check and repeat   is extremely small ("+
                          utilities::n2s(new_dt,__FILE__,__LINE__)+"). If this happens just few times, it can be ignored, but"
                                                                   "if this message it is frequent during the evoluton of a system, it "
                                                                   "can be an hint that something is broken in the bse.",__FILE__,__LINE__);

        set(new_dt);
        set_0(get());
        return true;
    }
    else
        return false;


}

bool BTimestep::check_repeat_OmegaRem_NS(Binstar *binstar) {

    double eps_tollerance=1e-10; //Maximum Tollerance between fraction of Mass exchanged through RLO and max_variation
    double new_dt        = max_timestep(binstar);
    double max_variation = binstar->get_svpar_num("ts_maximum_variation");
    bool shouldIrepeat=false;



    for (size_t star_ID=0; star_ID<2; star_ID++) {
        Star *s = binstar->getstar(star_ID);

        //Evaluate the total Omega variation due to all processes
        double DOmega=0;
        for (auto& proc : binstar->getprocesses()){
            DOmega+=proc->get_var(star_ID, OmegaRem::ID);
        }

        if (DOmega==0){
            continue;
        }

        double FOmega = std::abs(DOmega)/s->getp_0(OmegaRem::ID);
        //NB We do not check FOmega>max_variation, because some time FOmega - max_variation is extremely small close
        //to the machine precision so that FOmega > max_variation is true, but max_variation/FOmega is equal to 1
        //and new_dt is the same of the old one and we have a never ending loop.

        if ((FOmega - max_variation) > eps_tollerance) {
            new_dt = get()*max_variation/FOmega;
            shouldIrepeat = true;
            break; //If true does not check the other star we already know we have to repeat
        }
    }

    if(shouldIrepeat){
        //Check and modify if new_dt>max_dt or new_dt<min_dt
        check_dt_limits(new_dt, binstar);
        if(new_dt >= 0.8*get()) return false; // no need to change the timestep,
        // change it only if it is smaller than 80% of the current step


        if (new_dt<utilities::TINY)
            svlog.warning("New dt estimated in binary check and repeat   is extremely small ("+
                          utilities::n2s(new_dt,__FILE__,__LINE__)+"). If this happens just few times, it can be ignored, but"
                                                                   "if this message it is frequent during the evoluton of a system, it "
                                                                   "can be an hint that something is broken in the bse.",__FILE__,__LINE__);

        set(new_dt);
        set_0(get());
        return true;
    }
    else
        return false;

}

bool BTimestep::check_repeat_OmegaSpin(Binstar *binstar) {
    double eps_tollerance=1e-10; //Maximum Tollerance between fraction of Mass exchanged through RLO and max_variation
    double new_dt        = max_timestep(binstar);
    double max_variation = binstar->get_svpar_num("ts_maximum_variation");
    bool shouldIrepeat=false;

    for (size_t star_ID=0; star_ID<2; star_ID++) {
        Star *s = binstar->getstar(star_ID);


        //Notice, at this stage OmegaSpin has been already updated, while Spin that is a derived quantities not yet.
        //Special case, the initial velocity was 0, so let's start to gain angular velocity reaching
        //at most 0.1 of the critical velocity
        if (s->getp_0(OmegaSpin::ID)<eps_tollerance){

            double Spin0 = s->getp_0(AngMomSpin::ID);
            double Spin  = s->getp(AngMomSpin::ID);

            if (Spin>0.1){
                double Fspin = fabs(Spin-Spin0)/get();
                new_dt = get()*0.1/Fspin;
                shouldIrepeat = true;
                break; //If true does not check the other star we already know we have to repeat
            }

        } else{
            //Evaluate the total Omega variation due to all processes
            double DOmega=0;
            for (auto& proc : binstar->getprocesses()){
                DOmega+=proc->get_var(star_ID, AngMomSpin::ID);
            }

            if (DOmega==0){
                continue;
            }

            double FOmega = std::abs(DOmega)/s->getp_0(AngMomSpin::ID);
            //NB We do not check FOmega>max_variation, because some time FOmega - max_variation is extremely small close
            //to the machine precision so that FOmega > max_variation is true, but max_variation/FOmega is equal to 1
            //and new_dt is the same of the old one and we have a never ending loop.

            if ((FOmega - max_variation) > eps_tollerance) {
                new_dt = get()*max_variation/FOmega;
                shouldIrepeat = true;
                break; //If true does not check the other star we already know we have to repeat
            }
        }

    }

    if(shouldIrepeat){
        //Check and modify if new_dt>max_dt or new_dt<min_dt
        check_dt_limits(new_dt, binstar);
        if(new_dt >= 0.8*get()) return false; // no need to change the timestep,
        // change it only if it is smaller than 80% of the current step


        if (new_dt<utilities::TINY)
            svlog.warning("New dt estimated in binary check and repeat   is extremely small ("+
                          utilities::n2s(new_dt,__FILE__,__LINE__)+"). If this happens just few times, it can be ignored, but"
                                                                   "if this message it is frequent during the evoluton of a system, it "
                                                                   "can be an hint that something is broken in the bse.",__FILE__,__LINE__);

        set(new_dt);
        set_0(get());
        return true;
    }
    else
        return false;


}



void BTimestep::check_dt_limits(double &dt, Binstar *binstar) {

    if (binstar->get_svpar_num("ts_max_dt")>0)
        dt=std::min(dt, binstar->get_svpar_num("ts_max_dt"));

    if (binstar->get_svpar_num("ts_min_dt")>0)
        dt=std::max(dt, binstar->get_svpar_num("ts_min_dt"));

}

double BTimestep::tiny_dt(Binstar *binstar) {

    double tiny_dt;

    //Estimate the min dt to avoid to change phase due to this tiny_dt
    //If a star is remnant the min_dt will be negative therefore we set it to a large value
    double min_dt_star_0=binstar->getstar(0)->amiremnant() ? 1e30: binstar->getstar(0)->get_next_tphase()-binstar->getstar(0)->getp(Localtime::ID);
    double min_dt_star_1=binstar->getstar(1)->amiremnant() ? 1e30 : binstar->getstar(1)->get_next_tphase()-binstar->getstar(1)->getp(Localtime::ID);
    tiny_dt = std::min(0.99*std::min(min_dt_star_0,min_dt_star_1),utilities::TINY);


    //Check
    if (tiny_dt<0)
        svlog.critical("Tiny timestep is negative or zero: "+utilities::n2s(tiny_dt,__FILE__,__LINE__)
                , __FILE__, __LINE__,sevnstd::bse_error());

    return tiny_dt;
}


int BWorldtime::special_evolve(_UNUSED Binstar *binstar){

    set_0(get());
    set(get() + binstar->getp(BTimestep::ID));

    return 0;
}

int AngMom::evolve(Binstar *binstar) {

    set_0(get()); //save the previous angular momentum

    double M1 = binstar->getstar(0)->getp(Mass::ID);
    double M2 = binstar->getstar(1)->getp(Mass::ID);
    double a = binstar->getp(Semimajor::ID);
    double e = binstar->getp(Eccentricity::ID);

    double J = AngMom_from_orbit(a, e, M1, M2);

    set(J);



    return 0;

}

int Period::evolve(Binstar *binstar) {

    set_0(get());

    double a = binstar->getp(Semimajor::ID);
    double m1 = binstar->getstar(0)->getp(Mass::ID);
    double m2 = binstar->getstar(1)->getp(Mass::ID);

    set( calc_period(a, m1+m2) );

    return 0;

}

int GWtime::evolve(Binstar *binstar){

    set_0(get()); //save the previous GWtime

    double M1 = binstar->getstar(0)->getp(Mass::ID);
    double M2 = binstar->getstar(1)->getp(Mass::ID);
    double a = binstar->getp(Semimajor::ID);
    double e = binstar->getp(Eccentricity::ID);
    double a4= a*a*a*a;
    double q=1-e*e;
    double num=a4*std::sqrt(q*q*q*q*q*q*q); //a^4*(1-e^2)^(7/2)
    double den=M1*M2*(M1+M2);

    set(time_scaling*num/den);

    return 0;

}

double _RL::RL_Eg(Star* primary, Star* secondary, Binstar* b){

    return utilities::roche_lobe_Eg(primary->getp(Mass::ID),secondary->getp(Mass::ID),b->getp(Semimajor::ID));

}

int RL0::evolve(Binstar *binstar) {

    set_0(get());

    double RL_radius = RL_Eg(binstar->getstar(0), binstar->getstar(1), binstar);
    //double R = binstar->getstar(0)->getp(Radius::ID);
    //utilities::wait("RL1",RL_radius,R,__FILE__,__LINE__);

    set(RL_radius);

    return 0;

}

int RL1::evolve(Binstar *binstar) {

    set_0(get());

    double RL_radius = RL_Eg(binstar->getstar(1), binstar->getstar(0), binstar);
    //double R = binstar->getstar(1)->getp(Radius::ID);
    //utilities::wait("RL2",RL_radius,R,__FILE__,__LINE__);

    set(RL_radius);

    return 0;

}

int BEvent::special_evolve(Binstar *binstar) {

    double current_event=EventsList::NoEvent;


    ///Check events
    //If more than one event happens, the priority is given to the one with the highest value in the EventList enum
    //in lookup_and_phases

    //Check BSE Event
    current_event=binstar->get_event_from_processes();
    //Check SSE event
    if (binstar->getstar(0)->getp(Event::ID)!=-1){
        current_event=std::max(binstar->getstar(0)->getp(Event::ID),current_event);
    }
    else if(binstar->getstar(1)->getp(Event::ID)!=-1){
        current_event=std::max(binstar->getstar(1)->getp(Event::ID),current_event);
    }

    //Set current event
    set(current_event);
    //Reset events, if an event just happened reset all the events
    if (current_event!=EventsList::NoEvent)
        binstar->reset_events();

    return EXIT_SUCCESS;
}
