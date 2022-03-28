//
// Created by spera on 16/06/20.
//

//TODO SNKICKS and KICKS should be merged

#include <star/procs/kicks.h>
#include <utilities.h>
#include <star.h>
#include <supernova.h>


Kicks* Kicks::Instance(std::string const &name){
    auto it = GetStaticMap().find(name);
    if(it!= GetStaticMap().end()) GetUsed()[std::distance(GetStaticMap().begin(), it)] = 1;
    return it == GetStaticMap().end() ? nullptr : (it->second)->instance(); //when I call the instance create a new object with the non trivial constructor
}

void Kicks::check_and_correct_vkick(Star* s) {

    //Check if Mremant 0 (this happens after a PPSIN)
    //In this case set the kick to 0
    if(s->get_supernova()->get_Mremnant()==0){
        s->vkick[3] = s->vkick[2] = s->vkick[1] = s->vkick[0] = 0.0;
        return;
    }

    //Check if below Vmin and correct
    if (s->vkick[3]<s->get_svpar_num("sn_min_vkick")){
        double new_kick = s->get_svpar_num("sn_min_vkick");

        //if minimum vkick is 0 just set all the components to 0
        if (new_kick==0){
            s->vkick[3] = s->vkick[2] = s->vkick[1] = s->vkick[0] = 0.0;
        }
        //if minimum vkick>0 just but old vkick=0, set the new module equal to vkick and draw new isotropic direction for the components
        else if(s->vkick[3]==0){
            //Generate random teta angle
            double random = uniformRealDistribution(utilities::mtrand);
            double teta = 2.0 * M_PI * random;

            s->vkick[3] = new_kick; //Total velocity
            s->vkick[2] = 2.0 * new_kick * random - new_kick;  //Vz (random between -Vkick and Vkick
            double vel2d = sqrt(new_kick*new_kick - s->vkick[2]*s->vkick[2]);
            s->vkick[0] = vel2d * cos(teta);  // Vx
            s->vkick[1] = vel2d * sin(teta);  // Vy

        }
        //if minimum vkick>0 just and old vkick>0, set the new module equal to vkick and rescale all the components by a factor new_module/old_module
        else{
            double correction=new_kick/s->vkick[3];
            //Correct so that vkick=vmin
            for (auto& vkick : s->vkick)
                vkick*=correction;
        }
    }
}

void Kicks::apply(_UNUSED Star *s){
    //If star is empty after a SN kick always 0
    _apply(s);
    //Check and correct vkick
    check_and_correct_vkick(s);
    return;
}

void Hobbs::_apply(Star *s) {

    double hobbs_std = s -> get_svpar_num("sn_kick_velocity_stdev");

   //maxwellian velocities with sigma = get_svpar_num("sn_kick_velocity_sigma"); usual: 265 km/s, as in Hobbs et al., 2005, MNRAS 360 974
   //all kicks are rescaled accordingly to the fallback fraction (direct collapse = no kicks)
    s->vkick[0] = draw_from_gaussian(hobbs_std);
    s->vkick[1] = draw_from_gaussian(hobbs_std);
    s->vkick[2] = draw_from_gaussian(hobbs_std);
    s->vkick[3] = random_velocity_kick = sqrt( s->vkick[0]*s->vkick[0] + s->vkick[1]*s->vkick[1] + s->vkick[2]*s->vkick[2]);

    //Correct
    double f_fb_correction =  s->get_supernova()->get_fallback_frac() != -1 ? (1.0 - s->get_supernova()->get_fallback_frac()) : 1.0;
    for (auto& vkick : s->vkick)
        vkick*=f_fb_correction;


}

void Unified::_apply(Star *s) {

    ///Initial check
    if (s->get_supernova()->get_AverageEjected()<0)
        svlog.critical("Average Ejected Mass in Unified kick is negative. It is likely not initialised  in your chosen SN model ("
                       +s->get_supernova()->name()+").",__FILE__,__LINE__,sevnstd::sn_error());
    else if (s->get_supernova()->get_AverageRemnant()<0)
        svlog.critical("Average remnant Mass in Unified kick is negative. It is likely not initialised  in your chosen SN model ("
                       +s->get_supernova()->name()+").",__FILE__,__LINE__,sevnstd::sn_error());

    double unified_std = s -> get_svpar_num("sn_kick_velocity_stdev");

    // Neutrino mass loss is included in s->get_supernova()->get_Mejected() = M_star_final - Mremnant, i.e.  s->get_supernova()->get_Mejected() is ALWAYS != 0.
    // Note: for ECSNe, Mejected = MCO - Mremnant, i.e. kicks are != 0 because of neutrino mass loss (in principle, kicks != 0 should come from the ejection of tiny envelopes around ONe-WD)

    double ejected_mass = s->get_supernova()->get_fallback_frac() == 1 ? 0.0 : (s->get_supernova()->get_remnant_type() == Lookup::Remnants::NS_ECSN ? s->getp(MCO::ID) - s->get_supernova()->get_Mremnant() : s->get_supernova()->get_Mejected()); //do not consider neutrinos in case of direct collapse
    double correction = (ejected_mass/s->get_supernova()->get_AverageEjected()) * (s->get_supernova()->get_AverageRemnant()/s->get_supernova()->get_Mremnant());

    s->vkick[0] = draw_from_gaussian(unified_std);
    s->vkick[1] = draw_from_gaussian(unified_std);
    s->vkick[2] = draw_from_gaussian(unified_std);
    s->vkick[3] = sqrt( s->vkick[0]*s->vkick[0] + s->vkick[1]*s->vkick[1] + s->vkick[2]*s->vkick[2]);
    //Velocity kick before correction
    random_velocity_kick = s->vkick[3];
    //Correct
    for (auto& vkick : s->vkick)
        vkick*=correction;

}

void Zeros::_apply(Star *s) {
    random_velocity_kick = s->vkick[0] = s->vkick[1] = s->vkick[2] = s->vkick[3] = 0.;
}

void EC15CC265::_apply(Star *s) {

    //If ECSN explosion, use the ECSN Maxwellian, for CCSN use the CCSN Maxwellian
    //all kicks are rescaled accordingly to the fallback fraction (direct collapse = no kicks)

    if(s->get_supernova()->get_remnant_type()==Remnants::NS_ECSN){
        s->vkick[0] = gaussian_ecsn(utilities::mtrand);
        s->vkick[1] = gaussian_ecsn(utilities::mtrand);
        s->vkick[2] = gaussian_ecsn(utilities::mtrand);
        s->vkick[3] = random_velocity_kick =sqrt(s->vkick[0]*s->vkick[0] + s->vkick[1]*s->vkick[1] + s->vkick[2]*s->vkick[2]);
    } else{
        s->vkick[0] = gaussian_ccsn(utilities::mtrand);
        s->vkick[1] = gaussian_ccsn(utilities::mtrand);
        s->vkick[2] = gaussian_ccsn(utilities::mtrand);
        s->vkick[3] = random_velocity_kick =sqrt(s->vkick[0]*s->vkick[0] + s->vkick[1]*s->vkick[1] + s->vkick[2]*s->vkick[2]);
    }

    if ( (s->get_supernova()->get_remnant_type()==Remnants::BH)){
        double f_fb_correction =  s->get_supernova()->get_fallback_frac() != -1 ? (1.0 - s->get_supernova()->get_fallback_frac()) : 1.0;
        for (auto& v : s->vkick)
            v=v*f_fb_correction;
    }

}

void CC15::_apply(Star *s) {
    //maxwellian velocities with sigma = 265 km/s, as in Hobbs et al., 2005, MNRAS 360 974
    //all kicks are rescaled accordingly to the fallback fraction (direct collapse = no kicks)

    s->vkick[0] = gaussian15(utilities::mtrand);
    s->vkick[1] = gaussian15(utilities::mtrand);
    s->vkick[2] = gaussian15(utilities::mtrand);
    s->vkick[3] = random_velocity_kick =sqrt(s->vkick[0]*s->vkick[0] + s->vkick[1]*s->vkick[1] + s->vkick[2]*s->vkick[2]);

    if ( (s->get_supernova()->get_remnant_type()==Remnants::BH)){
        double f_fb_correction =  s->get_supernova()->get_fallback_frac() != -1 ? (1.0 - s->get_supernova()->get_fallback_frac()) : 1.0;
        for (auto& v : s->vkick)
            v=v*f_fb_correction;
    }
}


