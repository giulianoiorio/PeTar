//
// Created by spera on 16/06/20.
//

#ifndef SEVN_KICKS_H
#define SEVN_KICKS_H

#include <string>
#include <map>
#include <vector>
#include <sevnlog.h>
#include <random>
#include <utilities.h>

class Star;
using sevnstd::SevnLogging;

//TODO Kicks here should be restored, Kicks should be general enough to not contain gaussian_265
//then we can create a son class Maxwellian with methods that can be used by Hobbs Unified etc
class Kicks {

public:

    Kicks() : standard_gaussian(0,1), uniformRealDistribution(0,1){
    }

    virtual ~Kicks(){
    }

    static Kicks *Instance(std::string const &name);
    virtual Kicks *instance(){ svlog.critical("This is supposed to be pure virtual, do you have an uninitialized module?", __FILE__, __LINE__); return nullptr; }
    virtual inline std::string name() { return "Generic SN-kick model"; }

    /**
     * Wrapper for specified _apply functions
     * @param s
     */
    void apply(_UNUSED Star *s);

    virtual void _apply(_UNUSED Star *s) = 0;

    inline double get_random_kick(){
        if (random_velocity_kick==utilities::NULL_DOUBLE)
            svlog.critical("Trying to get random_kick before it is initiliased",__FILE__,__LINE__);

        return random_velocity_kick;
    }

    
protected:

    void Register(Kicks *ptr, const std::string &_name) {
        GetStaticMap().insert(std::make_pair(_name, ptr));
        GetUsed().resize(GetStaticMap().size());
        GetUsed()[GetUsed().size()-1] = 0;
    }

    std::normal_distribution<double> standard_gaussian;
    std::uniform_real_distribution<double> uniformRealDistribution;



    inline void set_random_kick(const double& a){
        random_velocity_kick =a;
        return;
    }
    void kick_initializer();
    SevnLogging svlog;
    double random_velocity_kick=utilities::NULL_DOUBLE;

    /**
     * Check if we have to make correction to the final Vkick (after all the fallback and similar correction)
     * It checks:
     *      - If Mremant=0 (e.g. after a PPISN)
     * So far it just check that the final Vkick is not lower than the parameter sn_min_vkick. If this is the case
     * it just sets the final velocity to the minimum value and rescale all the components by the factor min_vkick/old_vkick
     * if old_vkick is 0, new isotropic velocity components are randomly drawn
     * @param s
     */
    virtual void check_and_correct_vkick(Star* s);


    inline double draw_from_gaussian(double std, double mean=0.){
        return std*standard_gaussian(utilities::mtrand) + mean;
    }

private:

    static std::map<std::string, Kicks *> & GetStaticMap(){
        static std::map<std::string, Kicks *> _locmap;
        return _locmap;
    }
    static std::vector<int> & GetUsed(){
        static std::vector<int> _used;
        return _used;
    }


};


class Hobbs : public Kicks{

    Hobbs(bool reg = true){
        if(reg) {
            Register(this, name());
        }
    }

    static Hobbs _hobbs;

    void _apply(Star *s) override;

    //when I use the regist function I use the non trivial constructor of the Hermite6th class
    Hobbs* instance() {
        return (new Hobbs(false));
    }

    inline std::string name() override { return "hobbs"; }

};

/**
 * Kick as in models alpha in Giacobbo&Mapelli2018  https://ui.adsabs.harvard.edu/abs/2018MNRAS.480.2011G/abstract
 * The ECSN receives a kick from a Maxwellian with sigma= 15 km/s the CCSN from a Maxwellian with sigma= 265 km/s,
 * the BH kick is then corrected for the fallback fraction vkick_corrected =(1−ffb)vkick_maxwellian,
 */
class EC15CC265 : public Kicks{

    EC15CC265(bool reg = true) : gaussian_ecsn(0,sigma_ecsn), gaussian_ccsn(0,sigma_ccsn){
        if(reg) {
            Register(this, name());
        }
    }

    static EC15CC265 _ec15cc265;

    void _apply(Star *s) override;

    //when I use the regist function I use the non trivial constructor of the Hermite6th class
    EC15CC265* instance() {
        return (new EC15CC265(false));
    }

    inline std::string name() override { return "ec15cc265"; }

protected:

    const double sigma_ecsn=15.0; /*!< dispersion of the Maxwellian used to draw the natal kick of ECSN remnant*/
    const double sigma_ccsn=256.0; /*!< dispersion of the Maxwellian used to draw the natal kick of CCSN remnant*/
    //Make them static? Maybe
    std::normal_distribution<> gaussian_ecsn;
    std::normal_distribution<> gaussian_ccsn;
};

class Unified : public Kicks{

    Unified(bool reg = true){
        if(reg) {
            Register(this, name());
        }
    }

    static Unified _unified;

    void _apply(Star *s) override;

    Unified* instance() {
        return (new Unified(false));
    }

    inline std::string name() override { return "unified"; }

};

class Zeros : public Kicks{

public:
    Zeros(bool reg = true){
        if(reg) {
            Register(this, name());
        }
    }

    static Zeros _zeros;

    void _apply(Star *s) override;

    Zeros* instance() {
        return (new Zeros(false));
    }

    inline std::string name() override { return "zeros"; }

};

/**
 * Kick as in model CC15 by Giacobbo&Mapelli2018  https://ui.adsabs.harvard.edu/abs/2018MNRAS.480.2011G/abstract
 * Both ECSN and CCSN receive a kick drawn from a Maxwellian with sigma= 15 km/s,
 * the BH kick is then corrected for the fallback fraction vkick_corrected =(1−ffb)vkick_maxwellian,
 */
class CC15 : public Kicks{

public:

    CC15(bool reg = true) : gaussian15(0,sigma){
        if(reg) {
            Register(this, name());
        }
    }

    CC15* instance() {
        return (new CC15(false));
    }

    static CC15 _cc15;
    void _apply(Star *s) override;
    inline std::string name() override { return "cc15";}

protected:

    const double sigma=15.0; /*!< dispersion of the Maxwellian used to draw the natal kick*/
    std::normal_distribution<> gaussian15;
};



#endif //SEVN_KICKS_H
