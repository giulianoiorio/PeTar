//
// Created by spera on 30/12/18.
//

#ifndef SEVN_SUPERNOVA_H
#define SEVN_SUPERNOVA_H


#include <string>
#include <map>
#include <vector>
#include <sevnlog.h>
#include <lookup_and_phases.h>
#include <star/procs/kicks.h>
#include <random>

class Star;
using sevnstd::SevnLogging;

#define _UNUSED __attribute__ ((unused))

class supernova {

public:
    supernova(Star *s = nullptr);

    virtual ~supernova(){
        delete kick;
        kick = nullptr;
    }


    void main(Star *s); //TODO remove star s from here since it is passed in the constructor

    /**
     * Method used to initialise directy a remnant without passing through an explosion
     * @param s Pointer to the star, used to set remnant properties
     * @param Mass_remnant  Mass of the remnant
     * @param Remnant_type Type of the Remnant, it should be an elemnt of the enum Lookup::Remnant
     */
    void initialise_remnant(Star *s, double Mass_remnant, Lookup::Remnants Remnant_type);

    static supernova *Instance(std::string const &name, Star *s);

    virtual supernova *instance(_UNUSED Star *s) {
        svlog.critical("This is supposed to be pure virtual, do you have an uninitialized module?", __FILE__, __LINE__);
        return nullptr;
    }

    virtual inline std::string name() const { return "Generic Supernova"; }

    inline double get_fallback_frac() const { return fallback_frac; }

    inline double get_pisncorrection() const { return pisn_correction; }

    inline double get_M_neutrinos() const { return M_neutrinos; }

    inline double get_Mremnant() const { return Mremnant; }

    inline double get_AverageRemnant() const { return Average_remnant; }
    inline double get_AverageEjected() const { return Average_ejected; }

    inline double get_Mejected() const { return Mejected; }
    inline double get_remnant_type() const { return remnant_type; }


    /**
     * Trigger a SNIa explosion
     * @param s Pointer to star
     */
    virtual void explosion_SNI(Star *s);

    /**
     * Trigger a SNIa explosion during a BSE evolution
     * @param s Pointer to star
     * @param b Pointer to binary
     */
    virtual void explosion_SNI(Star *s, Binstar *b);

protected:
    double fallback_frac, M_neutrinos;
    double pisn_correction;
    double Mremnant, remnant_type, Mejected;
    double Average_ejected, Average_remnant;

    void Register(supernova *ptr, const std::string &_name) {
        //Register only if the name is not already in the map,
        //this is necessary to avoid multiple registration of the same model when the constructor is called from
        //inherited classes
        if (GetStaticMap().find(_name)==GetStaticMap().end()){
            GetStaticMap().insert(std::make_pair(_name, ptr));
            GetUsed().resize(GetStaticMap().size());
            GetUsed()[GetUsed().size() - 1] = 0;
        }
    }

    virtual void explosion(Star *s) = 0; //Pure virtual
    virtual void ECSN(Star *s);

    double corrections(const double mass, Star *s) {
        double mremnant = mass;
        mremnant = pisn(mremnant, s);
        mremnant = neutrino_mass_loss(mremnant, s);
        return mremnant;
    }

    double pisn(const double mass, Star *s);

    double neutrino_mass_loss(const double mass, Star *s);

    /**
     * Set the member Average_ejected and Average_remnant. These are used in the Unified Sn model
     * @param s Pointer to the star
     * @param default_Average_Mremnant  default Average_Mremnant for this SN model
     * @param default_Average_Mejected  default Average_Mejected for this SN model
     */
    void set_Average_for_Unified(Star* s, double default_Average_Mremnant, double default_Average_Mejected);


    SevnLogging svlog;



private:
    Kicks *kick;
    //std::unique_ptr<Kicks> *kick;

    static std::map<std::string, supernova *> &GetStaticMap() {
        static std::map<std::string, supernova *> _locmap;
        return _locmap;
    }

    static std::vector<int> &GetUsed() {
        static std::vector<int> _used;
        return _used;
    }


    void WDformation(Star *s);

    void remnant_properties(Star *s);

    /**
     * Use this function to properly initiliase the pointer to the object
     * It uses the protected member Mremnant
     * @param s Pointer to the star
     */
    void set_staremnant(Star *s);



};


/**
 * HOW TO ADD A NEW SN MODEL
 * A SN class model need to have the following feature:
 * - Be a class derived from class supernova or one of its child
 * - Override of the method class* instance(Star *s) is mandatory.
 * - Override of the method inline std::string name() is mandatory.
 * - Override of the method void explosion(Star *s) is mandatory (it is a pure virtual function in the base class).
 *
 * Implementation steps
 *   Assume we are implementing a new SN class  dummy
 *  -1: Define and Implement the class constructor
 *      in supernova.h:
 *          dummy(Star *s = nullptr)
 *      in supernova.cpp:
 *          dummy::dummy(Star *s) : supernova(s){
 *          if (s==nullptr)
 *              Register(this, name());
 *          else{
 *              ... *Initialise class member  parameters
 *              set Average_ejected and Average_remnant with
 *              set_Average_for_Unified(s,default_Average_Mremnant,default_Average_Mejected)
 *              where default_Average_Mremnant and default_Average_Mejected should be given.
 *              If a negative value is given, we are de-facto disabling the Unified kick model
 *              that will return an error if chosen and default_Average_Mremnatn or  default_Average_Mejected is negative or zero
 *          }
 *        }
 *      **NOTICE**:
 *          - It is important calling the parent constructor supernova
 *
 *  -2: initialise the class static object (supernova.h):
 *      static dummy _dummy;
 *      *the standard is to use as name of the member the name of the class (all lowercase)
 *      *with a an underscore as prefix.
 *
 *  -3: override the method instance (supernova.h):
 *      dummy* dummy(Star *s){
 *          return (new dummy(s));
 *       }
 *
 *  -4: override the method name (supernova.h)
 *      inline std::string name() override {return "dummy";}
 *      ** NOTICE: the chosen name is the one that will be used to set the SN model
 *                 in the input binary list.
 *
 *
 *  -5: override, define (supernova.h) and implement (supernova.cpp) the method explosion
 *      void explosion(Star*){
 *          *Here is mandatory to assign:
 *              *fallback_frac
 *              *Mremnant
 *              *remnant_type
 *      }
 *
 *   -6: Add the static class definition in static_main.h
 *      dummy dummy::_dummy;
 *
 *
 *
 *  NOTICE: The fallback_frac member is then used in the kick process to correct the natal kick. If it is not set, in the kick
 *  class it is assumed to be =0 (no correction).
 *
 */

/**
Auxiliary class to handle the random generation of NS masses
*/
class NSfromGau : virtual public supernova {
public:
    NSfromGau(_UNUSED Star *s){}
    void ECSN(Star *s) override;
    double get_NS_mass(Star *s);

protected:
    std::normal_distribution<double> normal_dist{0.,1.};
    inline double generate_random_gau(double mean, double std){
      return normal_dist(utilities::mtrand)*std + mean;
    }

};

/**
 * Delayed SN model by Fryer+12
 * BH/NS: agnostic, based on the max NS mass
 * NS mass: estimated.
 * max NS mass: from SEVN option sn_max_ns_mass
 * ECSN: from base supernova class
 * fallback: estimated
 */
class delayed : virtual public supernova{
public:
    delayed(Star *s = nullptr);

    static delayed _delayed;

    void explosion(Star *s) override;

    delayed* instance(Star *s) {
        return (new delayed(s));
    }

    inline std::string name() const override { return "delayed"; }


};
/**
 * Delayed SN model by Fryer+12 with NS mass drawn from a Maxwellian
 * BH/NS: agnostic, based on the max NS mass
 * NS mass: drawn from a Gaussian with average and std set by sn_Mremnant_average_NS and sn_Mremnant_std_NS.
 * but always larger than 1.1 Msun
 * max NS mass: from SEVN option sn_max_ns_mass
 * ECSN: NS mass drawn from a Gaussian with average and std set by sn_Mremnant_average_NS and sn_Mremnant_std_NS.
 * fallback: estimated as in delayed
 */
class delayed_gauNS :  public  delayed, NSfromGau{
public:

    delayed_gauNS(Star *s = nullptr);

    static delayed_gauNS _delayed_gauns;

    void explosion(Star *s) override;
    void ECSN(Star *s) override {NSfromGau::ECSN(s);}

    delayed_gauNS* instance(Star *s) {
        return (new delayed_gauNS(s));
    }

    inline std::string name() const override { return "delayed_gauNS"; }
};
/**
 * Rapid SN model by Fryer+12
 * BH/NS: agnostic, based on the max NS mass
 * NS mass: estimated.
 * max NS mass: from SEVN option sn_max_ns_mass
 * ECSN: from base supernova class
 * fallback: estimated
 */
class rapid : virtual public supernova{
public:
     rapid(Star *s = nullptr);

    static rapid _rapid;

    void explosion(Star *s) override;

    rapid* instance(Star *s){
        return (new rapid(s));
    }

    inline std::string name() const override { return "rapid"; }

};

/**
 * Rapid SN model by Fryer+12 with NS mass drawn from a Maxwellian
 * BH/NS: agnostic, based on the max NS mass
 * NS mass: drawn from a Gaussian with average and std set by sn_Mremnant_average_NS and sn_Mremnant_std_NS.
 * but always larger than 1.1 Msun
 * max NS mass: from SEVN option sn_max_ns_mass
 * ECSN: NS mass drawn from a Gaussian with average and std set by sn_Mremnant_average_NS and sn_Mremnant_std_NS.
 * fallback: estimated as in rapid
 */
class rapid_gauNS : public rapid, NSfromGau{
public:
    rapid_gauNS(Star *s = nullptr);

    static rapid_gauNS _rapid_gauns;

    void explosion(Star *s) override;
    void ECSN(Star *s) override {NSfromGau::ECSN(s);}

    rapid_gauNS* instance(Star *s){
        return (new rapid_gauNS(s));
    }

    inline std::string name() const override { return "rapid_gauNS"; }

};
/**
 * Compactness Model, bases on Mapelli+20 (https://arxiv.org/pdf/1909.01371.pdf).
 * BH/NS: estimated
 * NS mass: drawn from a Gaussian with average and std set by sn_Mremnant_average_NS and sn_Mremnant_std_NS.
 * but always larger than 1.1 Msun
 * max NS mass: from SEVN option sn_max_ns_mass
 * ECSN: from base supernova class
 * fallback: set to 0 for NS and get from SEVN option sn_compact_fallback for BHs
 */
class compactness : public supernova{

    compactness(Star *s = nullptr);

    static compactness _compactness;

    void explosion(Star *s) override ;

    compactness* instance(Star *s){
        return (new compactness(s));
    }

    inline std::string name() const  override { return "compact";}

protected:

    /**
     * Estimate the compactenss from the properties at the onset of the collpase
     * from Eq. 2 in Mapelli+20 (https://arxiv.org/pdf/1909.01371.pdf).
     * @param s  Pointer to star
     * @return  Compactness parameter (see Mapelli+20)
     */
    double csi25_mapelli20(Star *s);

    SevnLogging svlog;
    double csi25_explosion_tshold;
    double Average_Mremnant_NS;
    double Std_Mremnant_NS;
    std::string auxiliary_table_name;

    std::uniform_real_distribution<double> rand_unif_0_1{0.,1.};
    std::normal_distribution<double> normal_dist{0.,1.};
    inline double generate_random_gau(double mean, double std);

    static std::vector<std::vector<double>> csi25_vs_explosion_probability;

    /**
     * Estimate if the SN trigger an explosion or an implosion
     * @param csi25 Compacteness parameters csi25
     * @return true if an explosion is trigerred, false if an implosion is triggered
     * @Note the estimate method depends on the parameter csi25_explosion_tshold
     * if it is positive we just compare csi25 with csi25_explosion_tshold. If
     * csi25_explosion_tshold=-1, we estimate the likelihood of explosion at given csi25
     * (from Patton&Sukhbold20) and then we draw a random number between 0 aand 1, if it is lower
     * than per probability lkelihood a explosion is triggered.
     */
    bool triggering_explosion(double csi25);

    /**
     * Load the auxiliary table to esimate the explosion likelihood (from Patton&Sukhbold20),
     * @param s Pointer to the star (from which we get the io object)
     */
    void load_table(Star *s);


};

/**
 * Simple direct collapse model where Mrem=Mstar
 * BH/NS: agnostic, based on the max NS mass
 * NS mass: estimated
 * max NS mass: from SEVN option sn_max_ns_mass
 * ECSN: from basic supernova
 * fallback: always set to 1.0
 */
class directcollapse : public supernova{

    directcollapse(Star *s = nullptr);

    static directcollapse _directcollapse;

    void explosion(Star *s) override;

    directcollapse* instance(Star *s){
        //TODO here we use new but never delete, check
        return (new directcollapse(s));
    }

    inline std::string name() const override { return "directcollapse"; }

};

/**
 * SN model based on the Death Matrix by Woosley+20 (Tab.2, https://arxiv.org/pdf/2001.10492.pdf)
 * BH/NS: based on the max NS mass
 * NS mass: estimated
 * max NS mass: set to 2.3
 * ECSN: same as other NS
 * fallback: always set to 0.0 for NS and to 1.0 for BHs
*/
class DeathMatrix : public supernova{
public:
    DeathMatrix(Star *s = nullptr);

    static DeathMatrix _deathmatrix;

    void explosion(Star *s) override;
    void ECSN(Star *s) override;

    DeathMatrix* instance(Star *s){
        return (new DeathMatrix(s));
    }

    inline std::string name() const  override { return "deathmatrix";}

protected:

    //Some specific varaibles
    const std::string auxiliary_table_name="deathmatrix_Woosley+20.dat";  /**< name of the input file containing the death matrix values */
    const double NS_max_mass=2.30;   /**< Max NS mass (see Woosley+20) */
    const double NS_min_mass=1.24;   /**< Min NS mass (see Woosley+20) */
    const double preSN_MHE_min_NS=2.07;  /**< below this value of MHE set the mass of the remnant to NS_min_mass  */
    const double preSN_MHE_max_PISN=60.12;  /**< below this value of MHE the supernova leaves no remnant  */


    static std::vector<std::vector<double>> death_matrix;

    /**
     * Load the auxiliary table to estimate the Mrem given the preSN MHE (deathmatrix_Woosley+20.dat)
     * @param s Pointer to the star (from which we get the io object)
     */
    void load_table(Star *s);
};



/**
 * Do nothing
 */
class disabled : public supernova{

    disabled(Star *s = nullptr) : supernova(s) {
        if ( s== nullptr){
            Register(this, name());
        }
    }
    static disabled _disabled;

    void explosion(Star *s) override;

    disabled* instance(Star *s) {
        return (new disabled(s));
    }

    inline std::string name() const override { return "disabled"; }

};




#endif //SEVN_SUPERNOVA_H
