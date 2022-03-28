//
// Created by Giuliano Iorio on 2021-03-07.
//

#ifndef SEVN_QCRIT_H
#define SEVN_QCRIT_H

#include <utilities.h>
#include <map>


/**
 *
 * HOW TO ADD A NEW QCRIT IMPLEMENTATION
 *
 *
 *
 * 1- Write the definition of the new qcrit inside the struct Qcrit in qcrit.h.
 * The default template name has a  qcrit as prefix and separate names using _
 * e.g.  double qcrit_Hurley(Star* donor, _UNUSED Star* accretor);
 *
 *
 * 2- Implement the method in qcrit.cpp
 * e.g. double Qcrit::qcrit_Hurley(Star *donor, Star *accretor){...}
 *
 * 3- Add new elements in the map qcritmap in qcrit.cpp
 * The  first map element is a name or alias to be used to set this particular qcrit option in the parameters,
 * the second is the address to the method inside Qcrit. The default style is to use a long and short name version
 * for the same qcrit option.
 * e.g.
 *         {"Hurley", &Qcrit::qcrit_Hurley},
 *         {"H",&Qcrit::qcrit_Hurley},
 *
 * this particular qcrit can be called with -rlo_qcrit Hurley or -rlo_qcrit H.
 *
 *Notice that, considering the SEVN implementation, all the stellar properties inside the function should be called with getp_0 instead of getp.
 *
 *
 */



class Star;

struct Qcrit{

    ///Constructors
    Qcrit(){}
    Qcrit(double _qc_const) : qc_const(_qc_const){}


    ///Methods
    //Constants
    double qc_const;
    /**
     * Constant qc
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @return qc of the donor star
     */
    inline double qcrit_constant(_UNUSED Star* donor, _UNUSED Star* accretor){ return qc_const;}

    //Hurley
    /**
     * qc from Hurley+02
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @return qc of the donor star
     */
    double qcrit_Hurley(Star* donor, _UNUSED Star* accretor);


    //Hurley + Webbink
    /**
     * qc from Hurley+02 with Webbink (1988) equations.
     * Note: this function has been taken directly from BSE, the comments highlights
     * some difference (if present) in what is written in BSE and what is reported in Hurley+02
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @return qc of the donor star
     */
    double qcrit_Hurley_Webbink(Star* donor, _UNUSED Star* accretor);

    //Hurley + Webbink + Shao
    /**
     * qcrit from Hurley+02 with Webbink (1988) equations for stars in the RGB phase, plus the Shao+21 (https://arxiv.org/pdf/2107.03565.pdf)
     * formalism for Star+BH mass transfer stability
     * @param donor donor Pointer to the donor star
     * @param accretor Pointer to the accretor star
     * @return qc of the donor star
     */
    double qcrit_Hurley_Webbink_Shao(Star* donor, Star *accretor);

    //COMPAS-Neijssel
    /**
     * Use the COMPAS prescriptions taken from Neijssel+2019, section 2.3
     * This implementation is taken directly from COSMIC (https://github.com/COSMIC-PopSynth/COSMIC/blob/develop/cosmic/src/evolv2.f).
     * In the documentation they report: "We convert from radial response to qcrit for MS and HG,
     *  which assumes conservative mass transfer,
     *  Stable MT is always assumed for stripped stars,
     *  Assume standard qcrit from BSE for kstar>=10"
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @return qc of the donor star
     */
    double qcrit_Neijssel(Star* donor, _UNUSED Star* accretor);

    /**
     * Default qcrit used in COMPAS, we took these values direcly from the default OPTION in COMPAS
     * (https://github.com/TeamCOMPAS/COMPAS/blob/dev/src/Options.cpp, starting from row 352).
     * Notice in COMPAS they define q=Maccretor/Mdonor, while we define q=Mdonor/Maccretor so qcrit_SEVN=1/qcrit_COMPAS
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @return qc of the donor star
     */
    double qcrit_COMPAS(Star* donor, Star* accretor);

    /**
     * Claeys+14 qcrit as used in COSMIC (https://github.com/COSMIC-PopSynth/COSMIC/blob/develop/cosmic/src/evolv2.f row 2041)
     * Notice that their values are not exactly the same of Tab. 2 in Claeys (even considering that qc_Cosmic = 1/qc_Claeys).
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @return qc of the donor star
     */
    double qcrit_COSMIC_Claeys(Star* donor, Star* accretor);


    /** AUXILIARY **/
    /**
     * Qcrit for giants (BSE type 3,5,6) by Webbink (1988) from models
     * for condensed polytropes.
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @return qc of the donor star
     */
    double _Webbink_giant(Star* donor, _UNUSED Star* accretor);

    /**
     * Qcrit for giants (BSE type 3,5,6) by Hurley +02:
     * qc = (1.67-zpars(7)+2.0*pow(mcore1/m1),5.))/2.13 Eq.57
     * with zpar from Eq. 47 Hurley+00
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @return qc of the donor star
     */
    double _Hurley_giant(Star* donor, _UNUSED Star* accretor);


};

typedef double (Qcrit::*qcrit_method_t)(Star *donor, Star *accretor);
typedef std::map<std::string, qcrit_method_t> QCRITMAP;
extern const QCRITMAP qcritmap;








#endif //SEVN_QCRIT_H
