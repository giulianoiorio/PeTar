#include <static_main.h>
#include <Orbit.h> //at the very beginning so that all the properties listed are constructed before the main function
#include <property.h> //at the very beginning so that all the properties listed are constructed before the main function
#include <BinaryProperty.h>
#include <Processes.h> //at the very beginning so that all the properties listed are constructed before the main function
#include <IO.h>
#include <binstar.h>
#include <sevnlog.h>
#include <sstream>
#include <utilities.h>
#include <unistd.h>
#include <errhand.h>
using sevnstd::SevnLogging;

int main (int argc, char **argv){

    //TODO should print out a help function with the option -h

    //TODO to AVOID memory pressure the program can maybe load stars in chunks!!!
    // It will check the memory load and decide the optimal chunk of the initial conditions
    SevnLogging svlog; //TO initialise the debug level


    //TODO crea una cartella di output e puliscila ogni volta che inizia il programma
    svlog.info("Star properties = " + utilities::n2s(Property::all.size(),__FILE__,__LINE__));
    svlog.info("Binary properties = " + utilities::n2s(BinaryProperty::all.size(),__FILE__,__LINE__));
    svlog.info("Binary processes = " + utilities::n2s(Process::all.size(),__FILE__,__LINE__));


    IO sevnio(argc, argv);



    std::vector<Binstar> binaries;
    binaries.reserve(sevnio.STARS_MATRIX.size());


    for (size_t i = 0; i < sevnio.STARS_MATRIX.size(); i++)
        binaries.emplace_back(&sevnio, sevnio.STARS_MATRIX[i], i, sevnio.rand_seeds[i]);



    //TODO We use the same sevnio for all thestars. Is that a problem for parallelisation?
    //throw sevnstd::bse_error("ciao");

    /*
    try {
        svlog.critical("ciao", __FILE__, __LINE__);
    }
    catch (sevnstd::sevnerr &e){
        std::cout<<e.what()<<std::endl;
        utilities::wait();
    }
    */

    //svlog.pinfo("ciao",43,"we\n","File:",__FILE__,"\nLine:",__LINE__);
    //utilities::wait();





    //TODO handle the excpetion from the single binaries evolution
#pragma omp parallel num_threads(sevnio.nthreads) firstprivate(utilities::mtrand)
    {

#pragma omp for schedule(static)
        for (size_t i=0; i<binaries.size(); i++){


            utilities::mtrand.seed(sevnio.rand_seeds[i]);//Set random state for riproducibility


            //Star* s=binaries[i].getstar(0);
            //s->find_new_track();
            //utilities::wait();

            /*********************************
             * EVOLVE
             ********************************/
            binaries[i].recordstate(); //always record the initial state
            for(;;) {

                try {


                    binaries[i].evolve();
                    //std::cout<<binaries[i].getp(BWorldtime::ID)<<" "<<sevnio.rand_seeds[0]<<" "<< binaries[i].getstar(0)->get_rseed() <<  binaries[i].getstar(1)->get_rseed() <<std::endl;


                    if (binaries[i].breaktrigger()) //stopping condition for evolving the star
                        break; //go to evolve the next star


                    //if (binaries[i].notprint()) //Check if we have to print intermediate steps
                    //    continue;

                    if (binaries[i].isoutputtime() || binaries[i].printall()) {

                        binaries[i].recordstate_w_timeupdate();
                    }
                }
                catch(sevnstd::sevnerr& e){
                    //catch (std::exception &e){

                    //svlog.pwarning("Evolution of binary ", binaries[i].get_name(), "(ID",binaries[i].get_ID(),")",
                    //        "broken by an an error:",e.what());

                    svlog.error("Evolution of binary "+binaries[i].get_name() + "(ID "+ utilities::n2s(binaries[i].get_ID(),__FILE__,__LINE__)+
                                ") broken by an error: "+e.what(),__FILE__,__LINE__,false);

                    //std::cerr<< binaries[i].get_name()  <<" " <<e.what()<<std::endl << std::flush;
                    break;
                }

                //utilities::wait();


            }
            binaries[i].recordstate(); //always record the final state
            /********************************/


            binaries[i].print(); //print out all the recorded states for the star

            /*
            for (auto & cc: binaries[i].getstar(0)->getstate())
                svlog.debug("Star 0 state: "+utilities::n2s(cc,__FILE__,__LINE__));

            for (auto & cc: binaries[i].getstar(1)->getstate())
                svlog.debug("Star 1 state: "+utilities::n2s(cc,__FILE__,__LINE__));
            */

        }


    }


    /*
   std::cout<<svlog.get_Ndebug()<<std::endl;
   std::cout<<svlog.get_Nwarning()<<std::endl;
     */

/*
//evolve the loaded stars
    #pragma omp parallel num_threads(sevnio.nthreads)
{

    #pragma omp for schedule(static)
    for (size_t i = 0; i < binaries.size(); i++) {
        //#pragma omp critical
        //std::cout<<" Hello from thread "<<omp_get_thread_num()<<std::endl;
        //usleep(1e6);

        utilities::mtrand.seed(sevnio.rand_seeds[i]);

        binaries[i].recordstate(); //always record the initial state
        for (;;) {

            binaries[i].evolve();
            //worldtime = worldtime + stars[i].getp(Timestep::ID);

            if (binaries[i].breaktrigger()) //stopping condition for evolving the star
                break; //go to evolve the next star

            if (binaries[i].isoutputtime() || binaries[i].printall()) {
                binaries[i].recordstate();
            }
        }

        binaries[i].recordstate(); //always record the final state

        binaries[i].print(); //print out all the recorded states for the star
        //TODO maybe at this point we can also free some memory, someway, somehow
    }

}
*/


    //Print used params
    std::cerr<<std::flush;
    std::cout<<std::flush;
    std::cout<<sevnio.svpar;

    sevnio.print_params();

    return 0;

}


