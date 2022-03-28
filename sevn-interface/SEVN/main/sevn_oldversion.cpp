#include <static_main.h>
#include <property.h> //at the very beginning so that all the properties listed are constructed before the main function
#include <star.h>
using sevnstd::SevnLogging;

int main (int argc, char **argv){

    //TODO should print out a help function with the option -h

    //TODO to AVOID memory pressure the program can maybe load stars in chunks!!!
    // It will check the memory load and decide the optimal chunk of the initial conditions


    //TODO crea una cartella di output e puliscila ogni volta che inizia il programma
    cout<<"Main file = "<<Property::all.size()<<endl;
    SevnLogging svlog; //TO initialise the debug level
    IO sevnio(argc, argv);

    std::vector<Star> stars;
    stars.reserve(sevnio.STARS_MATRIX.size());




    for (size_t i = 0; i < sevnio.STARS_MATRIX.size(); i++)
        //stars.emplace_back(Star(&sevnio, sevnio.STARS_MATRIX[i], i));
        // GI 71119: with emplace_Back we can use the variadic property to put directly the values to initialise the object Star. In thiw way
        // the object will bi directly instanciated instead of create an obejct and then move it as happens with push_back().
        stars.emplace_back(&sevnio, sevnio.STARS_MATRIX[i], i);


//evolve the loaded stars
#pragma omp parallel num_threads(sevnio.nthreads) firstprivate(utilities::mtrand)
    {

#pragma omp for schedule(static)
        for (size_t i = 0; i < stars.size(); i++) {
            //#pragma omp critical
            //std::cout<<" Hello from thread "<<omp_get_thread_num()<<std::endl;
            //usleep(1e6);


            utilities::mtrand.seed(sevnio.rand_seeds[i]);


            //svlog.debug("Wordltime "+utilities::n2s(stars[i].getp(Worldtime::ID),__FILE__,__LINE__));
            //svlog.debug("Next "+utilities::n2s(stars[i].getp(NextOutput::ID),__FILE__,__LINE__));
            //utilities::wait();

            stars[i].recordstate(); //always record the initial state
            //svlog.debug("Wordltime "+utilities::n2s(stars[i].getp(Worldtime::ID),__FILE__,__LINE__));
            //svlog.debug("Next "+utilities::n2s(stars[i].getp(NextOutput::ID),__FILE__,__LINE__));
            //utilities::wait();

            for (;;) {

                stars[i].evolvestar();

                if (stars[i].breaktrigger()) //stopping condition for evolving the star
                    break; //go to evolve the next star

                svlog.debug("Wordltime "+utilities::n2s(stars[i].getp(Worldtime::ID),__FILE__,__LINE__));
                svlog.debug("Next "+utilities::n2s(stars[i].getp(NextOutput::ID),__FILE__,__LINE__));
                //utilities::wait("Waiting point");

                if (stars[i].notprint()) //Check if we have to print intermediate steps
                    continue;
                else if (stars[i].isoutputtime() || stars[i].printall()) {
                    stars[i].recordstate_w_timeupdate(); //Record state and estimate next output time
                }


            }
            stars[i].recordstate(); //always record the final state

            stars[i].print(); //print out all the recorded states for the star
            //TODO maybe at this point we can also free some memory, someway, somehow
        }

    }

    std::cout<<svlog.get_Ninfo()<<std::endl;

    return 0;

}

