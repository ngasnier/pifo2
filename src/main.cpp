#include "atlas/atlas.h"
#include "atlas/library/Library.h"
#include "atlas/runtime/Log.h"
#include "atlas/parallel/mpi/mpi.h"
#include "eckit/runtime/Tool.h"
#include "app/Application.h"
#include "app/DataProcessor.h"
#include "app/ModelRun.h"
#include "app/ApplicationFactory.h"

// #include <eckit/config/YAMLConfiguration.h>
namespace pifo {
    using namespace std;
    using namespace atlas;
    using namespace atlas::mpi;

    class Pifo : public eckit::Tool {
    public:
        Pifo( int argc, char** argv ) : eckit::Tool( argc, argv ) {};

        ~Pifo() override = default;

        void run() override
        {            
            // Init Atlas & affichage informations
            atlas::Library::instance().initialise( argc(), argv() );

            Log::info() << "PIFO2\n" << endl;
            Log::info() << Library::instance().information() << endl << endl;
            Log::info() << "  MPI tasks: " << mpi::comm().size() << endl;
            Log::info() << "  OpenMP threads per MPI task: " << atlas_omp_get_max_threads() << endl;
            Log::info() << endl;

            pifo::app::ApplicationFactory appFactory;
            appFactory.registerType<pifo::app::DataProcessor>("dataprocessor");
            appFactory.registerType<pifo::app::ModelRun>("run");

            auto app = appFactory.create("run");
            app->run();

            atlas::Library::instance().finalise();            
        }
    };
}

int main( int argc, char* argv[] ) {
    pifo::Pifo pifo(argc, argv);
    pifo.start();
    return 0;
}