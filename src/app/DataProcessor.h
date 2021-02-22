#pragma once

#include "Application.h"

namespace pifo {
    namespace app  {
        class DataProcessor : public Application {
        public:
            DataProcessor() : Application()
            {

            }

            virtual void run();
        };

    }
}