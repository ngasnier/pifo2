#pragma once

#include "util/Factory.h"
#include "app/Application.h"

namespace pifo {
    namespace app {
        class ApplicationFactory : public Factory<Application>
        {
        public:
            ApplicationFactory()
            {

            }

        };
    }
}