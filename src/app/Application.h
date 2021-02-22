#pragma once

namespace pifo {
    namespace app  {
        class Application {
        public:
            Application();

            virtual void run() = 0;
        };

    }
}