#pragma once

#include <string>
#include <vector>
#include <eccodes.h>

namespace pifo {
    class GribFile {
    public:
        struct GribFieldDescription  {
            std::string name;
            std::string shortName;
            std::string typeOfLevel;
            std::vector<long> levelList;
            long ni;
            long nj;
            long numberOfPoints;
        };

        GribFile(const std::string& fn);

        ~GribFile();
        
        std::vector<GribFieldDescription> getFieldList();

        void getData(const std::string& field, const std::string& levelType, long level, double* data);

        void getLatitudes(const std::string& field, const std::string& levelType, long level, double* data);

        void getLongitudes(const std::string& field, const std::string& levelType, long level, double* data);

    private:
        std::string gribfile;
        codes_context* context = nullptr;
        codes_index* fieldNamesIndex = nullptr;

        void loadFieldIndex();

        codes_handle* selectHandle(const std::string& field, const std::string& levelType, long level);

        void getDataOfKey(const std::string& field, const std::string& levelType, long level, const std::string& key, double* data);

    };
}