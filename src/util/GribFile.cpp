#include <stdexcept>
#include "GribFile.h"

namespace pifo {

    GribFile::GribFile(const std::string& fn)
    {
        gribfile = fn;
        context = codes_context_get_default();
        loadFieldIndex();
    }

    GribFile::~GribFile()
    {
        if (fieldNamesIndex) codes_index_delete(fieldNamesIndex);
    }
    
    std::vector<GribFile::GribFieldDescription> GribFile::getFieldList()
    {
        char name[256];
        int err;
        size_t name_size;
        size_t type_size;
        size_t level_size;
        char** names = nullptr;
        char** types = nullptr;
        long* levels = nullptr;
        std::string error_message;
        std::vector<GribFile::GribFieldDescription> fieldList;
        GribFieldDescription field;

        // First get the short names index 
        if ((err=codes_index_get_size (fieldNamesIndex, "shortName", &name_size)))
        {
            error_message = "codes_index_get_size error "+ std::to_string(err);
        }
        else
        {
            names = new char*[name_size];
            if ((err=codes_index_get_string(fieldNamesIndex, "shortName", names, &name_size)))
            {
                error_message = "codes_index_get_string error "+ std::to_string(err);
            }
        }
        
        // Get the typeOfLevel index
        if (!err)
        {
            if ((err=codes_index_get_size (fieldNamesIndex, "typeOfLevel", &type_size)))
            {
                error_message = "codes_index_get_size error "+ std::to_string(err);
            }
            else
            {
                types = new char*[type_size];
                if ((err=codes_index_get_string(fieldNamesIndex, "typeOfLevel", types, &type_size)))
                {
                    error_message = "codes_index_get_string error "+ std::to_string(err);
                }
            }
        }

        // Get the level index
        if (!err)
        {
            if ((err=codes_index_get_size (fieldNamesIndex, "level", &level_size)))
            {
                error_message = "codes_index_get_size error "+ std::to_string(err);
            }
            else
            {
                levels = new long[level_size];
                if ((err=codes_index_get_long(fieldNamesIndex, "level", levels, &level_size)))
                {
                    error_message = "codes_index_get_string error "+ std::to_string(err);
                }
            }
        }


        if (!err)
        {
            for (size_t i=0;i<name_size;i++)
            {
                if ((err=codes_index_select_string(fieldNamesIndex, "shortName", names[i])))
                {
                    error_message = "codes_index_select_string error "+ std::to_string(err);
                    break;
                }

                for (size_t j=0;j<type_size;j++)
                {
                    if ((err=codes_index_select_string(fieldNamesIndex, "typeOfLevel", types[j])))
                    {
                        error_message = "codes_index_select_string error "+ std::to_string(err);
                        break;
                    }

                    field = GribFieldDescription();
                    field.shortName = std::string(names[i]);
                    field.typeOfLevel = std::string(types[j]);


                    size_t k = 0;
                    do {

                        if ((err=codes_index_select_long(fieldNamesIndex, "level", levels[k])))
                        {
                            error_message = "codes_index_select_long error "+ std::to_string(err);
                            break;
                        }

                        codes_handle* handle = codes_handle_new_from_index(fieldNamesIndex, &err);
                        if (err && err!=CODES_END_OF_INDEX)
                        {
                            error_message = "codes_handle_new_from_index error "+std::to_string(err);
                            break;
                        }

                        if (!err)
                        {

                            size_t len;
                            if ((err = codes_get_length(handle, "name", &len)))
                            {
                                error_message = "codes_get_length error "+std::to_string(err);
                                break;
                            }
                            else
                            {
                                if ((err=codes_get_long(handle, "Ni", &field.ni)))
                                {
                                    error_message = "codes_get_long error "+std::to_string(err);
                                    break;
                                }

                                if ((err=codes_get_long(handle, "Nj", &field.nj)))
                                {
                                    error_message = "codes_get_long error "+std::to_string(err);
                                    break;
                                }

                                if ((err=codes_get_long(handle, "numberOfPoints", &field.numberOfPoints)))
                                {
                                    error_message = "codes_get_long error "+std::to_string(err);
                                    break;
                                }

                                if ((err=codes_get_string(handle, "name", name, &len)))
                                {
                                    error_message = "codes_get_string error "+std::to_string(err);
                                    break;
                                }
                                else
                                {
                                    field.name = std::string(name);
                                    field.levelList.push_back(levels[k]);
                                }                                        
                            }

                            codes_handle_delete(handle);                            
                        }

                        k++;
                    } while (k<level_size);

                    if (err && err!=CODES_END_OF_INDEX) break;
                    err = 0; // Need to reset error in order not to break on upper loop

                    if (field.levelList.size()>0) fieldList.push_back(field);
                }

                if (err) break;
            }
        }

        // Cleanup and error management

        if (names) delete names;
        if (types) delete types;
        if (levels) delete levels;

        if (err)
        {
            throw std::runtime_error(error_message);
        }

        return fieldList;
    }

    void GribFile::getData(const std::string& field, const std::string& levelType, long level, double* data)
    {
        getDataOfKey(field, levelType, level, "values", data);
    }    

    void GribFile::getLatitudes(const std::string& field, const std::string& levelType, long level, double* data)
    {
        getDataOfKey(field, levelType, level, "latitudes", data);
    }

    void GribFile::getLongitudes(const std::string& field, const std::string& levelType, long level, double* data)
    {
        getDataOfKey(field, levelType, level, "longitudes", data);
    }

    void GribFile::loadFieldIndex()
    {
        int err;
        fieldNamesIndex = codes_index_new_from_file (context, const_cast<char*>(gribfile.c_str()), "shortName,typeOfLevel,level", &err);
        if (err)
        {
            throw std::runtime_error("index error "+ std::to_string(err));
        }
    }

    codes_handle* GribFile::selectHandle(const std::string& field, const std::string& levelType, long level)
    {
        int err;
        std::string error_message;
        codes_handle* handle = nullptr;

        if ((err=codes_index_select_string(fieldNamesIndex, "shortName", const_cast<char*>(field.c_str()))))
        {
            error_message = "codes_index_select_string error "+ std::to_string(err);
        }
        else if ((err=codes_index_select_string(fieldNamesIndex, "typeOfLevel", const_cast<char*>(levelType.c_str()))))
        {
            error_message = "codes_index_select_string error "+ std::to_string(err);
        }
        else if ((err=codes_index_select_long(fieldNamesIndex, "level", level)))
        {
            error_message = "codes_index_select_long error "+ std::to_string(err);
        }
        else 
        {
            handle = codes_handle_new_from_index(fieldNamesIndex, &err);

            if (err && err!=CODES_END_OF_INDEX)
            {
                error_message = "codes_handle_new_from_index error "+std::to_string(err);
            }
        }

        if (err)
        {
            throw std::runtime_error(error_message);
        }

        return handle;
    }

    void GribFile::getDataOfKey(const std::string& field, const std::string& levelType, long level, const std::string& key, double* data)
    {
        int err;
        std::string error_message;
        codes_handle* handle;

        handle = selectHandle(field, levelType, level);
        size_t length;
        if ((err = codes_get_size (handle, const_cast<char*>(key.c_str()), &length)))
        {
            error_message = "codes_get_length error "+std::to_string(err);
        }
        else
        {
            if ((err = codes_get_double_array(handle, const_cast<char*>(key.c_str()), data, &length)))
            {
                error_message = "codes_get_length error "+std::to_string(err);
            }
        }
        
        codes_handle_delete(handle);                            
        
        if (err)
        {
            throw std::runtime_error(error_message);
        }
    }
}