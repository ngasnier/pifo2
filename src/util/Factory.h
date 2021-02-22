#pragma once

#include <string>
#include <map>
#include <functional>

// Ref : https://stackoverflow.com/questions/5120768/how-to-implement-the-factory-method-pattern-in-c-correctly

/**
 * A generic factory class.
 */
template <class Base, class... Params> class Factory {
public:
    Factory() {}
    Factory(const Factory &) = delete;
    Factory &operator=(const Factory &) = delete;

    /**
     * Create the object and return it as a unique pointer.
     * 
     * @param key key name of the type of object to instanciate.
     * @param args list of additionnal args used to create the object.
     */
    auto create(const std::string& key, Params... args)
    {
        std::unique_ptr<Base> obj{creators[key](args...)};
        return obj;
    }

    /**
     * Register a type, creates a default creator function, and associates it with key name.
     * 
     * @tparam TDerived the class to register
     * @param name the key name of the class to register
     */
    template <typename TDerived>
    void registerType(const std::string& name)
    {
        static_assert(std::is_base_of<Base, TDerived>::value, "Factory::registerType doesn't accept this type because it doesn't derive from base class");
        creators[name] = std::bind(&createFunc<TDerived>);
    }

    /**
     * Register a creator function for a given key.
     * 
     * @param name the key name of the class to register
     * @param creator a creator function
     */
    void registerCreator(const std::string& name,
                        std::function<Base *(Params...)> &&creator) 
    {
        creators[name] = std::move(creator);
    }

protected:
    std::map<std::string, std::function<Base *(Params... )>> creators;

private:
    template <typename TDerived>
    static Base* createFunc(Params... args)
    {
        return new TDerived(args...);
    }
};
