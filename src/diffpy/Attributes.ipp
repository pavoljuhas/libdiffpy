/*****************************************************************************
*
* libdiffpy         by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE_DANSE.txt for license information.
*
******************************************************************************
*
* class Attributes - interface for calling setter and getter methods using
*   their string names.
*
*****************************************************************************/

#ifndef ATTRIBUTES_IPP_INCLUDED
#define ATTRIBUTES_IPP_INCLUDED

#include <stdexcept>
#include <cassert>

namespace diffpy {
namespace attributes {

//////////////////////////////////////////////////////////////////////////////
// class DoubleAttribute
//////////////////////////////////////////////////////////////////////////////

template <class T>
class GetSetTypes
{
    public:

        // getter types
        typedef double(*GF)(const T*);
        typedef double(T::*GMF)() const;
        typedef const double&(T::*GMFR)() const;

        // setter types
        typedef void(*SF)(T*, const double&);
        typedef void(T::*SMF)(double);
        typedef void(T::*SMFR)(const double&);

};


class GetSetTools
{
    public:

        template <class T>
        static double
        callgetter(const T* obj, typename GetSetTypes<T>::GF getter)
        {
            return getter(obj);
        }

        template <class T>
        static double
        callgetter(const T* obj, typename GetSetTypes<T>::GMF getter)
        {
            return (obj->*getter)();
        }

        template <class T>
        static double
        callgetter(const T* obj, typename GetSetTypes<T>::GMFR getter)
        {
            return (obj->*getter)();
        }

        template <class T>
        static void
        callsetter(T* obj, typename GetSetTypes<T>::SF setter, const double& x)
        {
            setter(obj, x);
        }

        template <class T>
        static void
        callsetter(T* obj, typename GetSetTypes<T>::SMF setter, const double& x)
        {
            (obj->*setter)(x);
        }

        template <class T>
        static void
        callsetter(T* obj, typename GetSetTypes<T>::SMFR setter, const double& x)
        {
            (obj->*setter)(x);
        }

};


template <class T, class Getter, class Setter>
class DoubleAttribute : public BaseDoubleAttribute
{
    public:

        // constructor
        DoubleAttribute(Getter g, Setter s)
        {
            mgetter = g;
            msetter = s;
        }

        double getValue(const Attributes* obj) const
        {
            const T* tobj = dynamic_cast<const T*>(obj);
            assert(tobj);
            double rv = GetSetTools::callgetter(tobj, mgetter);
            return rv;
        }

        void setValue(Attributes* obj, double value)
        {
            if (this->isreadonly())  throwDoubleAttributeReadOnly();
            T* tobj = dynamic_cast<T*>(obj);
            assert(tobj);
            GetSetTools::callsetter(tobj, msetter, value);
        }

        bool isreadonly() const
        {
            return !msetter;
        }

    private:
        // data
        Getter mgetter;
        Setter msetter;

};

//////////////////////////////////////////////////////////////////////////////
// class Attributes
//////////////////////////////////////////////////////////////////////////////

// Template Protected Methods ------------------------------------------------

template <class T, class Getter>
void Attributes::registerDoubleAttribute(
        const std::string& name, T* obj, Getter g)
{
    typedef  void (T::*Setter)(double x);
    registerBaseDoubleAttribute(this, name,
            new DoubleAttribute<T, Getter, Setter>(g, NULL));
}


template <class T, class Getter, class Setter>
void Attributes::registerDoubleAttribute(
        const std::string& name, T* obj, Getter g, Setter s)
{
    registerBaseDoubleAttribute(this, name,
            new DoubleAttribute<T, Getter, Setter>(g, s));
}

}   // namespace attributes
}   // namespace diffpy

// vim:ft=cpp:

#endif  // ATTRIBUTES_IPP_INCLUDED
