/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2011 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class AtomRadiiTable -- storage of empirical atomic radii
*
* $Id$
*
*****************************************************************************/

#include <algorithm>
#include <cctype>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <boost/serialization/export.hpp>

#include <diffpy/srreal/AtomRadiiTable.hpp>
#include <diffpy/serialization.hpp>

namespace diffpy {
namespace srreal {

using namespace std;

// Public Methods ------------------------------------------------------------

double AtomRadiiTable::lookup(const string& smbl) const
{
    boost::unordered_map<string,double>::const_iterator tb;
    tb = mcustomradius.find(smbl);
    if (tb != mcustomradius.end())  return tb->second;
    return this->tableLookup(smbl);
}


double AtomRadiiTable::tableLookup(const string& smbl) const
{
    return 0.0;
}


void AtomRadiiTable::setCustom(const string& smbl, double radius)
{
    mcustomradius[smbl] = radius;
}


void AtomRadiiTable::fromString(const string& s)
{
    // create deblanked string
    string s1;
    remove_copy_if(s.begin(), s.end(), back_inserter(s1), ::isspace);
    // replace commas with space so we can use the split function
    replace(s1.begin(), s1.end(), ',', ' ');
    istringstream ss1(s1);
    boost::unordered_map<string,double> rds;
    for (string w; ss1 >> w;)
    {
        string::size_type p = w.find(':');
        if (p == string::npos)
        {
            ostringstream emsg;
            emsg << "Invalid radius specification, missing ':' in '" <<
                w << "'.";
            throw invalid_argument(emsg.str());
        }
        string smbl = w.substr(0, p);
        double value;
        istringstream ssv(w.substr(p + 1));
        if (!(ssv >> value))
        {
            ostringstream emsg;
            emsg << "Invalid floating point number in '" << w << "'.";
            throw invalid_argument(emsg.str());
        }
        rds[smbl] = value;
    }
    // everything worked up to here, we can do the assignment
    boost::unordered_map<string,double>::const_iterator kv = rds.begin();
    for (; kv != rds.end(); ++kv)  mcustomradius[kv->first] = kv->second;
}


void AtomRadiiTable::resetCustom(const string& smbl)
{
    mcustomradius.erase(smbl);
}


void AtomRadiiTable::resetAll()
{
    mcustomradius.clear();
}


const boost::unordered_map<string,double>&
AtomRadiiTable::getAllCustom() const
{
    return mcustomradius;
}


string AtomRadiiTable::toString(string separator) const
{
    boost::unordered_map<string,double>::const_iterator tb;
    vector<string> symbols;
    for (tb = mcustomradius.begin(); tb != mcustomradius.end(); ++tb)
    {
        symbols.push_back(tb->first);
    }
    sort(symbols.begin(), symbols.end());
    ostringstream rv;
    vector<string>::const_iterator smbi;
    for (smbi = symbols.begin(); smbi != symbols.end(); ++smbi)
    {
        if (smbi != symbols.begin())  rv << separator;
        rv << *smbi << ':' << mcustomradius.at(*smbi);
    }
    return rv.str();
}

}   // srreal
}   // diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZE(diffpy::srreal::AtomRadiiTable)
BOOST_CLASS_EXPORT(diffpy::srreal::AtomRadiiTable)

// End of file
