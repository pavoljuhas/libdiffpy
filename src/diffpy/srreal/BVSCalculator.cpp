/*****************************************************************************
*
* libdiffpy         by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2010 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE_DANSE.txt for license information.
*
******************************************************************************
*
* class BVSCalculator -- concrete counter of pairs in a structure.
*
*****************************************************************************/

#include <cmath>
#include <cassert>

#include <diffpy/validators.hpp>
#include <diffpy/serialization.ipp>
#include <diffpy/srreal/AtomUtils.hpp>
#include <diffpy/srreal/BVSCalculator.hpp>

using namespace std;
using namespace diffpy::validators;

namespace diffpy {
namespace srreal {

// Local Helpers -------------------------------------------------------------

namespace {

int flattrilindex(int i, int j)
{
    // use lower-tringular part of a symmetric matrix where j <= i.
    if (i < j)  swap(i, j);
    int k = (i * (i + 1)) / 2 + j;
    return k;
}

}   // namespace

// Constructor ---------------------------------------------------------------

BVSCalculator::BVSCalculator()
{
    // default configuration
    const double valence_precision = 1e-5;
    // use very large rmax, it will be cropped by rmaxFromPrecision
    this->setRmax(100);
    BVParametersTablePtr bvtb(new BVParametersTable);
    this->setBVParamTable(bvtb);
    this->setValencePrecision(valence_precision);
    this->setStructure(mstructure);
    // attributes
    this->registerDoubleAttribute("valenceprecision", this,
            &BVSCalculator::getValencePrecision,
            &BVSCalculator::setValencePrecision);
    this->registerDoubleAttribute("rmaxused", this,
            &BVSCalculator::getRmaxUsed);
}

// Public Methods ------------------------------------------------------------

// results

QuantityType BVSCalculator::valences() const
{
    QuantityType rv(mstructure_cache.valences.begin(),
            mstructure_cache.valences.end());
    return rv;
}


QuantityType BVSCalculator::bvdiff() const
{
    QuantityType vobs = this->valences();
    assert(vobs.size() == this->value().size());
    int cntsites = this->countSites();
    QuantityType rv(cntsites);
    const QuantityType& vsim = this->value();
    for (int i = 0; i < cntsites; ++i)
    {
        rv[i] = fabs(vobs[i]) - fabs(vsim[i]);
    }
    return rv;
}


double BVSCalculator::bvmsdiff() const
{
    QuantityType bd = this->bvdiff();
    int cntsites = this->countSites();
    assert(int(bd.size()) == cntsites);
    double sumofsquares = 0.0;
    for (int i = 0; i < cntsites; ++i)
    {
        sumofsquares += mstructure->siteMultiplicity(i) *
            mstructure->siteOccupancy(i) * bd[i] * bd[i];
    }
    double totocc = mstructure->totalOccupancy();
    double rv = (totocc > 0.0) ? (sumofsquares / totocc) : 0.0;
    return rv;
}


double BVSCalculator::bvrmsdiff() const
{
    double rvsq = this->bvmsdiff();
    return sqrt(rvsq);
}


void BVSCalculator::setBVParamTable(BVParametersTablePtr bvtb)
{
    ensureNonNull("BVParametersTable", bvtb);
    mbvptable = bvtb;
}


BVParametersTablePtr& BVSCalculator::getBVParamTable()
{
    return mbvptable;
}


const BVParametersTablePtr& BVSCalculator::getBVParamTable() const
{
    return mbvptable;
}


void BVSCalculator::setValencePrecision(double eps)
{
    ensureEpsilonPositive("valenceprecision", eps);
    mvalenceprecision = eps;
}


double BVSCalculator::getValencePrecision() const
{
    return mvalenceprecision;
}


double BVSCalculator::getRmaxUsed() const
{
    double rv = min(this->getRmax(),
            this->rmaxFromPrecision(this->getValencePrecision()));
    return rv;
}

// Protected Methods ---------------------------------------------------------

// PairQuantity overloads

void BVSCalculator::resetValue()
{
    // structure data need to be cached for rmaxFromPrecision
    this->cacheStructureData();
    this->resizeValue(this->countSites());
    this->PairQuantity::resetValue();
}


void BVSCalculator::configureBondGenerator(BaseBondGenerator& bnds) const
{
    bnds.setRmax(this->getRmaxUsed());
}


void BVSCalculator::addPairContribution(const BaseBondGenerator& bnds,
        int summationscale)
{
    const int tp0 = mstructure_cache.typeofsite[bnds.site0()];
    const int tp1 = mstructure_cache.typeofsite[bnds.site1()];
    const int k = flattrilindex(tp0, tp1);
    assert(0 <= k && k < int(mstructure_cache.bpused.size()));
    const BVParam& bp = mstructure_cache.bpused[k];
    int v0 = mstructure_cache.valences[bnds.site0()];
    int v1 = mstructure_cache.valences[bnds.site1()];
    double valencehalf = bp.bondvalence(bnds.distance()) / 2.0;
    // do nothing if there are no bond parameters for this pair
    if (0 == valencehalf)    return;
    int pm0 = (v0 >= 0) ? 1 : -1;
    int pm1 = (v1 >= 0) ? 1 : -1;
    const double& o0 = mstructure->siteOccupancy(bnds.site0());
    const double& o1 = mstructure->siteOccupancy(bnds.site1());
    mvalue[bnds.site0()] += summationscale * pm0 * valencehalf * o1;
    mvalue[bnds.site1()] += summationscale * pm1 * valencehalf * o0;
}

// Private Methods -----------------------------------------------------------

void BVSCalculator::cacheStructureData()
{
    typedef boost::unordered_map<string,int> TypeToIndexMap;
    int cntsites = this->countSites();
    const BVParametersTable& bvtb = *(this->getBVParamTable());
    // index the unique atom types occuring in the structure
    mstructure_cache.valences.resize(cntsites);
    mstructure_cache.typeofsite.resize(cntsites);
    TypeToIndexMap atomtypeidx;
    for (int siteidx = 0; siteidx < cntsites; ++siteidx)
    {
        const string& smbl = mstructure->siteAtomType(siteidx);
        mstructure_cache.valences[siteidx] = bvtb.getAtomValence(smbl);
        TypeToIndexMap::iterator tpii;
        tpii = atomtypeidx.emplace(smbl, int(atomtypeidx.size())).first;
        mstructure_cache.typeofsite[siteidx] = tpii->second;
    }
    // build a flattened half-diagonal matrix of types vs bond parameters
    const int ntps = atomtypeidx.size();
    const int nbvpars = (ntps * (ntps + 1)) / 2;
    mstructure_cache.bpused.assign(nbvpars, BVParametersTable::none());
    TypeToIndexMap::const_iterator tpii, tpjj;
    for (tpii = atomtypeidx.begin(); tpii != atomtypeidx.end(); ++tpii)
    {
        for (tpjj = tpii; tpjj != atomtypeidx.end(); ++tpjj)
        {
            int k = flattrilindex(tpii->second, tpjj->second);
            assert(0 <= k && k < int(mstructure_cache.bpused.size()));
            mstructure_cache.bpused[k] = bvtb.lookup(tpii->first, tpjj->first);
        }
    }
}


double BVSCalculator::rmaxFromPrecision(double eps) const
{
    const BVParametersTable& bvtb = *(this->getBVParamTable());
    vector<BVParam>::const_iterator bpit = mstructure_cache.bpused.begin();
    double rv = 0.0;
    for (; bpit != mstructure_cache.bpused.end(); ++bpit)
    {
        // lookup parameters here again in case BVParametersTable changed
        // after building the cache of bpused.
        const BVParam& bp = bvtb.lookup(*bpit);
        rv = max(rv, bp.bondvalenceToDistance(eps));
    }
    return rv;
}

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::BVSCalculator)
BOOST_CLASS_EXPORT_IMPLEMENT(diffpy::srreal::BVSCalculator)

// End of file
