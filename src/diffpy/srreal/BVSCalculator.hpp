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
* class BVSCalculator -- bond valence sums calculator
*
*****************************************************************************/

#ifndef BVSCALCULATOR_HPP_INCLUDED
#define BVSCALCULATOR_HPP_INCLUDED

#include <diffpy/srreal/PairQuantity.hpp>
#include <diffpy/srreal/BVParametersTable.hpp>

namespace diffpy {
namespace srreal {

class BVSCalculator : public PairQuantity
{
    public:

        // constructor
        BVSCalculator();

        // results
        /// expected valence per each site
        QuantityType valences() const;
        /// difference between expected and calculated absolute valence at
        /// each site.  Positive for underbonding, negative for overbonding.
        QuantityType bvdiff() const;
        /// mean square difference of BVS from the expected values
        double bvmsdiff() const;
        /// root mean square difference of BVS from the expected values
        double bvrmsdiff() const;

        // access and configuration of BVS parameters
        void setBVParamTable(BVParametersTablePtr);
        BVParametersTablePtr& getBVParamTable();
        const BVParametersTablePtr& getBVParamTable() const;

        // R-range configuration using the valence precision
        /// set cutoff value for bond valence contributions
        void setValencePrecision(double);
        /// return cutoff value for bond valence contributions
        double getValencePrecision() const;
        /// effective rmax value, where valence contributions become smaller
        /// than the valence precision cutoff.  Always less or equal to rmax.
        double getRmaxUsed() const;

    protected:

        // PairQuantity overloads
        virtual void resetValue();
        virtual void configureBondGenerator(BaseBondGenerator&) const;
        virtual void addPairContribution(const BaseBondGenerator&, int);

    private:

        // methods
        void cacheStructureData();
        /// rmax necessary for achieving the specified valence precision
        double rmaxFromPrecision(double) const;

        // data
        // configuration
        BVParametersTablePtr mbvptable;
        double mvalenceprecision;
        // cache
        struct {
            // TODO unused baresymbols, consider removing in libdiffpy-2.0.
            std::vector<std::string> baresymbols;
            std::vector<int> valences;
            std::vector<int> typeofsite;
            std::vector<BVParam> bpused;
        } mstructure_cache;

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            ar & boost::serialization::base_object<PairQuantity>(*this);
            ar & mbvptable;
            ar & mvalenceprecision;
            ar & mstructure_cache.baresymbols;
            ar & mstructure_cache.valences;
            if (version >= 1) {
                ar & mstructure_cache.typeofsite;
                ar & mstructure_cache.bpused;
            }
        }

};  // class BVSCalculator

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_CLASS_VERSION(diffpy::srreal::BVSCalculator, 1)
BOOST_CLASS_EXPORT_KEY(diffpy::srreal::BVSCalculator)

#endif  // BVSCALCULATOR_HPP_INCLUDED
