// ESPP_CLASS
#ifndef _FIXEDTUPLELIST_HPP
#define _FIXEDTUPLELIST_HPP

#include "log4espp.hpp"
//#include "types.hpp"

#include "Particle.hpp"
#include "esutil/ESPPIterator.hpp"
#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>
//#include "FixedListComm.hpp"

namespace espresso {
    class FixedTupleList: public TupleList  {
        protected:
            boost::signals2::connection con1, con2, con3, con4;
            shared_ptr<storage::Storage> storage;
            typedef std::vector<longint> tuple;
            typedef boost::unordered_map<longint, tuple> GlobalTuples;
            GlobalTuples globalTuples;
            using TupleList::add;

        public:
            FixedTupleList(shared_ptr<storage::Storage> _storage);
            ~FixedTupleList();

            void add(longint pid) { tmppids.push_back(pid); } // add particle id (called from python)
            void addTs() { addT(tmppids); tmppids.clear(); } // add tuple (called from python)
            void beforeSendParticles(ParticleList& pl, class OutBuffer& buf);
            void afterRecvParticles(ParticleList& pl, class InBuffer& buf);
            void onParticlesChanged();

            int getNumPart(longint pid); // get number of particles in globalmap for given pid

            // this signals the AT particles used to rebuild fixed pair, triple, quadruple AT bonds
            // for now FixedPairListAdress connects to it
            boost::signals2::signal2 <void, std::vector<longint>&, class OutBuffer&>
                   beforeSendATParticles;

            static void registerPython();

        private:
            tuple tmppids;
            bool addT(tuple pids); // add tuple
            static LOG4ESPP_DECL_LOGGER(theLogger);
    };
}

#endif

