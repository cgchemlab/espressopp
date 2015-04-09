/*
  Copyright (C) 2015
      Jakub Krajniak (jkrajniak at gmail.com)

  This file is part of ESPResSo++.

  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _FIXEDPAIRLISTTYPE_HPP
#define _FIXEDPAIRLISTTYPE_HPP

#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>

#include "log4espp.hpp"

#include "FixedPairList.hpp"
#include "esutil/ESPPIterator.hpp"
#include "types.hpp"

namespace espressopp {
/**
 * This is a subclass of FixedPairList. This fixed tuple list
 * only holds the bonds of defined types. Whenever particles change
 * type then the bonds are recalculated. This is done locally.
 */
class FixedPairListType : public FixedPairList {
 public:
  FixedPairListType(shared_ptr<storage::Storage> _storage, longint type1, longint type2);

  bool add(longint pid1, longint pid2);
  void onParticlesChanged();

  static void registerPython();

 private:
  longint type1_, type2_;

  inline bool check_type(longint p1_type, longint p2_type) {
    return ((type1_ == p1_type && type2_ == p2_type) || (type1_ == p2_type && type2_ == p1_type));
  }

  using PairList::add;
  static LOG4ESPP_DECL_LOGGER(theLogger);
};
}  // end namespace espressopp
#endif

