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

// ESPP_CLASS
#ifndef _ANALYSIS_KINETICENERGY_HPP
#define _ANALYSIS_KINETICENERGY_HPP

#include "types.hpp"
#include "Observable.hpp"
#include "Temperature.hpp"
#include "NPart.hpp"
#include "interaction/Interaction.hpp"

namespace espressopp {
namespace analysis {

class KineticEnergy : public Observable {
 public:
  explicit KineticEnergy(shared_ptr<System> system) : Observable(system) {
    result_type = real_scalar;
    temperature_ = shared_ptr<Temperature>(new Temperature(system));
    npart_ = shared_ptr<NPart>(new NPart(system));
  }
  KineticEnergy(shared_ptr<System> system, shared_ptr<Temperature> temperature)
      : Observable(system) {
    result_type = real_scalar;
    temperature_ = temperature;
    npart_ = shared_ptr<NPart>(new NPart(system));
  }
  ~KineticEnergy() {}
  real compute_real() const;

  static void registerPython();
 private:
  shared_ptr<interaction::Interaction> interaction_;
  shared_ptr<Temperature> temperature_;
  shared_ptr<NPart> npart_;
  bool compute_at_;  // set to true then computeEnergyAA, otherwise computeEnergyCG
  bool compute_global_;
};
}  // end namespace analysis
}  // end namespace espressopp

#endif
