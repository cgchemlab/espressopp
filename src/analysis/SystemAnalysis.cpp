/*
  Copyright (c) 2015
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

#include <string>
#include <utility>
#include <vector>
#include "python.hpp"
#include "SystemAnalysis.hpp"
#include "integrator/MDIntegrator.hpp"

namespace espressopp {
namespace analysis {

void SystemAnalysis::perform_action() {
  current_step_ = integrator_->getStep();
  compute_potential_energy();
  compute_kinetic_energy();
  write();
}

void SystemAnalysis::compute_potential_energy() {
  pot_energy_.clear();
  real e = 0.0;
  epot_ = 0.0;
  for (int i = 0; i < system_->getNumberOfInteractions(); i++) {
    e = system_->getInteraction(i)->computeEnergy();
    epot_ += e;
    pot_energy_.push_back(e);
  }
}

void SystemAnalysis::compute_kinetic_energy() {
  T_ = temp_->computeRaw();
  ekin_ = (3.0/2.0) * npart_->compute_real() * T_;
}

void SystemAnalysis::write() {
  // That's a text file, let's only write on one node.
  if (system_->comm->rank() == 0) {
    std::ofstream output_file;
    if (!header_written_) {
      output_file.open(file_name_.c_str(), std::fstream::out);
      last_line_ = prepare_line();
      header_written_ = true;
    } else {
      output_file.open(file_name_.c_str(), std::ofstream::out | std::ofstream::app);
      last_line_ = prepare_line();
    }
    output_file << last_line_ << std::endl;
    output_file.close();
  }
}

std::string SystemAnalysis::prepare_header() {
  std::stringstream line;
  int size_observables = observables_.size();
  line << "#step" << delimiter_
              << "T" << delimiter_ << "Ek" << delimiter_ << "Epot"
              << delimiter_;
  // Write columns for interactions.
  for (int i = 0; i < system_->getNumberOfInteractions(); i++)
    line << "e" << i << delimiter_;

  // Write columns for additional observables.
  for (int i = 0; i < size_observables - 1; i++) {
    line << observables_[i].first << delimiter_;
  }
  if (size_observables - 1 >= 0)
    line << observables_[size_observables-1].first;

  // Close wile.
  line << "\n";

  return line.str();
}

std::string SystemAnalysis::prepare_line() {
  if (system_->comm->rank() == 0) {
    std::stringstream line;
    int size_observables = observables_.size();
    if (!header_written_) {  // First run;
      line << prepare_header();
      header_written_ = true;
    }

    line << current_step_ << delimiter_
      << T_ << delimiter_
      << ekin_ << delimiter_
      << epot_ << delimiter_;
    int number_of_interactions = pot_energy_.size();
    for (int i = 0; i < number_of_interactions; i++) {
      line << pot_energy_[i];
      if (i != number_of_interactions - 1 || size_observables > 0)
        line << delimiter_;
    }
    // Add columns from observable.
    for (int i = 0; i < size_observables - 1; i++) {
      line << observables_[i].second->compute_real() << delimiter_;
    }
    if (size_observables - 1 >= 0)
      line << observables_[size_observables-1].second->compute_real();

    // Close file.
    return line.str();
  }
}

void SystemAnalysis::info() {
  if (system_->comm->rank() == 0) {
    if (!show_header_) {
      std::cout << prepare_header();
      show_header_ = true;
    }
    std::cout << last_line_ << std::endl;
  }
}

void SystemAnalysis::add_observable(std::string name, shared_ptr<Observable> obs) {
  std::cout << "adding observable " << name;
  observables_.push_back(std::make_pair(name, obs));
}

void SystemAnalysis::registerPython() {
  using namespace espressopp::python;  // NOLINT
  class_<SystemAnalysis, bases< ParticleAccess > >
      ("analysis_SystemAnalysis", init<
          shared_ptr<System>,
          shared_ptr<integrator::MDIntegrator>,
          std::string,
          std::string
          >())
      .def("add_observable", &SystemAnalysis::add_observable)
      .def("info", &SystemAnalysis::info);
}

}  // end namespace analysis
}  // end namespace espressopp
