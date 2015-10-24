/*
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
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

#ifndef _LBOUTPUT_SCREEN_HPP
#define _LBOUTPUT_SCREEN_HPP

#include "LBOutput.hpp"
#include "esutil/Timer.hpp"

namespace espressopp {
  namespace analysis {
		class LBOutputScreen : public LBOutput {
		public:
			LBOutputScreen(shared_ptr<System> _system,
                          shared_ptr< integrator::LatticeBoltzmann > _latticeboltzmann);

        /** Destructor for output. */
/*        ~LBOutputScreen ();
*/
			void writeOutput();
			void findLBMom();
			
			void setTimerOld(time_t _value);
			time_t getTimerOld();

			void setTimerNew(time_t _value);
			time_t getTimerNew();

			void setLBTimerOld(real _lbTime_old);
			real getLBTimerOld();
			
			void setLBTimerNew(real _lbTime_new);
			real getLBTimerNew();
			
			void setOldStepNum(long int _oldStepNum);
			long int getOldStepNum();
			
			static void registerPython();
			
		private:
			time_t timer_old, timer_new;
			real	lbTime_old, lbTime_new;
			esutil::WallTimer timeLBtoMD;  //!< used for timing
			int oldStepNum;
    };
  }
}

#endif
