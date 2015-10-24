/*
  Copyright (C) 2012-2015 Max Planck Institute for Polymer Research
  Copyright (C) 2008-2011 Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
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

#include "python.hpp"
#include "LatticeBoltzmann.hpp"
#include <iomanip>
#include <iostream>
#include <fstream>

#include "boost/serialization/vector.hpp"
#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "esutil/RNG.hpp"
#include "esutil/Grid.hpp"
#include "bc/BC.hpp"
#include "mpi.hpp"

#define REQ_HALO_SPREAD 501
#define COMM_DIR_0 700
#define COMM_DIR_1 701
#define COMM_DIR_2 702
#define COMM_DIR_3 703
#define COMM_DIR_4 704
#define COMM_DIR_5 705
#define COMM_FORCE_0 706
#define COMM_FORCE_1 707
#define COMM_FORCE_2 708
#define COMM_FORCE_3 709
#define COMM_FORCE_4 710
#define COMM_FORCE_5 711
#define COMM_DEN_0 712
#define COMM_DEN_1 713
#define COMM_DEN_2 714
#define COMM_DEN_3 715
#define COMM_DEN_4 716
#define COMM_DEN_5 717

using namespace boost;

namespace espressopp {

	using namespace iterator;
	namespace integrator {
		LOG4ESPP_LOGGER(LatticeBoltzmann::theLogger, "LatticeBoltzmann");

		/* LB Constructor; expects 3 reals, 1 vector and 5 integers */
		LatticeBoltzmann::LatticeBoltzmann(shared_ptr<System> _system,
																			 Int3D _nodeGrid, Int3D _Ni,
																			 real _a, real _tau,
																			 int _numDims, int _numVels)
		: Extension(_system), nodeGrid(_nodeGrid), Ni(_Ni), a(_a), tau(_tau), numDims(_numDims), numVels(_numVels)
		{
			/* create storage for variables that are static (not changing on the run) and live on a lattice node */
			c_i	= std::vector<Real3D>(_numVels, Real3D(0.,0.,0.));
			eqWeight = std::vector<real>(_numVels, 0.);
			inv_b = std::vector<real>(_numVels, 0.);
			phi = std::vector<real>(_numVels, 0.);
			myNeighbour = std::vector<int>(2*_numDims, 0);

			/* set flags for extended simulations to their defaults */
			setLBTempFlag(0);																// no fluctuations

			/* setup simulation parameters */
			setStart(0);																			// set coupling start flag to 0
			setStepNum(0);																		// set step number to 0

			/* setup random numbers generator */
			if (!_system->rng) {
				throw std::runtime_error("system has no RNG");
			}
			rng = _system->rng;

			/* setup default external (and coupling) forces parameters */
			setExtForceFlag(0);															// no external forces
			setCouplForceFlag(0);														// no LB to MD coupling
			setNSteps(1);																		// # MD steps between LB update
			
			int _Npart = _system->storage->getNRealParticles();
			int _totNPart = 0;
			mpi::all_reduce(*getSystem()->comm, _Npart, _totNPart, std::plus<int>());
			setTotNPart(_totNPart);
			
			fOnPart = std::vector<Real3D>(_totNPart+1,Real3D(0.));	// +1 as id starts with 1

			/* if coupling is present initialise related flags, coefficients and arrays */
			if (_totNPart != 0) {
				setCouplForceFlag(1);													// make LB to MD coupling
				setFricCoeff(5.);															// friction coeffitient
				
				for(int _id = 0; _id <= _totNPart; _id++) {
					setFOnPart(_id, Real3D(0.));
				}
			}

			/* setup domain decompositions for LB */
			setHaloSkin(1);
			findMyNeighbours();
			
			Int3D _numSites = Int3D(0,0,0);
			Real3D _myLeft = Real3D(0.,0.,0.);
			real _L = 0.;
			
			longint _myRank = getSystem()->comm->rank();
			Int3D _myPosition = getMyPosition();

			for (int _dim = 0; _dim < 3; ++_dim) {
				_L = getSystem()->bc->getBoxL().getItem(_dim);
				_numSites[_dim] = floor((_myPosition[_dim]+1)*_L/(_nodeGrid[_dim]*getA())) -
													floor(_myPosition[_dim]*_L/(_nodeGrid[_dim]*getA())) +
													2 * getHaloSkin();
				_myLeft[_dim] = floor(_myPosition[_dim]*_L/_nodeGrid[_dim]);
			}
#warning: probably one needs to eliminate Ni as input parameter and calculate it through Li and a
//			setNi(_L / getA());
			setMyNi (_numSites);
			setMyLeft (_myLeft);
			
			/* initialise general lattice parameters on a site */
			LatticePar(_system, getNumVels(),getA(),getTau());
//			LatticePar::initEqWeights();
			
			/* stretch lattices resizing them in 3 dimensions */
			lbfluid = new lblattice;
			ghostlat = new lblattice;
			lbmom = new lbmoments;
			lbfor = new lbforces;

			(*lbfluid).resize(_numSites[0]);
			(*ghostlat).resize(_numSites[0]);
			(*lbmom).resize(_numSites[0]);
			(*lbfor).resize(_numSites[0]);
			
			for (int i = 0; i < _numSites[0]; i++) {
				(*lbfluid)[i].resize(_numSites[1]);
				(*ghostlat)[i].resize(_numSites[1]);
				(*lbmom)[i].resize(_numSites[1]);
				(*lbfor)[i].resize(_numSites[1]);
				for (int j = 0; j < _numSites[1]; j++) {
					(*lbfluid)[i][j].resize(_numSites[2]);
					(*ghostlat)[i][j].resize(_numSites[2]);
					(*lbmom)[i][j].resize(_numSites[2]);
					(*lbfor)[i][j].resize(_numSites[2]);
				}
			}

			/* initialise global weights and coefficients from the local ones */
			initLatticeModel();

			// reset timers
			colstream.reset();
			comm.reset();
			swapping.reset();
			time_colstr = 0.;
			time_comm = 0.;
			time_sw = 0.;
		}
		
/*******************************************************************************************/
		
		void LatticeBoltzmann::disconnect() {
			_recalc2.disconnect();
			_befIntV.disconnect();
			
			delete (lbfluid);
			delete (ghostlat);
			delete (lbmom);
			delete (lbfor);
		}

		void LatticeBoltzmann::connect() {
#warning: need to correct zeroMDCMVel if we not at the VERY start of the simulation. It zeros CM Vel EVERY time when integrator.run(X_steps) starts!!!
			_recalc2 = integrator->recalc2.connect ( boost::bind(&LatticeBoltzmann::zeroMDCMVel, this));
			_befIntV = integrator->befIntV.connect ( boost::bind(&LatticeBoltzmann::makeLBStep, this));
		}
		
/*******************************************************************************************/
		
		/* Setter and getter for the parallelisation things */
		void LatticeBoltzmann::setMyNeighbour (int _dir, int _rank) { myNeighbour[_dir] = _rank;}
		int LatticeBoltzmann::getMyNeighbour (int _dir) { return myNeighbour[_dir];}

		void LatticeBoltzmann::setMyPosition (Int3D _myPosition) { myPosition = _myPosition;}
		Int3D LatticeBoltzmann::getMyPosition () { return myPosition;}
		
		void LatticeBoltzmann::setNodeGrid (Int3D _nodeGrid) { nodeGrid = _nodeGrid;}
		Int3D LatticeBoltzmann::getNodeGrid () {return nodeGrid;}
		
		void LatticeBoltzmann::setHaloSkin (int _haloSkin) { haloSkin = _haloSkin;}
		int LatticeBoltzmann::getHaloSkin () {return haloSkin;}
	
		void LatticeBoltzmann::setMyNi (Int3D _myNi) { myNi = _myNi;}
		Int3D LatticeBoltzmann::getMyNi () {return myNi;}
		
		void LatticeBoltzmann::setMyLeft (Real3D _myLeft) { myLeft = _myLeft;}
		Real3D LatticeBoltzmann::getMyLeft () {return myLeft;}
		
		/* Setter and getter for the lattice model */
		void LatticeBoltzmann::setNi (Int3D _Ni) { Ni = _Ni;}
		Int3D LatticeBoltzmann::getNi () {return Ni;}

		void LatticeBoltzmann::setA (real _a) { a = _a;
			printf ("Lattice spacing %4.2f\n", a);}
		real LatticeBoltzmann::getA () { return a;}

		void LatticeBoltzmann::setTau (real _tau) { tau = _tau;
			printf ("lattice time step %4.2f\n", tau);}
		real LatticeBoltzmann::getTau () { return tau;}

		void LatticeBoltzmann::setNumVels (int _numVels) { numVels = _numVels;
			printf ("Number of Velocities %2d; ", numVels);}
		int LatticeBoltzmann::getNumVels () { return numVels;}

		void LatticeBoltzmann::setNumDims (int _numDims) { numDims = _numDims;
			printf ("Number of Dimensions %2d; ", numDims);}
		int LatticeBoltzmann::getNumDims () { return numDims;}
		
		void LatticeBoltzmann::setEqWeight (int _l, real _value) { eqWeight[_l] = _value;}
		real LatticeBoltzmann::getEqWeight (int _l) {return eqWeight[_l];}
		
		void LatticeBoltzmann::setCi (int _l, Real3D _vec) {c_i[_l] = _vec;}
		Real3D LatticeBoltzmann::getCi (int _l) {return c_i[_l];}
		
		void LatticeBoltzmann::setCs2 (real _cs2) { cs2 = _cs2;}
		real LatticeBoltzmann::getCs2 () { return cs2;}
		
		void LatticeBoltzmann::setInvB (int _l, real _value) {inv_b[_l] = _value;}
		real LatticeBoltzmann::getInvB (int _l) {return inv_b[_l];}

		/* Setter and getter for simulation parameters */
		void LatticeBoltzmann::setStepNum (int _step) { stepNum = _step;}
		int LatticeBoltzmann::getStepNum () { return stepNum;}

		void LatticeBoltzmann::setCopyTimestep (real _copyTimestep) { copyTimestep = _copyTimestep;}
		real LatticeBoltzmann::getCopyTimestep () { return copyTimestep;}
		
		void LatticeBoltzmann::setStart (int _start) { start = _start;}
		int LatticeBoltzmann::getStart () { return start;}

		/* Setter and getter for LB temperature things */
		void LatticeBoltzmann::setLBTemp (real _lbTemp) { lbTemp = _lbTemp; initFluctuations();}
		real LatticeBoltzmann::getLBTemp () { return lbTemp;}

		void LatticeBoltzmann::setLBTempFlag (int _lbTempFlag) {lbTempFlag = _lbTempFlag;}
		int LatticeBoltzmann::getLBTempFlag () {return lbTempFlag;}

		void LatticeBoltzmann::setPhi (int _l, real _value) {phi[_l] = _value;}
		real LatticeBoltzmann::getPhi (int _l) {return phi[_l];}

		/* Setter and getting for LB viscosity control */
		void LatticeBoltzmann::setViscB (real _visc_b) {visc_b = _visc_b;
			setGammaB((getNumDims()*_visc_b-getCs2()*getTau()*convTimeMDtoLB()*convLenMDtoLB()/convMassMDtoLB())/(getNumDims()*_visc_b+getCs2()*getTau()*convTimeMDtoLB()*convLenMDtoLB()/convMassMDtoLB()));}
		real LatticeBoltzmann::getViscB () { return visc_b;}
		
		void LatticeBoltzmann::setViscS (real _visc_s) {visc_s = _visc_s;
			setGammaS((2.*_visc_s-getCs2()*getTau()*convTimeMDtoLB()*convLenMDtoLB()/convMassMDtoLB())/(2.*_visc_s+getCs2()*getTau()*convTimeMDtoLB()*convLenMDtoLB()/convMassMDtoLB()));}
		real LatticeBoltzmann::getViscS () { return visc_s;}
		
		void LatticeBoltzmann::setGammaB (real _gamma_b) {gamma_b = _gamma_b; initGammas(0);}
		real LatticeBoltzmann::getGammaB () { return gamma_b;}

		void LatticeBoltzmann::setGammaS (real _gamma_s) {gamma_s = _gamma_s; initGammas(1);}
		real LatticeBoltzmann::getGammaS () { return gamma_s;}

		void LatticeBoltzmann::setGammaOdd (real _gamma_odd) {gamma_odd = _gamma_odd; initGammas(2);}
		real LatticeBoltzmann::getGammaOdd () { return gamma_odd;}

		void LatticeBoltzmann::setGammaEven (real _gamma_even) {gamma_even = _gamma_even; initGammas(3);}
		real LatticeBoltzmann::getGammaEven () { return gamma_even;}

		/* Setter and getter for external and coupling force control */
		void LatticeBoltzmann::setExtForceFlag (int _extForceFlag) {extForceFlag = _extForceFlag;}
		int LatticeBoltzmann::getExtForceFlag () {return extForceFlag;}
		
		void LatticeBoltzmann::setCouplForceFlag (int _couplForceFlag) {couplForceFlag = _couplForceFlag;}
		int LatticeBoltzmann::getCouplForceFlag () {return couplForceFlag;}

		void LatticeBoltzmann::setExtForceLoc (Int3D _Ni, Real3D _extForceLoc) {
			return (*lbfor)[_Ni.getItem(0)][_Ni.getItem(1)][_Ni.getItem(2)].setExtForceLoc(_extForceLoc);	}
		Real3D LatticeBoltzmann::getExtForceLoc (Int3D _Ni) {
			return (*lbfor)[_Ni.getItem(0)][_Ni.getItem(1)][_Ni.getItem(2)].getExtForceLoc();	}
		void LatticeBoltzmann::addExtForceLoc (Int3D _Ni, Real3D _extForceLoc) {
			return (*lbfor)[_Ni.getItem(0)][_Ni.getItem(1)][_Ni.getItem(2)].addExtForceLoc(_extForceLoc);	}

		void LatticeBoltzmann::setFricCoeff (real _fricCoeff) { fricCoeff = _fricCoeff;}
		real LatticeBoltzmann::getFricCoeff () { return fricCoeff;}
		
		void LatticeBoltzmann::setNSteps (int _nSteps) { nSteps = _nSteps;}
		int LatticeBoltzmann::getNSteps () { return nSteps;}

		void LatticeBoltzmann::setTotNPart (int _totNPart) { totNPart = _totNPart;}
		int LatticeBoltzmann::getTotNPart () { return totNPart;}
		
		void LatticeBoltzmann::setFOnPart (int _id, Real3D _fOnPart) {fOnPart[_id] = _fOnPart;}
		Real3D LatticeBoltzmann::getFOnPart (int _id) {return fOnPart[_id];}
		void LatticeBoltzmann::addFOnPart (int _id, Real3D _fOnPart) {fOnPart[_id] += _fOnPart;}
		
		/* Setter and getter for access to population values */
		void LatticeBoltzmann::setLBFluid (Int3D _Ni, int _l, real _value) {
			(*lbfluid)[_Ni.getItem(0)][_Ni.getItem(1)][_Ni.getItem(2)].setF_i(_l, _value);	}
		real LatticeBoltzmann::getLBFluid (Int3D _Ni, int _l) {
			return (*lbfluid)[_Ni.getItem(0)][_Ni.getItem(1)][_Ni.getItem(2)].getF_i(_l);	}
		
		void LatticeBoltzmann::setGhostFluid (Int3D _Ni, int _l, real _value) {
			(*ghostlat)[_Ni.getItem(0)][_Ni.getItem(1)][_Ni.getItem(2)].setF_i(_l, _value);	}
	
		void LatticeBoltzmann::setLBMom (Int3D _Ni, int _l, real _value) {
			(*lbmom)[_Ni.getItem(0)][_Ni.getItem(1)][_Ni.getItem(2)].setMom_i(_l, _value);	}
		real LatticeBoltzmann::getLBMom (Int3D _Ni, int _l) {
			return (*lbmom)[_Ni.getItem(0)][_Ni.getItem(1)][_Ni.getItem(2)].getMom_i(_l);	}
		
		/* Helpers for MD to LB (and vice versa) unit conversion */
		real LatticeBoltzmann::convMassMDtoLB() {return 1.;}
#warning: need a foolproof in case there is no access to the integrator (and getTimeStep) yet
		real LatticeBoltzmann::convTimeMDtoLB() {return 1. / (integrator->getTimeStep() * getTau());}
		real LatticeBoltzmann::convLenMDtoLB() {
			return getNi().getItem(0) / (getSystem()->bc->getBoxL().getItem(0) * getA());}
		
		/* Profiling definitions */
		void LatticeBoltzmann::setProfStep (int _profStep) { profStep = _profStep;}
		int LatticeBoltzmann::getProfStep () { return profStep;}
		
/*******************************************************************************************/

		/* INITIALIZATION OF THE LATTICE MODEL: EQ.WEIGHTS, CI'S, ... */
		void LatticeBoltzmann::initLatticeModel () {
			using std::setprecision;
			using std::fixed;
			using std::setw;

			setCs2(1. / 3. * getA() * getA() / (getTau() * getTau()));

			std::cout << setprecision(4); std::cout << fixed;

			setCi( 0, Real3D(0.,  0.,  0.));
			setCi( 1, Real3D(1.,  0.,  0.)); setCi( 2, Real3D(-1.,  0.,  0.));
			setCi( 3, Real3D(0.,  1.,  0.)); setCi( 4, Real3D( 0., -1.,  0.));
			setCi( 5, Real3D(0.,  0.,  1.)); setCi( 6, Real3D( 0.,  0., -1.));
			setCi( 7, Real3D(1.,  1.,  0.)); setCi( 8, Real3D(-1., -1.,  0.));
			setCi( 9, Real3D(1., -1.,  0.)); setCi(10, Real3D(-1.,  1.,  0.));
			setCi(11, Real3D(1.,  0.,  1.)); setCi(12, Real3D(-1.,  0., -1.));
			setCi(13, Real3D(1.,  0., -1.)); setCi(14, Real3D(-1.,  0.,  1.));
			setCi(15, Real3D(0.,  1.,  1.)); setCi(16, Real3D( 0., -1., -1.));
			setCi(17, Real3D(0.,  1., -1.)); setCi(18, Real3D( 0., -1.,  1.));
			
/*			Int3D _myNi = getMyNi();
			int _numVels = getNumVels();
			for (int i = 0; i < _myNi[0]; i++) {
				for (int j = 0; j < _myNi[1]; j++) {
					for (int k = 0; k < _myNi[2]; k++) {
						for (int l = 0; l < _numVels; l++) {
							(*lbfluid)[i][j][k].initLatticeModelLoc();
						}
					}
				}
			}
*/
			longint _myRank = getSystem()->comm->rank();
			if (_myRank == 0) {
				std::cout << "-------------------------------------" << std::endl;
				std::cout << "Lattice Boltzmann lattice model D3Q19" << std::endl;
				std::cout << "-------------------------------------" << std::endl;
				std::cout << "Lattice of " << getNi().getItem(0);
				std::cout << " x " << getNi().getItem(1);
				std::cout << " x " << getNi().getItem(2) << " size " << std::endl;
				std::cout << "-------------------------------------" << std::endl;
			}
			
			// pass local eq. weights and inversed coeff. to the global ones
			for (int l = 0; l < getNumVels(); l++) {
				setEqWeight(l, LatticePar::getEqWeightLoc(l));
				setInvB(l, LatticePar::getInvBLoc(l));
			}
		}

/*******************************************************************************************/
		
    /* (RE)INITIALIZATION OF GAMMAS */
    void LatticeBoltzmann::initGammas (int _idGamma) {
      using std::setprecision;
      using std::fixed;
      using std::setw;

      // (re)set values of gammas depending on the id of the gamma that was changed
      for (int i = 0; i < getMyNi().getItem(0); i++) {
        for (int j = 0; j < getMyNi().getItem(1); j++) {
          for (int k = 0; k < getMyNi().getItem(2); k++) {
            for (int l = 0; l < getNumVels(); l++) {		// ich glaube die Schleife ist überflüssig hier
              if (_idGamma == 0) (*lbfluid)[i][j][k].setGammaBLoc(getGammaB());
              if (_idGamma == 1) (*lbfluid)[i][j][k].setGammaSLoc(getGammaS());
              if (_idGamma == 2) (*lbfluid)[i][j][k].setGammaOddLoc(getGammaOdd());
              if (_idGamma == 3) (*lbfluid)[i][j][k].setGammaEvenLoc(getGammaEven());
            }
          }
        }
      }
			
      // print for control
			longint _myRank = getSystem()->comm->rank();
			if (_myRank == 0) {
				std::cout << setprecision(8);
				std::cout << "One of the gamma's controlling viscosities has been changed:\n";
				if (_idGamma == 0) std::cout << "  gammaB is " << (*lbfluid)[0][0][0].getGammaBLoc() << "\n";
				if (_idGamma == 1) std::cout << "  gammaS is " << (*lbfluid)[0][0][0].getGammaSLoc() << "\n";
				if (_idGamma == 2) std::cout << ", gammaOdd is " << (*lbfluid)[0][0][0].getGammaOddLoc() << "\n";
				if (_idGamma == 3) std::cout << ", gammaEven is " << (*lbfluid)[0][0][0].getGammaEvenLoc() << "\n";
				std::cout << "-------------------------------------\n";
			}
    }
		
/*******************************************************************************************/
		
    /* (RE)INITIALIZATION OF THERMAL FLUCTUATIONS */
    void LatticeBoltzmann::initFluctuations () {
      using std::setprecision;
      using std::fixed;
      using std::setw;

      /* set amplitudes of local fluctuations */
      real _lbTemp;
      real mu, a3;

			longint _myRank = getSystem()->comm->rank();

			if (_myRank == 0) {
				std::cout << "Mass conversion coeff. MD->LB: " << convMassMDtoLB() << "\n";
				std::cout << "Length conversion coeff. MD->LB: " << convLenMDtoLB() << "\n";
				std::cout << "Time conversion coeff. MD->LB: " << convTimeMDtoLB() << "\n";
				std::cout << "-------------------------------------\n";
			}
			
			_lbTemp = getLBTemp() * convMassMDtoLB() * pow(convLenMDtoLB() / convTimeMDtoLB(), 2.);
      a3 = getA() * getA() * getA();    // a^3
      mu = _lbTemp / (getCs2() * a3);   // thermal mass density

			if (_lbTemp == 0.) {
        // account for fluctuations being turned off
        setLBTempFlag(0);
				if (_myRank == 0) {
					std::cout << "The temperature of the LB-fluid is 0. The fluctuations are turned off!\n";
				}
      } else {
        // account for fluctuations being turned on!
        setLBTempFlag(1);

				if (_myRank == 0) {
					std::cout << setprecision(8);
					std::cout << fixed;   // some output tricks
					std::cout << "The fluctuations have been introduced into the system:\n";
					std::cout << "lbTemp = " << getLBTemp() << "\n";
				}
				
        setPhi(0, 0.);
        setPhi(1, 0.); setPhi(2, 0.); setPhi(3, 0.);
        setPhi(4, sqrt(mu / getInvB(4) * (1. - getGammaB() * getGammaB())));
        for (int l = 5; l < 10; l++) {
          setPhi(l, sqrt(mu / getInvB(l) * (1. - getGammaS() * getGammaS())));
        }
        for (int l = 10; l < getNumVels(); l++) {
          setPhi(l, sqrt(mu / getInvB(l)));
        }

        for (int i = 0; i < getMyNi().getItem(0); i++) {
          for (int j = 0; j < getMyNi().getItem(1); j++) {
            for (int k = 0; k < getMyNi().getItem(2); k++) {
              for (int l = 0; l < getNumVels(); l++) {
                (*lbfluid)[i][j][k].setPhiLoc(l,getPhi(l));    // set amplitudes of local fluctuations
              }
            }
          }
        }
				
				if (_myRank == 0) {
					std::cout << "The amplitudes phi_i of the fluctuations have been redefined.\n";
					std::cout << "-------------------------------------\n";
				}
      }
    }
		
/*******************************************************************************************/
		
    /* MAKE ONE LB STEP. PUSH-PULL SCHEME IS USED */
    void LatticeBoltzmann::makeLBStep () {
			setStepNum(integrator->getStep());
			int _stepNum = getStepNum();
			int _nSteps = getNSteps();
			int _profStep = getProfStep();
			
			if (getCouplForceFlag() == 1) {
				makeDecompose();
				coupleLBtoMD ();
			}
			
			if (_stepNum % _profStep == 0 && _stepNum!=0) {
				printf ("CPU %d: colstr takes %f sec, comm % f, swapping %f\n",
								getSystem()->comm->rank(), time_colstr, time_comm, time_sw);

				colstream.reset();
				comm.reset();
				swapping.reset();
				time_colstr = 0.;
				time_comm = 0.;
				time_sw = 0.;
			}
			
			if (_stepNum % _nSteps == 0) {
				/* PUSH-scheme (first collide then stream) */
				collideStream ();
			}
		}
		
/*******************************************************************************************/

		/* REAL MD-PARTICLES DECOMPOSITION IF THEY MOVED TO THE GHOST REGION */
		///* it is needed as we couple LB to MD at half-timestep and, most importantly, because real particles can still
		/// leave the real regions and reside in ghost layer. Tricky. Have to understand it better. */
		void LatticeBoltzmann::makeDecompose() {
			int _offset = getHaloSkin();
			real _a = getA();
			Int3D _myNi = getMyNi();
			Real3D _myLeft = getMyLeft();
			
			System& system = getSystemRef();

			CellList realCells = system.storage->getRealCells();
			
			int makeDecompose = 0;
			int totalDecompose = 0;
			for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
				if (cit->position().getItem(0)< _myLeft[0] ||
						cit->position().getItem(1)< _myLeft[1] ||
						cit->position().getItem(2)< _myLeft[2] ||
						cit->position().getItem(0) + 2 * _offset > _myLeft[0] + _myNi[0]*_a ||
						cit->position().getItem(1) + 2 * _offset > _myLeft[1] + _myNi[1]*_a ||
						cit->position().getItem(2) + 2 * _offset > _myLeft[2] + _myNi[2]*_a) {
					makeDecompose = 1;
				}
			}
			
			boost::mpi::all_reduce(*getSystem()->comm, makeDecompose, totalDecompose, std::plus<int>());
			
			if (totalDecompose != 0) {
				system.storage->decompose();
			}
			
		}
		
/*******************************************************************************************/
		
		/* FIND AND OUTPUT CENTER-OF-MASS VELOCITY OF MD-PARTICLES */
		Real3D LatticeBoltzmann::findCMVelMD (int _id) {
			System& system = getSystemRef();
			
			CellList realCells = system.storage->getRealCells();
			
			int _totPart = getTotNPart();
			Real3D myVelCM = Real3D(0.,0.,0.);
			Real3D specVelCM = Real3D(0.,0.,0.);
			
			// loop over all particles in the current CPU
			for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
				Real3D& vel = cit->velocity();
				myVelCM += vel;
			}

			Real3D velCM = Real3D(0.,0.,0.);
			mpi::all_reduce(*getSystem()->comm, myVelCM, velCM, std::plus<Real3D>());
			
			if (getSystem()->comm->rank() == 0) {
				// output of if needed
				if (_id == 1) {
					printf("findCMVelMD: cmV(t+ 1/2dt) of LJ system is %18.14f %18.14f %18.14f \n",
								 velCM.getItem(0), velCM.getItem(1), velCM.getItem(2));
				} else if (_id == 2) {
					printf("findCMVelMD: cmV(t + dt) of LJ system is   %18.14f %18.14f %18.14f \n",
								 velCM.getItem(0), velCM.getItem(1), velCM.getItem(2));
				} else {
				}
			}
			
			// calculate specific center of mass to be subtracted from particle's velocities
			specVelCM = velCM / _totPart;
			
			return specVelCM;
		}
		
/*******************************************************************************************/
		
		/* SET CM VELOCITY OF THE MD TO ZERO AT THE START OF COUPLING */
		void LatticeBoltzmann::zeroMDCMVel () {
			int _myRank = getSystem()->comm->rank();
			setCopyTimestep(integrator->getTimeStep());	// set a copy of a timestep to the real MD timestep

			setStepNum(integrator->getStep());
			if (getStepNum()!=0) setStart(1);
			
			if (getStart() == 0 && getCouplForceFlag() != 0) {
				Real3D specCmVel = findCMVelMD(0);
				// output reporting on subtraction of drift's vel
				if (_myRank == 0) {
					printf("subtracting drift velocity from MD's center of mass (if any)\n");
					printf("cm velocity per particle is %18.14f %18.14f %18.14f \n",
								 specCmVel[0], specCmVel[1], specCmVel[2]);
				}
				
				galileanTransf(specCmVel);
		 
				specCmVel = findCMVelMD(0);
				// check if everything worked correctly
				if (_myRank == 0) {
					printf("cm velocity per particle after Galilean transformation is %18.14f %18.14f %18.14f \n",
								 specCmVel[0], specCmVel[1], specCmVel[2]);
					printf("-------------------------------------\n");
				}

				setStart(1);
			} else if (getStart() == 1 && getCouplForceFlag() != 0) {
				readCouplForces();
				restoreLBForces();
			} else {
			}
		}
		
/*******************************************************************************************/
		
		/* PERFORM GALILEAN TRANSFORMATION */
		void LatticeBoltzmann::galileanTransf (Real3D _specCmVel) {
			System& system = getSystemRef();
			
			CellList realCells = system.storage->getRealCells();
			
			for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
				cit->velocity() -= _specCmVel;
			}
		}
		
/*******************************************************************************************/
		
		/* COLLIDE-STREAM STEP */
    void LatticeBoltzmann::collideStream () {
			int _offset = getHaloSkin();
			int _extForceFlag = getExtForceFlag();
			int _couplForceFlag = getCouplForceFlag();
			int _lbTempFlag = getLBTempFlag();
			Int3D _myNi = getMyNi();
			
			/* copy forces from halo region to the real one */
			if (getCouplForceFlag() == 1) {
				copyForcesFromHalo();
			}
			
			/* collision-streaming */
			real time1 = colstream.getElapsedTime();
			for (int i = _offset; i < _myNi[0]-_offset; i++) {
        for (int j = _offset; j < _myNi[1]-_offset; j++) {
          for (int k = _offset; k < _myNi[2]-_offset; k++) {
						Real3D _f = (*lbfor)[i][j][k].getExtForceLoc() + (*lbfor)[i][j][k].getCouplForceLoc();
						(*lbfluid)[i][j][k].collision(_lbTempFlag, _extForceFlag, _couplForceFlag, _f);

						streaming (i,j,k);
					}
        }
      }
			time_colstr += (colstream.getElapsedTime()-time1);
	
			/* halo communication */
			real time2 = comm.getElapsedTime();
			commHalo();
			time_comm += (comm.getElapsedTime()-time2);
			
			/* swapping of the pointers to the lattices */
			real time3 = swapping.getElapsedTime();
			lblattice *tmp = lbfluid;
			lbfluid = ghostlat;
			ghostlat = tmp;
			time_sw += (swapping.getElapsedTime()-time3);
			
			/* calculate den and j at the lattice sites in real region and copy them to halo */
#warning: should one cancel this condition if pure lb is in use? or move setCouplForceLoc into the collision loop?
			if (getCouplForceFlag() == 1) {
				/* set to zero coupling forces if the coupling exists */
				for (int i = 0; i < _myNi[0]; i++) {
					for (int j = 0; j < _myNi[1]; j++) {
						for (int k = 0; k < _myNi[2]; k++) {
							(*lbfor)[i][j][k].setCouplForceLoc(Real3D(0.0));
						}
					}
				}
				
				/* calculate den and j at the lattice sites in real region */
				calcDenMom();
				
				/* copy den and j from a real region to halo nodes */
				copyDenMomToHalo();
				
			}
    }
		
/*******************************************************************************************/
		
    /* STREAMING ALONG THE VELOCITY VECTORS. SERIAL */
    void LatticeBoltzmann::streaming(int _i, int _j, int _k) {
      // periodic boundaries are handled separately in commHalo() //
			
      // assign iterations
      int _ip = _i + 1; int _im = _i - 1;
      int _jp = _j + 1; int _jm = _j - 1;
      int _kp = _k + 1; int _km = _k - 1;

      // do not move the staying populations
      (*ghostlat)[_i][_j][_k].setF_i(0,(*lbfluid)[_i][_j][_k].getF_i(0));

      // move populations to the nearest neighbors
      (*ghostlat)[_ip][_j][_k].setF_i(1,(*lbfluid)[_i][_j][_k].getF_i(1));
      (*ghostlat)[_im][_j][_k].setF_i(2,(*lbfluid)[_i][_j][_k].getF_i(2));
      (*ghostlat)[_i][_jp][_k].setF_i(3,(*lbfluid)[_i][_j][_k].getF_i(3));
      (*ghostlat)[_i][_jm][_k].setF_i(4,(*lbfluid)[_i][_j][_k].getF_i(4));
      (*ghostlat)[_i][_j][_kp].setF_i(5,(*lbfluid)[_i][_j][_k].getF_i(5));
      (*ghostlat)[_i][_j][_km].setF_i(6,(*lbfluid)[_i][_j][_k].getF_i(6));

      // move populations to the next-to-the-nearest neighbors
      (*ghostlat)[_ip][_jp][_k].setF_i(7,(*lbfluid)[_i][_j][_k].getF_i(7));
      (*ghostlat)[_im][_jm][_k].setF_i(8,(*lbfluid)[_i][_j][_k].getF_i(8));
      (*ghostlat)[_ip][_jm][_k].setF_i(9,(*lbfluid)[_i][_j][_k].getF_i(9));
      (*ghostlat)[_im][_jp][_k].setF_i(10,(*lbfluid)[_i][_j][_k].getF_i(10));
      (*ghostlat)[_ip][_j][_kp].setF_i(11,(*lbfluid)[_i][_j][_k].getF_i(11));
      (*ghostlat)[_im][_j][_km].setF_i(12,(*lbfluid)[_i][_j][_k].getF_i(12));
      (*ghostlat)[_ip][_j][_km].setF_i(13,(*lbfluid)[_i][_j][_k].getF_i(13));
      (*ghostlat)[_im][_j][_kp].setF_i(14,(*lbfluid)[_i][_j][_k].getF_i(14));
      (*ghostlat)[_i][_jp][_kp].setF_i(15,(*lbfluid)[_i][_j][_k].getF_i(15));
      (*ghostlat)[_i][_jm][_km].setF_i(16,(*lbfluid)[_i][_j][_k].getF_i(16));
      (*ghostlat)[_i][_jp][_km].setF_i(17,(*lbfluid)[_i][_j][_k].getF_i(17));
      (*ghostlat)[_i][_jm][_kp].setF_i(18,(*lbfluid)[_i][_j][_k].getF_i(18));
    }

/*******************************************************************************************/
		
    /* SCHEME OF MD TO LB COUPLING */
    void LatticeBoltzmann::coupleLBtoMD() {
			setExtForceFlag(1);

			System& system = getSystemRef();

			CellList realCells = system.storage->getRealCells();

			// loop over all real particles in the current CPU
			for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
				calcRandForce(*cit);				// calculate random force exerted by the fluid onto a particle
				calcViscForce(*cit);			// calculate viscous drag force exerted by the fluid onto a particle
			}
    }
		
/*******************************************************************************************/
		
		void LatticeBoltzmann::calcRandForce (Particle& p) {
			real _fricCoeff = getFricCoeff();							// friction coefficient of LB-coupling
			real _timestep = integrator->getTimeStep();		// timestep of MD
			real _tempLB = getLBTemp();
			
			real prefactor = sqrt(24. * _fricCoeff * _tempLB / _timestep);		// amplitude of the noise
			Real3D ranval((*rng)() - 0.5, (*rng)() - 0.5, (*rng)() - 0.5);		// 3d random number
			setFOnPart (p.id(), prefactor * ranval);			// set force on a particle to the random one
		}
		
/*******************************************************************************************/
		
		void LatticeBoltzmann::calcViscForce (Particle& p) {
			int _offset = getHaloSkin();
			real _a = getA();
			real _invA = 1. / _a;
			real _fricCoeff = getFricCoeff();
			Real3D Li = getSystem()->bc->getBoxL();

			// account for particle's positions with respect to CPU's left border
			Real3D _pos = p.position() - getMyLeft();

			Real3D _posLB = _pos;
			_posLB += (double)_offset;
			_posLB *= _invA;
			
			Int3D bin;
			bin[0] = floor (_posLB[0]); bin[1] = floor (_posLB[1]); bin[2] = floor (_posLB[2]);

			// weight factors, dimensionless
			std::vector<real> delta = std::vector<real>(6, 0.);
			delta[0] = _posLB[0] - bin[0];
			delta[1] = _posLB[1] - bin[1]; 
			delta[2] = _posLB[2] - bin[2];
			delta[3] = _a - delta[0];
			delta[4] = _a - delta[1];
			delta[5] = _a - delta[2];
			
			real _convTimeMDtoLB = convTimeMDtoLB();
			real _convLenMDtoLB = convLenMDtoLB();
			real _convMassMDtoLB = convMassMDtoLB();
			int _numVels = getNumVels();
			real _invDenLoc;
			Real3D _jLoc;
			Real3D _u;
			
			Real3D interpVel = Real3D (0.);
			real _convCoeff = _convTimeMDtoLB / _convLenMDtoLB;
			// loop over neighboring LB nodes
			int _ip, _jp, _kp;
			for (int _i = 0; _i < 2; _i++) {
				for (int _j = 0; _j < 2; _j++) {
					for (int _k = 0; _k < 2; _k++) {
						// assign iterations
						_ip = bin[0] + _i; _jp = bin[1] + _j; _kp = bin[2] + _k;

 						_invDenLoc = 1. / (*lbmom)[_ip][_jp][_kp].getMom_i(0);
						_jLoc[0] = (*lbmom)[_ip][_jp][_kp].getMom_i(1);
						_jLoc[1] = (*lbmom)[_ip][_jp][_kp].getMom_i(2);
						_jLoc[2] = (*lbmom)[_ip][_jp][_kp].getMom_i(3);

						_u = _jLoc;
						_u *= _invDenLoc;
						_u *= _convCoeff;

						interpVel += _u * delta[3 * _i] * delta[3 * _j + 1] * delta[3 * _k + 2];
					}
				}
			}

			// add viscous force to the buffered random force acting onto particle p.id()
			addFOnPart(p.id(), -_fricCoeff * (p.velocity() - interpVel));
			
			// apply buffered force to the MD-particle p.id()
			p.force() += getFOnPart(p.id());

			// convert coupling force (LJ units) to the momentum change on a lattice (LB units)
			Real3D deltaJLoc = Real3D(0.);
			deltaJLoc -= getFOnPart(p.id()) * _convMassMDtoLB / (_convCoeff * _convTimeMDtoLB);
			
			// loop over neighboring LB nodes
			for (int _i = 0; _i < 2; _i++) {
				for (int _j = 0; _j < 2; _j++) {
					for (int _k = 0; _k < 2; _k++) {
						// periodic boundaries on the right side
						_ip = bin[0] + _i; _jp = bin[1] + _j; _kp = bin[2] + _k;
					
						// converting momentum into coupling force with weights delta[i]
						Real3D _fLoc = deltaJLoc;
						_fLoc *= delta[3*_i]; _fLoc *= delta[3*_j+1]; _fLoc *= delta[3*_k+2];
						
						// add coupling force to the correspondent lattice cite
						(*lbfor)[_ip][_jp][_kp].addCouplForceLoc(_fLoc);
					}
				}
			}
		}
		
/*******************************************************************************************/
		
//		void LatticeBoltzmann::addLBForces (Particle& p) {
//			// apply random and viscous force to MD the particle p.id()
//			p.force() += getFOnPart(p.id());
//		}
		
/*******************************************************************************************/
		
		/* RESTORING FORCES ACTING FROM LB FLUID ONTO MD PARTICLES */
		void LatticeBoltzmann::restoreLBForces () {
			System& system = getSystemRef();
			
			CellList realCells = system.storage->getRealCells();
	
			// loop over all particles in the current CPU
			for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
				cit->force() += getFOnPart(cit->id());
			}
		}
		
/*******************************************************************************************/
		
		/* CALCULATE DENSITY AND J AT THE LATTICE SITES IN REAL REGION */
		void LatticeBoltzmann::calcDenMom () {
			Int3D _myNi = getMyNi();
			int _numVels = getNumVels();
			int _offset = getHaloSkin();
			
			for (int i = _offset; i<_myNi[0]-_offset; ++i) {
				for (int j = _offset; j<_myNi[1]-_offset; ++j) {
					for (int k = _offset; k<_myNi[2]-_offset; ++k) {
						real denLoc = 0.;
						Real3D jLoc = Real3D(0.);
						for (int l = 0; l < _numVels; l++) {
							denLoc += (*lbfluid)[i][j][k].getF_i(l);
							jLoc += (*lbfluid)[i][j][k].getF_i(l)*getCi(l);
						}
						(*lbmom)[i][j][k].setMom_i(0,denLoc);
						(*lbmom)[i][j][k].setMom_i(1,jLoc[0]);
						(*lbmom)[i][j][k].setMom_i(2,jLoc[1]);
						(*lbmom)[i][j][k].setMom_i(3,jLoc[2]);
					}
				}
			}
		}
		
/*******************************************************************************************/
		
		////////////////////////////////////////
		/* HANDLING I/O AND RESTARTS */
		///////////////////////////////////////
		
/*******************************************************************************************/
		
		void LatticeBoltzmann::readCouplForces () {
			timeRead.reset();
			real timeStart = timeRead.getElapsedTime();
			
			/* fill in the coupling forces acting on MD-particles with zeros */
			int _totNPart = getTotNPart();
			for(int _id = 0; _id <= _totNPart; _id++) {
				setFOnPart(_id, Real3D(0.));
			}
			
			/* create filename for the input file */
			std::string filename;
			std::ostringstream convert;
			std::ostringstream _myRank;
			
			convert << getStepNum();
			_myRank << getSystem()->comm->rank();
			
			filename = "couplForces";
			filename.append(convert.str());
			filename.append(".");
			filename.append(_myRank.str());
			filename.append(".dat");
			
			/* access particles' data and open a file to read coupling forces from */
			long int _id;
			real _fx, _fy, _fz;
			
			System& system = getSystemRef();
			CellList realCells = system.storage->getRealCells();
			FILE * couplForcesFile = fopen(filename.c_str(),"r");
			
			if (couplForcesFile == NULL) {
				if (getStepNum() != 0) {
					std::cout << "!!! Attention !!! no file with coupling forces acting onto MD particles \
					found for step " << convert.str() << "\n";
				}
			} else {
				/* read forces acting onto MD-particles */
				for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
					fscanf (couplForcesFile, "%ld %lf %lf %lf \n", &_id, &_fx, &_fy, &_fz);
					setFOnPart(_id, Real3D(_fx,_fy,_fz));
				}
				/* read forces acting onto LB-sites */
				int _it, _jt, _kt;
				Int3D _myNi = getMyNi();

				for (int _i = 0; _i < _myNi[0]; _i++) {
					for (int _j = 0; _j < _myNi[1]; _j++) {
						for (int _k = 0; _k < _myNi[2]; _k++) {
//							printf("CPU %d: (*lbfluid)[_i][_j][_k].getCouplForceLoc() is %20.14f\n", getSystem()->comm->rank(), (*lbfluid)[_i][_j][_k].getCouplForceLoc().getItem(2));
						(*lbfor)[_i][_j][_k].setCouplForceLoc(Real3D(0.));
						}
					}
				}

				while (fscanf (couplForcesFile, "%d %d %d %lf %lf %lf \n", &_it, &_jt, &_kt, &_fx, &_fy, &_fz) == 6) {
					(*lbfor)[_it][_jt][_kt].setCouplForceLoc(Real3D(_fx,_fy,_fz));
				}
				
			}
			fclose (couplForcesFile);
			
			real timeEnd = timeRead.getElapsedTime() - timeStart;
			printf("CPU %d: read LB-to-MD coupling forces in %8.4f seconds\n",
						 getSystem()->comm->rank(), timeEnd);
		}
		
/*******************************************************************************************/
		
		void LatticeBoltzmann::saveCouplForces () {
			timeSave.reset();
			real timeStart = timeSave.getElapsedTime();
			
			/* create filename for the output file */
			std::string filename;
			std::ostringstream convert;
			std::ostringstream _myRank;
			
			convert << integrator->getStep();
			_myRank << getSystem()->comm->rank();

			filename = "couplForces";
			filename.append(convert.str());
			filename.append(".");
			filename.append(_myRank.str());
			filename.append(".dat");
			
			/* access particles' data and open a file to write coupling forces to */
			System& system = getSystemRef();
			CellList realCells = system.storage->getRealCells();

			/* write forces acting onto MD-particles */
			FILE * couplForcesFile = fopen(filename.c_str(),"w");

			for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
				long int _id = cit->id();
				fprintf (couplForcesFile, "%9ld %20.16lf %20.16lf %20.16lf \n", _id,
								 getFOnPart(_id).getItem(0),
								 getFOnPart(_id).getItem(1),
								 getFOnPart(_id).getItem(2));
			}

			/* write forces acting onto LB-sites */
			real threshold = 1e-30;													// threshold for output forces
			Int3D _myNi = getMyNi();
			
			for (int _i = 0; _i < _myNi[0]; _i++) {
				for (int _j = 0; _j < _myNi[1]; _j++) {
					for (int _k = 0; _k < _myNi[2]; _k++) {
						Real3D _couplForceLoc = (*lbfor)[_i][_j][_k].getCouplForceLoc();
						if (_couplForceLoc.sqr() < threshold) {
							// do not output zero-like forces
						} else {
							fprintf (couplForcesFile, "%5d %5d %5d %20.16f %20.16f %20.16f \n",
											 _i, _j, _k, _couplForceLoc[0],
											 _couplForceLoc[1], _couplForceLoc[2]);
						}
					}
				}
			}
			fclose (couplForcesFile);
			
			real timeEnd = timeSave.getElapsedTime() - timeStart;
			printf("CPU %d: saved LB-to-MD coupling forces in %8.4f seconds\n",
						 getSystem()->comm->rank(), timeEnd);
		}
		
/*******************************************************************************************/
		
		///////////////////////////
		/* PARALLELISATION */
		///////////////////////////
		
/*******************************************************************************************/
		
		/* FIND RANKS OF NEIGHBOURUNG CPU IN 6 DIRECTIONS */
		void LatticeBoltzmann::findMyNeighbours () {
			
			/* calculate dimensionality of the processors grid arrangement */
			int nodeGridDim = 0;
			Int3D _nodeGrid = getNodeGrid();
			for (int i = 0; i < 3; ++i) {
				if (_nodeGrid[i] != 1) {
					++nodeGridDim;
				}
			}
			
			/* define myRank and myPosition in the nodeGrid */
			longint _myRank = getSystem()->comm->rank();
			Int3D _myPosition = Int3D(0,0,0);
			
			esutil::Grid grid(_nodeGrid);
			grid.mapIndexToPosition(_myPosition, _myRank);
			setMyPosition(_myPosition);
			
			/* calculate ranks of neighbouring processors in every direction */
			Int3D _myNeighbourPos = Int3D(0,0,0);
			for (int _dim = 0; _dim < nodeGridDim; ++_dim) {
				for (int _j = 0; _j < 3; ++_j) {
					_myNeighbourPos[_j] = _myPosition[_j];
				}
				
				// left neighbor in direction _dim (x, y or z)
				_myNeighbourPos[_dim] = _myPosition[_dim] - 1;
				if (_myNeighbourPos[_dim] < 0) {
					_myNeighbourPos[_dim] += _nodeGrid[_dim];
				}
				setMyNeighbour(2*_dim, grid.mapPositionToIndex(_myNeighbourPos));
				
				// right neighbor in direction _dim (x, y or z)
				_myNeighbourPos[_dim] = _myPosition[_dim] + 1;
				if (_myNeighbourPos[_dim] >= _nodeGrid[_dim]) {
					_myNeighbourPos[_dim] -= _nodeGrid[_dim];
				}
				setMyNeighbour(2*_dim + 1, grid.mapPositionToIndex(_myNeighbourPos));
			}
			
			if (_myRank == 0) {
				printf ("Number of CPUs in use is %d\n", mpiWorld->size());
			}
			
			for (int _dim = nodeGridDim; _dim < 3; ++_dim) {
				setMyNeighbour(2*_dim, grid.mapPositionToIndex(_myPosition));
				setMyNeighbour(2*_dim+1, grid.mapPositionToIndex(_myPosition));
			}
			
		}
		
/*******************************************************************************************/
		
		/* COMMUNICATE POPULATIONS IN HALO REGIONS TO THE NEIGHBOURING CPUs */
		void LatticeBoltzmann::commHalo() {
			int i, j, k, index;											// running indices and index of the node to be copied
			static int numPopTransf = 5;						// number of populations and hydrod. moments to transfer
			int rnode, snode;												// rank of the node to receive from and to send to
			std::vector<real> bufToSend, bufToRecv;	// buffers used to send and to receive the data
			
			int _offset = getHaloSkin();
			Int3D _myNi = getMyNi();
			Int3D _myPosition = getMyPosition();
			
			mpi::environment env;
			mpi::communicator world;
			
			//////////////////////
			//// X-direction /////
			//////////////////////
			// number of reals to transfer
			int numDataTransf = numPopTransf * _myNi[1] * _myNi[2];

			bufToSend.resize(numDataTransf);				// resize bufToSend
			bufToRecv.resize(numDataTransf);				// resize bufToRecv
			
			/* send to right, recv from left i = 1, 7, 9, 11, 13 */
			snode = getMyNeighbour(1);
			rnode = getMyNeighbour(0);
			
			// prepare message for sending
			i = _myNi[0] - _offset;
			index = 0;
			for (k=0; k<_myNi[2]; k++) {
				for (j=0; j<_myNi[1]; j++) {
					bufToSend[index] = (*ghostlat)[i][j][k].getF_i(1);
					bufToSend[index+1] = (*ghostlat)[i][j][k].getF_i(7);
					bufToSend[index+2] = (*ghostlat)[i][j][k].getF_i(9);
					bufToSend[index+3] = (*ghostlat)[i][j][k].getF_i(11);
					bufToSend[index+4] = (*ghostlat)[i][j][k].getF_i(13);
					index += numPopTransf;
				}
			}
			
			// send and receive data or use memcpy if number of CPU in x-dir is 1
			if (getNodeGrid().getItem(0) > 1) {
				if (_myPosition[0] % 2 == 0) {
					world.send(snode, COMM_DIR_0, bufToSend);
					world.recv(rnode, COMM_DIR_0, bufToRecv);
				} else {
					world.recv(rnode, COMM_DIR_0, bufToRecv);
					world.send(snode, COMM_DIR_0, bufToSend);
				}
			} else {
				bufToRecv = bufToSend;
			}
			
			// unpack message
			i = _offset;
			index = 0;
			for (k=0; k<_myNi[2]; k++) {
				for (j=0; j<_myNi[1]; j++) {
					(*ghostlat)[i][j][k].setF_i(1, bufToRecv[index]);
					(*ghostlat)[i][j][k].setF_i(7, bufToRecv[index+1]);
					(*ghostlat)[i][j][k].setF_i(9, bufToRecv[index+2]);
					(*ghostlat)[i][j][k].setF_i(11, bufToRecv[index+3]);
					(*ghostlat)[i][j][k].setF_i(13, bufToRecv[index+4]);
					index += numPopTransf;
				}
			}
			
			/* send to left, recv from right i = 2, 8, 10, 12, 14 */
			snode = getMyNeighbour(0);
			rnode = getMyNeighbour(1);
			
			// prepare message for sending
			i = 0;
			index = 0;
			for (k=0; k<_myNi[2]; k++) {
				for (j=0; j<_myNi[1]; j++) {
					bufToSend[index] = (*ghostlat)[i][j][k].getF_i(2);
					bufToSend[index+1] = (*ghostlat)[i][j][k].getF_i(8);
					bufToSend[index+2] = (*ghostlat)[i][j][k].getF_i(10);
					bufToSend[index+3] = (*ghostlat)[i][j][k].getF_i(12);
					bufToSend[index+4] = (*ghostlat)[i][j][k].getF_i(14);
					index += numPopTransf;
				}
			}
			
			// send and receive data or use memcpy if number of CPU in x-dir is 1
			if (getNodeGrid().getItem(0) > 1) {
				if (_myPosition[0] % 2 == 0) {
					world.send(snode, COMM_DIR_1, bufToSend);
					world.recv(rnode, COMM_DIR_1, bufToRecv);
				} else {
					world.recv(rnode, COMM_DIR_1, bufToRecv);
					world.send(snode, COMM_DIR_1, bufToSend);
				}
			} else {
				bufToRecv = bufToSend;
			}
			
			// unpack message
			i = _myNi[0] - 2 * _offset;
			index = 0;
			for (k=0; k<_myNi[2]; k++) {
				for (j=0; j<_myNi[1]; j++) {
					(*ghostlat)[i][j][k].setF_i(2, bufToRecv[index]);
					(*ghostlat)[i][j][k].setF_i(8, bufToRecv[index+1]);
					(*ghostlat)[i][j][k].setF_i(10, bufToRecv[index+2]);
					(*ghostlat)[i][j][k].setF_i(12, bufToRecv[index+3]);
					(*ghostlat)[i][j][k].setF_i(14, bufToRecv[index+4]);
					index += numPopTransf;
				}
			}
			
			//////////////////////
			//// Y-direction /////
			//////////////////////
			numDataTransf = numPopTransf * _myNi[0] * _myNi[2];
			bufToSend.resize(numDataTransf);				// resize bufToSend
			bufToRecv.resize(numDataTransf);				// resize bufToRecv
			
			/* send to right, recv from left i = 3, 7, 10, 15, 17 */
			snode = getMyNeighbour(3);
			rnode = getMyNeighbour(2);
			
			// prepare message for sending
			j = _myNi[1] - _offset;
			index = 0;
			for (k=0; k<_myNi[2]; k++) {
				for (i=0; i<_myNi[0]; i++) {
					bufToSend[index] = (*ghostlat)[i][j][k].getF_i(3);
					bufToSend[index+1] = (*ghostlat)[i][j][k].getF_i(7);
					bufToSend[index+2] = (*ghostlat)[i][j][k].getF_i(10);
					bufToSend[index+3] = (*ghostlat)[i][j][k].getF_i(15);
					bufToSend[index+4] = (*ghostlat)[i][j][k].getF_i(17);
					index += numPopTransf;
				}
			}
			
			// send and receive data or use memcpy if number of CPU in y-dir is 1
			if (getNodeGrid().getItem(1) > 1) {
				if (_myPosition[1] % 2 == 0) {
					world.send(snode, COMM_DIR_2, bufToSend);
					world.recv(rnode, COMM_DIR_2, bufToRecv);
				} else {
					world.recv(rnode, COMM_DIR_2, bufToRecv);
					world.send(snode, COMM_DIR_2, bufToSend);
				}
			} else {
				bufToRecv = bufToSend;
			}
			
			// unpack message
			j = _offset;
			index = 0;
			for (k=0; k<_myNi[2]; k++) {
				for (i=0; i<_myNi[0]; i++) {
					(*ghostlat)[i][j][k].setF_i(3, bufToRecv[index]);
					(*ghostlat)[i][j][k].setF_i(7, bufToRecv[index+1]);
					(*ghostlat)[i][j][k].setF_i(10, bufToRecv[index+2]);
					(*ghostlat)[i][j][k].setF_i(15, bufToRecv[index+3]);
					(*ghostlat)[i][j][k].setF_i(17, bufToRecv[index+4]);
					index += numPopTransf;
				}
			}
			
			/* send to left, recv from right i = 4, 8, 9, 16, 18 */
			snode = getMyNeighbour(2);
			rnode = getMyNeighbour(3);
			
			// prepare message for sending
			j = 0;
			index = 0;
			for (k=0; k<_myNi[2]; k++) {
				for (i=0; i<_myNi[0]; i++) {
					bufToSend[index] = (*ghostlat)[i][j][k].getF_i(4);
					bufToSend[index+1] = (*ghostlat)[i][j][k].getF_i(8);
					bufToSend[index+2] = (*ghostlat)[i][j][k].getF_i(9);
					bufToSend[index+3] = (*ghostlat)[i][j][k].getF_i(16);
					bufToSend[index+4] = (*ghostlat)[i][j][k].getF_i(18);
					index += numPopTransf;
				}
			}
			
			// send and receive data or use memcpy if number of CPU in y-dir is 1
			if (getNodeGrid().getItem(1) > 1) {
				if (_myPosition[1] % 2 == 0) {
					world.send(snode, COMM_DIR_3, bufToSend);
					world.recv(rnode, COMM_DIR_3, bufToRecv);
				} else {
					world.recv(rnode, COMM_DIR_3, bufToRecv);
					world.send(snode, COMM_DIR_3, bufToSend);
				}
			} else {
				bufToRecv = bufToSend;
			}
			
			// unpack message
			j = _myNi[1] - 2 * _offset;
			index = 0;
			for (k=0; k<_myNi[2]; k++) {
				for (i=0; i<_myNi[0]; i++) {
					(*ghostlat)[i][j][k].setF_i(4, bufToRecv[index]);
					(*ghostlat)[i][j][k].setF_i(8, bufToRecv[index+1]);
					(*ghostlat)[i][j][k].setF_i(9, bufToRecv[index+2]);
					(*ghostlat)[i][j][k].setF_i(16, bufToRecv[index+3]);
					(*ghostlat)[i][j][k].setF_i(18, bufToRecv[index+4]);
					index += numPopTransf;
				}
			}
			
			//////////////////////
			//// Z-direction /////
			//////////////////////
			numDataTransf = numPopTransf * _myNi[0] * _myNi[1];
			bufToSend.resize(numDataTransf);				// resize bufToSend
			bufToRecv.resize(numDataTransf);				// resize bufToRecv
			
			/* send to right, recv from left i = 5, 11, 14, 15, 18 */
			snode = getMyNeighbour(5);
			rnode = getMyNeighbour(4);
			
			// prepare message for sending
			k = _myNi[2] - _offset;
			index = 0;
			for (j=0; j<_myNi[1]; j++) {
				for (i=0; i<_myNi[0]; i++) {
					bufToSend[index] = (*ghostlat)[i][j][k].getF_i(5);
					bufToSend[index+1] = (*ghostlat)[i][j][k].getF_i(11);
					bufToSend[index+2] = (*ghostlat)[i][j][k].getF_i(14);
					bufToSend[index+3] = (*ghostlat)[i][j][k].getF_i(15);
					bufToSend[index+4] = (*ghostlat)[i][j][k].getF_i(18);
					index += numPopTransf;
				}
			}
			
			// send and receive data or use memcpy if number of CPU in z-dir is 1
			if (getNodeGrid().getItem(2) > 1) {
				if (_myPosition[2] % 2 == 0) {
					world.send(snode, COMM_DIR_4, bufToSend);
					world.recv(rnode, COMM_DIR_4, bufToRecv);
				} else {
					world.recv(rnode, COMM_DIR_4, bufToRecv);
					world.send(snode, COMM_DIR_4, bufToSend);
				}
			} else {
				bufToRecv = bufToSend;
			}
			
			// unpack message
			k = _offset;
			index = 0;
			for (j=0; j<_myNi[1]; j++) {
				for (i=0; i<_myNi[0]; i++) {
					(*ghostlat)[i][j][k].setF_i(5, bufToRecv[index]);
					(*ghostlat)[i][j][k].setF_i(11, bufToRecv[index+1]);
					(*ghostlat)[i][j][k].setF_i(14, bufToRecv[index+2]);
					(*ghostlat)[i][j][k].setF_i(15, bufToRecv[index+3]);
					(*ghostlat)[i][j][k].setF_i(18, bufToRecv[index+4]);
					index += numPopTransf;
				}
			}
			
			/* send to left, recv from right i = 6, 12, 13, 16, 17 */
			snode = getMyNeighbour(4);
			rnode = getMyNeighbour(5);
			
			// prepare message for sending
			k = 0;
			index = 0;
			for (j=0; j<_myNi[1]; j++) {
				for (i=0; i<_myNi[0]; i++) {
					bufToSend[index] = (*ghostlat)[i][j][k].getF_i(6);
					bufToSend[index+1] = (*ghostlat)[i][j][k].getF_i(12);
					bufToSend[index+2] = (*ghostlat)[i][j][k].getF_i(13);
					bufToSend[index+3] = (*ghostlat)[i][j][k].getF_i(16);
					bufToSend[index+4] = (*ghostlat)[i][j][k].getF_i(17);
					index += numPopTransf;
				}
			}
			
			// send and receive data or use memcpy if number of CPU in z-dir is 1
			if (getNodeGrid().getItem(2) > 1) {
				if (_myPosition[2] % 2 == 0) {
					world.send(snode, COMM_DIR_5, bufToSend);
					world.recv(rnode, COMM_DIR_5, bufToRecv);
				} else {
					world.recv(rnode, COMM_DIR_5, bufToRecv);
					world.send(snode, COMM_DIR_5, bufToSend);
				}
			} else {
				bufToRecv = bufToSend;
			}
			
			// unpack message
			k = _myNi[2] - 2 * _offset;
			index = 0;
			for (j=0; j<_myNi[1]; j++) {
				for (i=0; i<_myNi[0]; i++) {
					(*ghostlat)[i][j][k].setF_i(6, bufToRecv[index]);
					(*ghostlat)[i][j][k].setF_i(12, bufToRecv[index+1]);
					(*ghostlat)[i][j][k].setF_i(13, bufToRecv[index+2]);
					(*ghostlat)[i][j][k].setF_i(16, bufToRecv[index+3]);
					(*ghostlat)[i][j][k].setF_i(17, bufToRecv[index+4]);
					index += numPopTransf;
				}
			}
			
			bufToSend.resize(0);
			bufToRecv.resize(0);
		}
		
/*******************************************************************************************/
		
		/* COPY COUPLING FORCES FROM HALO REGIONS TO THE REAL ONES */
		void LatticeBoltzmann::copyForcesFromHalo () {
			int i, j, k, index;											// running indices and index of the node to be copied
			int numForceComp = 3;										// number of force components to transfer
			int numDataTransf;											// number of data to transfer
			int rnode, snode;												// rank of the node to receive from and to send to
			std::vector<real> bufToSend, bufToRecv;	// buffers used to send and to receive the data
			Real3D _addForce;
			
			int _offset = getHaloSkin();
			Int3D _myNi = getMyNi();
			Int3D _myPosition = getMyPosition();
			
			mpi::environment env;
			mpi::communicator world;
			
			//////////////////////
			//// X-direction /////
			//////////////////////
			numDataTransf = numForceComp * _myNi[1] * _myNi[2];
			bufToSend.resize(numDataTransf);				// resize bufToSend
			bufToRecv.resize(numDataTransf);				// resize bufToRecv
			
			/* send to right, recv from left */
			snode = getMyNeighbour(1);
			rnode = getMyNeighbour(0);
			
			// prepare message for sending
			i = _myNi[0] - _offset;
			for (k=0; k<_myNi[2]; k++) {
				for (j=0; j<_myNi[1]; j++) {
					index = numForceComp*_myNi[1]*k + j*numForceComp;
					
					bufToSend[index] = (*lbfor)[i][j][k].getCouplForceLoc().getItem(0);
					bufToSend[index+1] = (*lbfor)[i][j][k].getCouplForceLoc().getItem(1);
					bufToSend[index+2] = (*lbfor)[i][j][k].getCouplForceLoc().getItem(2);
				}
			}
			
			// send and receive data or use memcpy if number of CPU in x-dir is 1
			if (getNodeGrid().getItem(0) > 1) {
				if (_myPosition[0] % 2 == 0) {
					world.send(snode, COMM_FORCE_0, bufToSend);
					world.recv(rnode, COMM_FORCE_0, bufToRecv);
				} else {
					world.recv(rnode, COMM_FORCE_0, bufToRecv);
					world.send(snode, COMM_FORCE_0, bufToSend);
				}
			} else {
				bufToRecv = bufToSend;
			}
			
			// unpack message
			i = _offset;
			for (k=0; k<_myNi[2]; k++) {
				for (j=0; j<_myNi[1]; j++) {
					index = numForceComp*_myNi[1]*k + j*numForceComp;
					_addForce = Real3D(bufToRecv[index], bufToRecv[index+1], bufToRecv[index+2]);
					(*lbfor)[i][j][k].addCouplForceLoc(_addForce);
				}
			}
			
			/* send to left, recv from right */
			snode = getMyNeighbour(0);
			rnode = getMyNeighbour(1);
			
			// prepare message for sending
			i = 0;
			for (k=0; k<_myNi[2]; k++) {
				for (j=0; j<_myNi[1]; j++) {
					index = numForceComp*_myNi[1]*k + j*numForceComp;
					
					bufToSend[index] = (*lbfor)[i][j][k].getCouplForceLoc().getItem(0);
					bufToSend[index+1] = (*lbfor)[i][j][k].getCouplForceLoc().getItem(1);
					bufToSend[index+2] = (*lbfor)[i][j][k].getCouplForceLoc().getItem(2);
				}
			}
			
			// send and receive data or use memcpy if number of CPU in x-dir is 1
			if (getNodeGrid().getItem(0) > 1) {
				if (_myPosition[0] % 2 == 0) {
					world.send(snode, COMM_FORCE_1, bufToSend);
					world.recv(rnode, COMM_FORCE_1, bufToRecv);
				} else {
					world.recv(rnode, COMM_FORCE_1, bufToRecv);
					world.send(snode, COMM_FORCE_1, bufToSend);
				}
			} else {
				bufToRecv = bufToSend;
			}
			
			// unpack message
			i = _myNi[0] - 2 * _offset;
			for (k=0; k<_myNi[2]; k++) {
				for (j=0; j<_myNi[1]; j++) {
					index = numForceComp*_myNi[1]*k + j*numForceComp;
					_addForce = Real3D(bufToRecv[index], bufToRecv[index+1], bufToRecv[index+2]);
					(*lbfor)[i][j][k].addCouplForceLoc(_addForce);
				}
			}

			//////////////////////
			//// Y-direction /////
			//////////////////////
			numDataTransf = numForceComp * _myNi[0] * _myNi[2];
			bufToSend.resize(numDataTransf);				// resize bufToSend
			bufToRecv.resize(numDataTransf);				// resize bufToRecv
			
			/* send to right, recv from left */
			snode = getMyNeighbour(3);
			rnode = getMyNeighbour(2);
			
			// prepare message for sending
			j = _myNi[1] - _offset;
			for (k=0; k<_myNi[2]; k++) {
				for (i=0; i<_myNi[0]; i++) {
					index = numForceComp*_myNi[0]*k + i*numForceComp;
					
					bufToSend[index] = (*lbfor)[i][j][k].getCouplForceLoc().getItem(0);
					bufToSend[index+1] = (*lbfor)[i][j][k].getCouplForceLoc().getItem(1);
					bufToSend[index+2] = (*lbfor)[i][j][k].getCouplForceLoc().getItem(2);
				}
			}
			
			// send and receive data or use memcpy if number of CPU in y-dir is 1
			if (getNodeGrid().getItem(1) > 1) {
				if (_myPosition[1] % 2 == 0) {
					world.send(snode, COMM_FORCE_2, bufToSend);
					world.recv(rnode, COMM_FORCE_2, bufToRecv);
				} else {
					world.recv(rnode, COMM_FORCE_2, bufToRecv);
					world.send(snode, COMM_FORCE_2, bufToSend);
				}
			} else {
				bufToRecv = bufToSend;
			}
			
			// unpack message
			j = _offset;
			for (k=0; k<_myNi[2]; k++) {
				for (i=0; i<_myNi[0]; i++) {
					index = numForceComp*_myNi[0]*k + i*numForceComp;
					_addForce = Real3D(bufToRecv[index], bufToRecv[index+1], bufToRecv[index+2]);
					(*lbfor)[i][j][k].addCouplForceLoc(_addForce);
				}
			}
			
			/* send to left, recv from right */
			snode = getMyNeighbour(2);
			rnode = getMyNeighbour(3);
			
			// prepare message for sending
			j = 0;
			for (k=0; k<_myNi[2]; k++) {
				for (i=0; i<_myNi[0]; i++) {
					index = numForceComp*_myNi[0]*k + i*numForceComp;
					
					bufToSend[index] = (*lbfor)[i][j][k].getCouplForceLoc().getItem(0);
					bufToSend[index+1] = (*lbfor)[i][j][k].getCouplForceLoc().getItem(1);
					bufToSend[index+2] = (*lbfor)[i][j][k].getCouplForceLoc().getItem(2);
				}
			}
			
			// send and receive data or use memcpy if number of CPU in y-dir is 1
			if (getNodeGrid().getItem(1) > 1) {
				if (_myPosition[1] % 2 == 0) {
					world.send(snode, COMM_FORCE_3, bufToSend);
					world.recv(rnode, COMM_FORCE_3, bufToRecv);
				} else {
					world.recv(rnode, COMM_FORCE_3, bufToRecv);
					world.send(snode, COMM_FORCE_3, bufToSend);
				}
			} else {
				bufToRecv = bufToSend;
			}
			
			// unpack message
			j = _myNi[1] - 2 * _offset;
			for (k=0; k<_myNi[2]; k++) {
				for (i=0; i<_myNi[0]; i++) {
					index = numForceComp*_myNi[0]*k + i*numForceComp;
					_addForce = Real3D(bufToRecv[index], bufToRecv[index+1], bufToRecv[index+2]);
					(*lbfor)[i][j][k].addCouplForceLoc(_addForce);
				}
			}
			
			//////////////////////
			//// Z-direction /////
			//////////////////////
			numDataTransf = numForceComp * _myNi[0] * _myNi[1];
			bufToSend.resize(numDataTransf);				// resize bufToSend
			bufToRecv.resize(numDataTransf);				// resize bufToRecv
			
			/* send to right, recv from left */
			snode = getMyNeighbour(5);
			rnode = getMyNeighbour(4);
			
			// prepare message for sending
			k = _myNi[2] - _offset;
			for (j=0; j<_myNi[1]; j++) {
				for (i=0; i<_myNi[0]; i++) {
					index = numForceComp*_myNi[0]*j + i*numForceComp;
					
					bufToSend[index] = (*lbfor)[i][j][k].getCouplForceLoc().getItem(0);
					bufToSend[index+1] = (*lbfor)[i][j][k].getCouplForceLoc().getItem(1);
					bufToSend[index+2] = (*lbfor)[i][j][k].getCouplForceLoc().getItem(2);
				}
			}
			
			// send and receive data or use memcpy if number of CPU in z-dir is 1
			if (getNodeGrid().getItem(2) > 1) {
				if (_myPosition[2] % 2 == 0) {
					world.send(snode, COMM_FORCE_4, bufToSend);
					world.recv(rnode, COMM_FORCE_4, bufToRecv);
				} else {
					world.recv(rnode, COMM_FORCE_4, bufToRecv);
					world.send(snode, COMM_FORCE_4, bufToSend);
				}
			} else {
				bufToRecv = bufToSend;
			}
			
			// unpack message
			k = _offset;
			for (j=0; j<_myNi[1]; j++) {
				for (i=0; i<_myNi[0]; i++) {
					index = numForceComp*_myNi[0]*j + i*numForceComp;
					_addForce = Real3D(bufToRecv[index], bufToRecv[index+1], bufToRecv[index+2]);
					(*lbfor)[i][j][k].addCouplForceLoc(_addForce);
				}
			}
			
			/* send to left, recv from right */
			snode = getMyNeighbour(4);
			rnode = getMyNeighbour(5);
			
			// prepare message for sending
			k = 0;
			for (j=0; j<_myNi[1]; j++) {
				for (i=0; i<_myNi[0]; i++) {
					index = numForceComp*_myNi[0]*j + i*numForceComp;
					
					bufToSend[index] = (*lbfor)[i][j][k].getCouplForceLoc().getItem(0);
					bufToSend[index+1] = (*lbfor)[i][j][k].getCouplForceLoc().getItem(1);
					bufToSend[index+2] = (*lbfor)[i][j][k].getCouplForceLoc().getItem(2);
				}
			}
			
			// send and receive data or use memcpy if number of CPU in z-dir is 1
			if (getNodeGrid().getItem(2) > 1) {
				if (_myPosition[2] % 2 == 0) {
					world.send(snode, COMM_FORCE_5, bufToSend);
					world.recv(rnode, COMM_FORCE_5, bufToRecv);
				} else {
					world.recv(rnode, COMM_FORCE_5, bufToRecv);
					world.send(snode, COMM_FORCE_5, bufToSend);
				}
			} else {
				bufToRecv = bufToSend;
			}
			
			// unpack message
			k = _myNi[2] - 2 * _offset;
			for (j=0; j<_myNi[1]; j++) {
				for (i=0; i<_myNi[0]; i++) {
					index = numForceComp*_myNi[0]*j + i*numForceComp;
					_addForce = Real3D(bufToRecv[index], bufToRecv[index+1], bufToRecv[index+2]);
					(*lbfor)[i][j][k].addCouplForceLoc(_addForce);
				}
			}
			
			bufToSend.resize(0);
			bufToRecv.resize(0);
		}

/*******************************************************************************************/
		
		/* COPY DEN AND J FROM A REAL REGION TO HALO NODES */
		void LatticeBoltzmann::copyDenMomToHalo() {
			int i, j, k, index;											// running indices and index of the node to be copied
			int numPopTransf = 4;										// number of populations or hydrod. moments to transfer
			int numDataTransf;											// number of data to transfer
			int rnode, snode;												// rank of the node to receive from and to send to
			int _numVels = getNumVels();
			std::vector<real> bufToSend, bufToRecv;	// buffers used to send and to receive the data

			int _offset = getHaloSkin();
			Int3D _myNi = getMyNi();
			Int3D _myPosition = getMyPosition();
			
			mpi::environment env;
			mpi::communicator world;
			
			//////////////////////
			//// X-direction /////
			//////////////////////
			numDataTransf = numPopTransf * _myNi[1] * _myNi[2];
			bufToSend.resize(numDataTransf);				// resize bufToSend
			bufToRecv.resize(numDataTransf);				// resize bufToRecv
			
			/* send to right, recv from left */
			snode = getMyNeighbour(1);
			rnode = getMyNeighbour(0);
			
			// prepare message for sending
			i = _myNi[0] - 2*_offset;
			for (k=0; k<_myNi[2]; k++) {
				for (j=0; j<_myNi[1]; j++) {
					index = numPopTransf*_myNi[1]*k + j*numPopTransf;
					for (int l = 0; l<numPopTransf; ++l) {
						bufToSend[index+l] = (*lbmom)[i][j][k].getMom_i(l);
					}
				}
			}
			
			// send and receive data or use memcpy if number of CPU in x-dir is 1
			if (getNodeGrid().getItem(0) > 1) {
				if (_myPosition[0] % 2 == 0) {
					world.send(snode, COMM_DEN_0, bufToSend);
					world.recv(rnode, COMM_DEN_0, bufToRecv);
				} else {
					world.recv(rnode, COMM_DEN_0, bufToRecv);
					world.send(snode, COMM_DEN_0, bufToSend);
				}
			} else {
				bufToRecv = bufToSend;
			}
			
			// unpack message
			i = 0;
			for (k=0; k<_myNi[2]; k++) {
				for (j=0; j<_myNi[1]; j++) {
					index = numPopTransf*_myNi[1]*k + j*numPopTransf;
					for (int l = 0; l<numPopTransf; ++l) {
						(*lbmom)[i][j][k].setMom_i(l, bufToRecv[index+l]);
					}
				}
			}
			
			/* send to left, recv from right i = 2, 8, 10, 12, 14 */
			snode = getMyNeighbour(0);
			rnode = getMyNeighbour(1);
			
			// prepare message for sending
			i = _offset;
			for (k=0; k<_myNi[2]; k++) {
				for (j=0; j<_myNi[1]; j++) {
					index = numPopTransf*_myNi[1]*k + j*numPopTransf;
					for (int l = 0; l<numPopTransf; ++l) {
						bufToSend[index+l] = (*lbmom)[i][j][k].getMom_i(l);
					}
				}
			}
			
			// send and receive data or use memcpy if number of CPU in x-dir is 1
			if (getNodeGrid().getItem(0) > 1) {
				if (_myPosition[0] % 2 == 0) {
					world.send(snode, COMM_DEN_1, bufToSend);
					world.recv(rnode, COMM_DEN_1, bufToRecv);
				} else {
					world.recv(rnode, COMM_DEN_1, bufToRecv);
					world.send(snode, COMM_DEN_1, bufToSend);
				}
			} else {
				bufToRecv = bufToSend;
			}
			
			// unpack message
			i = _myNi[0]-_offset;
			for (k=0; k<_myNi[2]; k++) {
				for (j=0; j<_myNi[1]; j++) {
					index = numPopTransf*_myNi[1]*k + j*numPopTransf;
					for (int l = 0; l<numPopTransf; ++l) {
						(*lbmom)[i][j][k].setMom_i(l, bufToRecv[index+l]);
					}
				}
			}
			
			//////////////////////
			//// Y-direction /////
			//////////////////////
			numDataTransf = numPopTransf * _myNi[0] * _myNi[2];
			bufToSend.resize(numDataTransf);				// resize bufToSend
			bufToRecv.resize(numDataTransf);				// resize bufToRecv
			
			/* send to right, recv from left */
			snode = getMyNeighbour(3);
			rnode = getMyNeighbour(2);
			
			// prepare message for sending
			j = _myNi[1] - 2*_offset;
			for (k=0; k<_myNi[2]; k++) {
				for (i=0; i<_myNi[0]; i++) {
					index = numPopTransf*_myNi[0]*k + i*numPopTransf;
					for (int l = 0; l<numPopTransf; ++l) {
						bufToSend[index+l] = (*lbmom)[i][j][k].getMom_i(l);
					}
				}
			}
			
			// send and receive data or use memcpy if number of CPU in y-dir is 1
			if (getNodeGrid().getItem(1) > 1) {
				if (_myPosition[1] % 2 == 0) {
					world.send(snode, COMM_DEN_2, bufToSend);
					world.recv(rnode, COMM_DEN_2, bufToRecv);
				} else {
					world.recv(rnode, COMM_DEN_2, bufToRecv);
					world.send(snode, COMM_DEN_2, bufToSend);
				}
			} else {
				bufToRecv = bufToSend;
			}
			
			// unpack message
			j = 0;
			for (k=0; k<_myNi[2]; k++) {
				for (i=0; i<_myNi[0]; i++) {
					index = numPopTransf*_myNi[0]*k + i*numPopTransf;
					for (int l = 0; l<numPopTransf; ++l) {
						(*lbmom)[i][j][k].setMom_i(l, bufToRecv[index+l]);
					}
				}
			}
			
			/* send to left, recv from right i = 4, 8, 9, 16, 18 */
			snode = getMyNeighbour(2);
			rnode = getMyNeighbour(3);
			
			// prepare message for sending
			j = _offset;
			for (k=0; k<_myNi[2]; k++) {
				for (i=0; i<_myNi[0]; i++) {
					index = numPopTransf*_myNi[0]*k + i*numPopTransf;
					for (int l = 0; l<numPopTransf; ++l) {
						bufToSend[index+l] = (*lbmom)[i][j][k].getMom_i(l);
					}
				}
			}
			
			// send and receive data or use memcpy if number of CPU in y-dir is 1
			if (getNodeGrid().getItem(1) > 1) {
				if (_myPosition[1] % 2 == 0) {
					world.send(snode, COMM_DEN_3, bufToSend);
					world.recv(rnode, COMM_DEN_3, bufToRecv);
				} else {
					world.recv(rnode, COMM_DEN_3, bufToRecv);
					world.send(snode, COMM_DEN_3, bufToSend);
				}
			} else {
				bufToRecv = bufToSend;
			}
			
			// unpack message
			j = _myNi[1] - _offset;
			for (k=0; k<_myNi[2]; k++) {
				for (i=0; i<_myNi[0]; i++) {
					index = numPopTransf*_myNi[0]*k + i*numPopTransf;
					for (int l = 0; l<numPopTransf; ++l) {
						(*lbmom)[i][j][k].setMom_i(l, bufToRecv[index+l]);
					}
				}
			}
			
			//////////////////////
			//// Z-direction /////
			//////////////////////
			numDataTransf = numPopTransf * _myNi[0] * _myNi[1];
			bufToSend.resize(numDataTransf);				// resize bufToSend
			bufToRecv.resize(numDataTransf);				// resize bufToRecv
			
			/* send to right, recv from left i = 5, 11, 14, 15, 18 */
			snode = getMyNeighbour(5);
			rnode = getMyNeighbour(4);
			
			// prepare message for sending
			k = _myNi[2] - 2 * _offset;
			for (j=0; j<_myNi[1]; j++) {
				for (i=0; i<_myNi[0]; i++) {
					index = numPopTransf*_myNi[0]*j + i*numPopTransf;
					for (int l = 0; l<numPopTransf; ++l) {
						bufToSend[index+l] = (*lbmom)[i][j][k].getMom_i(l);
					}
				}
			}
			
			// send and receive data or use memcpy if number of CPU in z-dir is 1
			if (getNodeGrid().getItem(2) > 1) {
				if (_myPosition[2] % 2 == 0) {
					world.send(snode, COMM_DEN_4, bufToSend);
					world.recv(rnode, COMM_DEN_4, bufToRecv);
				} else {
					world.recv(rnode, COMM_DEN_4, bufToRecv);
					world.send(snode, COMM_DEN_4, bufToSend);
				}
			} else {
				bufToRecv = bufToSend;
			}
			
			// unpack message
			k = 0;
			for (j=0; j<_myNi[1]; j++) {
				for (i=0; i<_myNi[0]; i++) {
					index = numPopTransf*_myNi[0]*j + i*numPopTransf;
					for (int l = 0; l<numPopTransf; ++l) {
						(*lbmom)[i][j][k].setMom_i(l, bufToRecv[index+l]);
					}
				}
			}
			
			/* send to left, recv from right i = 6, 12, 13, 16, 17 */
			snode = getMyNeighbour(4);
			rnode = getMyNeighbour(5);
			
			// prepare message for sending
			k = _offset;
			for (j=0; j<_myNi[1]; j++) {
				for (i=0; i<_myNi[0]; i++) {
					index = numPopTransf*_myNi[0]*j + i*numPopTransf;
					for (int l = 0; l<numPopTransf; ++l) {
						bufToSend[index+l] = (*lbmom)[i][j][k].getMom_i(l);
					}
				}
			}
			
			// send and receive data or use memcpy if number of CPU in z-dir is 1
			if (getNodeGrid().getItem(2) > 1) {
				if (_myPosition[2] % 2 == 0) {
					world.send(snode, COMM_DEN_5, bufToSend);
					world.recv(rnode, COMM_DEN_5, bufToRecv);
				} else {
					world.recv(rnode, COMM_DEN_5, bufToRecv);
					world.send(snode, COMM_DEN_5, bufToSend);
				}
			} else {
				bufToRecv = bufToSend;
			}
			
			// unpack message
			k = _myNi[2] - _offset;
			for (j=0; j<_myNi[1]; j++) {
				for (i=0; i<_myNi[0]; i++) {
					index = numPopTransf*_myNi[0]*j + i*numPopTransf;
					for (int l = 0; l<numPopTransf; ++l) {
						(*lbmom)[i][j][k].setMom_i(l, bufToRecv[index+l]);
					}
				}
			}
			
			bufToSend.resize(0);
			bufToRecv.resize(0);
		}
		
/*******************************************************************************************/

    /* Destructor of the LB */
    LatticeBoltzmann::~LatticeBoltzmann() {
			disconnect();
    }
		
/*******************************************************************************************/
		
		/******************************
		 ** REGISTRATION WITH PYTHON **
		 ******************************/
		
/*******************************************************************************************/
		
		void LatticeBoltzmann::registerPython() {

			using namespace espressopp::python;
			class_<LatticeBoltzmann, shared_ptr<LatticeBoltzmann>, bases<Extension> >

			("integrator_LatticeBoltzmann", init<	shared_ptr< System >,
																					Int3D, Int3D, real, real, int, int >())
			.add_property("nodeGrid", &LatticeBoltzmann::getNodeGrid, &LatticeBoltzmann::setNodeGrid)
			.add_property("Ni", &LatticeBoltzmann::getNi, &LatticeBoltzmann::setNi)
			.add_property("a", &LatticeBoltzmann::getA, &LatticeBoltzmann::setA)
			.add_property("tau", &LatticeBoltzmann::getTau, &LatticeBoltzmann::setTau)
			.add_property("numDims", &LatticeBoltzmann::getNumDims, &LatticeBoltzmann::setNumDims)
			.add_property("numVels", &LatticeBoltzmann::getNumVels, &LatticeBoltzmann::setNumVels)
			.add_property("visc_b", &LatticeBoltzmann::getViscB, &LatticeBoltzmann::setViscB)
			.add_property("visc_s", &LatticeBoltzmann::getViscS, &LatticeBoltzmann::setViscS)
			.add_property("gamma_b", &LatticeBoltzmann::getGammaB, &LatticeBoltzmann::setGammaB)
			.add_property("gamma_s", &LatticeBoltzmann::getGammaS, &LatticeBoltzmann::setGammaS)
			.add_property("gamma_odd", &LatticeBoltzmann::getGammaOdd, &LatticeBoltzmann::setGammaOdd)
			.add_property("gamma_even", &LatticeBoltzmann::getGammaEven, &LatticeBoltzmann::setGammaEven)
			.add_property("lbTemp", &LatticeBoltzmann::getLBTemp, &LatticeBoltzmann::setLBTemp)
			.add_property("fricCoeff", &LatticeBoltzmann::getFricCoeff, &LatticeBoltzmann::setFricCoeff)
			.add_property("nSteps", &LatticeBoltzmann::getNSteps, &LatticeBoltzmann::setNSteps)
			.add_property("profStep", &LatticeBoltzmann::getProfStep, &LatticeBoltzmann::setProfStep)
			.def("readCouplForces", &LatticeBoltzmann::readCouplForces)
			.def("saveCouplForces", &LatticeBoltzmann::saveCouplForces)
			.def("connect", &LatticeBoltzmann::connect)
			.def("disconnect", &LatticeBoltzmann::disconnect)
			;
    }
  }
}
