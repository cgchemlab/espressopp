#ifndef _STORAGE_DOMAINDECOMPOSITION_HPP
#define _STORAGE_DOMAINDECOMPOSITION_HPP
#include "Storage.hpp"
#include "CellGrid.hpp"
#include "NodeGrid.hpp"

namespace espresso {
  namespace storage {
    class NodeGridMismatch: public std::runtime_error
    {
    public:
      NodeGridMismatch();
    };

    class DomainDecomposition: public Storage {
    public:
      DomainDecomposition(shared_ptr< System >,
			  const boost::mpi::communicator &,
			  const int nodeGrid[3],
			  const int cellGrid[3],
			  bool useVList);

      virtual ~DomainDecomposition() {}

      virtual Cell *mapPositionToCellClipped(const real pos[3]);
      virtual Cell *mapPositionToCellChecked(const real pos[3]);

      const NodeGrid &getNodeGrid() const { return nodeGrid; }
      const CellGrid &getCellGrid() const { return cellGrid; }

    protected:
      virtual void resortRealParticles();
      virtual void exchangeGhosts();
      virtual void updateGhosts();
      virtual void collectGhostForces();

      void doGhostCommunication(bool sizesFirst,
				bool realToGhosts,
				const int dataElements = 0);

      void prepareGhostCommunication();

      /// init global Verlet list
      void initCellInteractions();
      /// set the grids and allocate space accordingly
      void createCellGrid(const int nodeGrid[3], const int cellGrid[3]);
      /// sort cells into local/ghost cell arrays
      void markCells();
      /// fill a list of cells with the cells from a certain region of the domain grid
      void fillCells(std::vector<Cell *> &,
		     const int leftBoundary[3],
		     const int rightBoundary[3]);
      /** read particles from a temporary buffer into the local cell structure.
	  The direction determines in which direction to fold the particles position.
	  Returns true if one of the given particles did not belong to this processors
	  domain.
      */
      bool appendParticles(ParticleList &, int dir);

      /// spatial domain decomposition of nodes
      NodeGrid nodeGrid;

      /// spatial domain decomposition on node in cells
      CellGrid cellGrid;

      /// expected capacity of send/recv buffers for neighbor communication
      size_t exchangeBufferSize;

      /** which cells to send and receive during one communication step.
	  In case this is a communication with ourselves, the send-cells
	  are transferred to the recv-cells. */
      struct CommCells {
	std::vector<Cell *> reals;
	std::vector<Cell *> ghosts;
      };
      /** which cells to send left, right, up, down, ...
	  For the order, see NodeGrid.
      */
      CommCells commCells[6];

      static LOG4ESPP_DECL_LOGGER(logger);
    };
  }
}
#endif
