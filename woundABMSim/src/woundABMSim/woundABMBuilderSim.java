package woundABMSim;

import repast.simphony.context.Context;
import repast.simphony.dataLoader.ContextBuilder;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridPoint;

public class woundABMBuilderSim implements ContextBuilder<Object>{

	@Override
/* Constructor */
	public Context<Object> build(Context<Object> context) {
		
		// Cell seeding parameters
		int initialCellCount = 576;				// cells (1cell/400um^2)
		int endTick = 1009;						// hr
		
		// Build subcontext "woundABMContext"
		woundABMContextSim woundabmspace = new woundABMContextSim();
		context.addSubContext(woundabmspace);
		context.add(woundabmspace);

		// Add initial cell population
		for (int i = 0; i < initialCellCount; i++) {
			CellAgentSim cellagent = new CellAgentSim();
			woundabmspace.add(cellagent);

			// Check for and delete overlapping cells
			Grid<?> grid = (Grid<?>) woundabmspace.getProjection("Cell Grid");
			GridPoint pt = grid.getLocation(cellagent);
			int x = pt.getX();
			int y = pt.getY();			
			if (cellagent.findNumNeighbor(grid, x, y) > 1) {
				woundabmspace.remove(cellagent);
				initialCellCount++;										// Maintain cell seeding density
			}
		}
		
		// End simulation
		RunEnvironment.getInstance().endAt(endTick);
		return context;
	}
}