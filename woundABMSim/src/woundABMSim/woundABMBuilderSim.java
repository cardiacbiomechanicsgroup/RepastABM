package woundABMSim;

import repast.simphony.context.Context;
import repast.simphony.dataLoader.ContextBuilder;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.parameter.Parameters;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridPoint;

public class woundABMBuilderSim implements ContextBuilder<Object>{
	
	@Override
/* Constructor */
	public Context<Object> build(Context<Object> context) {
		
		// Reload static values based off of parameters in CellAgentSim (necessary for proper batch execution)
		CellAgentSim.InitializeStaticValues();
		
		// Pull parameters
		Parameters p = RunEnvironment.getInstance().getParameters();
	
		// Cell seeding parameters
		int initialCellCount = (Integer) p.getValue("initialCellCount");	// cells (1cell/400um^2)
		int endTick = (Integer) p.getValue("endTick");						// hr
		
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