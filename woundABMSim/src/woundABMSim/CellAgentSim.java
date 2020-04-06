package woundABMSim;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.parameter.Parameters;
import repast.simphony.random.RandomHelper;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridPoint;
import repast.simphony.util.ContextUtils;
import repast.simphony.valueLayer.GridValueLayer;

public class CellAgentSim {

	// Geometry
	private static double sampleWidth; 					// um
	private static double sampleHeight;					// um
	private static double gridUnitSize;					// um
	private static int gridWidth;						// grids
	private static int gridHeight;						// grids
	private static int gridDiameter;					// Cell diameter in grids
	private static int gridRadius;						// Cell radius in grids
	private static int cellGridArea;					// Cell area in grids
	private static boolean cellRemoval;
	private final static double cellRadius = 5.0;		// um
	private final static double woundRadius = 113.0;	// um
	
	// Agent Tracking
	private double mitosisTime;			// hr
	private double mitosisAge;			// hr
	private double apoptosisAge;		// hr
	private double cellSpeed;			// um/hr	
	public double angleSelection;		// deg
	private int adjCellSpeed;
	private double migrationDistance;	// um
	
	// Weighting factors
	private static double Wp;	// Persistance cue weight
	private static double Ws; 	// Structural cue weight
	private static double Wm;	// Mechanics cue weight
	private static double Wc;	// Chemokine cue weight

	// Migration
	private final static double alphaGeometric = woundRadius/2828.0;	// Cryo infarct radius = 2828 um;
	private final static double alphaInfiltration = woundRadius/500.0; 	// Cryo infarct thickness = 500 um;
	private final static double cellSpeedMin = alphaInfiltration*1.0;	// um/hr
	private final static double cellSpeedMax = alphaInfiltration*10.0;	// um/hr

	// Discretization
	private final static int binSize = 5;				// deg
	private final static int binNum360 = 360/binSize;	// bins
	private final static int binNum180 = 180/binSize;	// bins

	// Collagen
	private static boolean colRotation;
	private static String depoType;
	private final static double fiberDensity = 178.25;		// fibers/um^2/7um thickness from cryoinfarct data
	private final static double kColRotMax = 3.0;			// deg/hr
	private final static double kColRotMin = 0.3;			// deg/hr
	private final static double kColGenMax = 17.9874;		// fiber/cell/hr
	private final static double kColGenMin = 0.1799;		// fiber/cell/hr
	private final static double kColDegMax = 0.0038;   		// unitless scaling factor
	private final static double kColDegMin = 3.7676e-04;	// unitless scaling factor

	// Cell cycle timing
	private final static double timeStep = 0.5;	// hr
	private static double gMitosisTime;			// hr
	private static double apoptosisTime;		// hr
	private static double depositionTime;		// hr
	private static double degradationTime;		// hr
	
	// Persistence
	private final static double pMin = 0.25;     // min persistence time (hr)
	private final static double pNoCue = 1.2;    // persistence time in the absence of external directional cues
	private static double rhoTuning;

	// Chemokine
	private final static double concMin = 0.0;
	private final static double concMax = 1.0;
	private final static double kcdeg = 0.001;         		// 1/s
	private final static double Dc = 100;              		// um2/s
	private final static double L = 2828;              		// um
	private final static double lam = Math.sqrt(kcdeg/Dc);  // 1/um
	private final static double maxConcGrad = (lam*L*Math.sinh(lam*L))/(Math.cosh(lam*L)+Math.sinh(lam*L))/(alphaGeometric*L);    // 1/um

	// Mechanics
	private final static double eii = 0.05;	// Maximum strain anisotropy based on cryoinfarct data
	private final static double ejj = 0;
	private final static double eij = 0;
	
	// Guidance cue normalization factors based on maximum values
	private final static double Mc = 0.5*cellRadius*maxConcGrad;
	private final static double Mm = Math.sqrt(Math.pow(eii, 2)-2*eii*ejj+Math.pow(ejj, 2)+4*Math.pow(eij,2))/4; // Ms = 0.0125	 
	private final static double Mp = cellSpeedMax;
	private static boolean scaleStructuralCue;

	/* Constructors */
	public CellAgentSim() {
		
		// Assign parameters
		this.mitosisTime = gMitosisTime;
		this.mitosisAge = RandomHelper.nextDoubleFromTo(0, gMitosisTime);
		this.apoptosisAge = mitosisAge;
		this.angleSelection = RandomHelper.nextDoubleFromTo(-180, 180);
		this.adjCellSpeed = 1;
		this.migrationDistance = RandomHelper.nextDoubleFromTo(0, gridUnitSize);
		
		List<GridPoint> coveredSites = cellCoverage(0, 0);
		cellGridArea = coveredSites.size();
	}

	public CellAgentSim(double mitosisTime,	double apoptosisAge, double mitosisAge, double angleSelection, 
			int adjCellSpeed) {

		// Assign parameters
		this.mitosisTime = mitosisTime;
		this.mitosisAge = mitosisAge;
		this.apoptosisAge = apoptosisAge;
		this.angleSelection = angleSelection;
		this.adjCellSpeed = adjCellSpeed;
		this.migrationDistance = 0;
		
		List<GridPoint> coveredSites = cellCoverage(0, 0);
		cellGridArea = coveredSites.size();
	}
	
	public static void InitializeStaticValues() {
		
		// Load parameters from GUI
		Parameters p = RunEnvironment.getInstance().getParameters();
		sampleWidth = (Double) p.getValue("sampleWidth");			// um
		sampleHeight = (Double) p.getValue("sampleHeight");			// um
		gridUnitSize = (Double) p.getValue("gridUnitSize");			// um
		Wp = (Double) p.getValue("Wp");								// Persistance cue weight
		Ws = (Double) p.getValue("Ws"); 							// Structural cue weight
		Wm = (Double) p.getValue("Wm");								// Mechanics cue weight
		Wc = (Double) p.getValue("Wc");								// Chemokine cue weight
		scaleStructuralCue = (Boolean) p.getValue("scaleStructuralCue");
		cellRemoval = (Boolean) p.getValue("cellRemoval");
		colRotation = (Boolean) p.getValue("colRotation");
		depoType = (String) p.getValue("depoType"); 				// "Aligned" or "Random"
		gMitosisTime = (Double) p.getValue("mitosisTime");			// hr
		apoptosisTime = (Double) p.getValue("apoptosisTime");		// hr
		depositionTime = (Double) p.getValue("depositionTime");		// hr
		degradationTime = (Double) p.getValue("degradationTime");	// hr
		
		// Parameter dependent values
		gridWidth = (int) Math.floor(sampleWidth/gridUnitSize);	// grids
		gridHeight = (int) Math.floor(sampleHeight/gridUnitSize);	// grids
		gridDiameter = (int) Math.round(cellRadius*2/gridUnitSize);	// Cell diameter in grids
		gridRadius = (int) Math.ceil(cellRadius/gridUnitSize);	// Cell radius in grids
		rhoTuning = 2*Wp*pMin/(pNoCue-pMin); // rhoTuning = 0.1754
	}

	/* Primary scheduled methods in order of execution */
	// Infarction method
	@ScheduledMethod(start = 0, priority = 3)
 	public void cellRemoval() {
		
		// Remove cells in the wound area
		if (cellRemoval == true) {
			final woundABMContextSim woundabmspace = (woundABMContextSim) ContextUtils.getContext(this);
			final Grid<?> grid = (Grid<?>) woundabmspace.getProjection("Cell Grid");
			final GridPoint pt = grid.getLocation(this);	

			final double removalDist = (woundRadius/gridUnitSize);
			final double dist = Math.sqrt(Math.pow((pt.getX()-gridWidth/2), 2)+
					Math.pow((pt.getY()-gridHeight/2), 2));
			if (dist <= removalDist) {
				woundabmspace.remove(this);
			}
		}
	}

	// Agent decision tree
	@ScheduledMethod(start = 0.5, interval = 0.5, priority = 1)
	public void step() {
		final woundABMContextSim woundABMSpace = (woundABMContextSim) ContextUtils.getContext(this);
		final Grid<?> grid = (Grid<?>) woundABMSpace.getProjection("Cell Grid");
		final GridValueLayer chemokine = (GridValueLayer) woundABMSpace.getValueLayer("Chemokine");
		final GridPoint pt = grid.getLocation(this);
		final double effectiveChemokine = (chemokine.get(pt.getX(), pt.getY())-concMin)/(concMax-concMin);

		// Define chemokine dependent mitosis to ensure cells are always coming from the outside of wound
		this.mitosisTime = effectiveChemokine*(gMitosisTime-apoptosisTime)+apoptosisTime;

		// Define chemokine dependent migration to speed up wound closure
		this.cellSpeed = effectiveChemokine*(cellSpeedMax-cellSpeedMin)+cellSpeedMin;

		// Get location
		final int x = pt.getX();
		final int y = pt.getY();

		// Check for mitosis
		if (mitosisAge >= mitosisTime) {	
			mitosisAge = randomClock();	// Quiescent cell sets mitosis age to random number
			mitose(woundABMSpace, grid, x, y);

		} else {

			// Update cell orientation
			this.angleSelection = guidanceCue(woundABMSpace, x, y);

			// Get grids covered by the cell
			List<GridPoint> coveredSites = cellCoverage(x, y);

			// Deposit, degrade, and rotate collagen in the grids beneath the cell
			for (GridPoint site : coveredSites) {

				// Get site x and y coordinates
				final int siteX = site.getX();
				final int siteY = site.getY();

				// Retrieve structural value layers
				final ArrayList<GridValueLayer> collagenLayers = (ArrayList<GridValueLayer>)woundABMSpace.GridValueLayerList();
				final ArrayList<GridValueLayer> fibrinLayers = (ArrayList<GridValueLayer>)woundABMSpace.fibrinGridValueLayerList();	
				
				// Check for rotation, degradation, deposition
				if (colRotation == true) {
					collagenRotation(siteX, siteY, chemokine, collagenLayers);
				}
				matrixDegradation(siteX, siteY, chemokine, collagenLayers, fibrinLayers, degradationTime);
				collagenDeposition(siteX, siteY, chemokine, collagenLayers, depositionTime);
			}
			
			// Check for migration
			if (migrationDistance >= gridUnitSize) {
				final int grids = (int) Math.floor(migrationDistance/gridUnitSize);
				this.migrationDistance = migrationDistance-gridUnitSize*grids;
				moveTowards(grid, x, y, grids);
			}
		}

		// Check for apoptosis
		if (apoptosisAge >= apoptosisTime) {
			woundABMSpace.remove(this);
		}

		// Increment the counters
		apoptosisAge = apoptosisAge+timeStep;
		mitosisAge = mitosisAge+timeStep;
		migrationDistance = migrationDistance+cellSpeed*timeStep;	
	}

	/* Secondary methods */
	// Mitosis method
	@SuppressWarnings({ "rawtypes", "unchecked" })
	private boolean mitose(woundABMContextSim woundabmspace, Grid grid, int x, int y) {

		// Create list of vacant neighboring sites
		final List<GridPoint> emptySites = getEmptySites(grid, x, y, gridDiameter);
		if (emptySites.size() > 0) {

			// Select a random site from list of empty sites
			final GridPoint site = emptySites.get(RandomHelper.nextIntFromTo(0,(emptySites.size()-1)));
			final int xCoor = site.getX();
			final int yCoor = site.getY();

			// Add daughter (angle same as mother's)
			final double dAngle = this.angleSelection;
			final double dMitosisClock = randomClock();
			CellAgentSim cellagent = new CellAgentSim(mitosisTime, dMitosisClock, dMitosisClock, dAngle, 1);
			woundabmspace.add(cellagent);
			grid.moveTo(cellagent, xCoor, yCoor);
			
			return true;
			
		} else {
			
			return false;
		}
	}

	// Returns an angle based on guidance cues
	private double guidanceCue(woundABMContextSim woundabmspace, int x, int y) {

		// Calculate structural guidance cue normalization factors based on maximum values
		double Ms = 1;
		if (scaleStructuralCue == true) {
			final GridValueLayer colTotal = (GridValueLayer) woundabmspace.getValueLayer("Collagen Sum");
			final List<GridPoint> coveredSites = cellCoverage(x, y);
			double colSum = 0;
			for (GridPoint site : coveredSites ) {
				final int siteX = site.getX();
				final int siteY = site.getY();
				colSum = colSum+colTotal.get(siteX, siteY);
			}
			Ms = (colSum/(cellGridArea*fiberDensity))/0.3;
		}

		// Acquire chemokine, strain, fiber, and persistence cues
		final double[] chemoCue = getChemoCue(woundabmspace, x, y);
		double CM = chemoCue[0]; 
		double CX = chemoCue[1];
		double CY = chemoCue[2];
		final double[] strainCue = getStrainCue(woundabmspace, x, y); 
		double MM = strainCue[0]; 
		double MX = strainCue[1];
		double MY = strainCue[2];
		final double[] fiberCue = getColFiberCue(woundabmspace, x, y);
		double SM = fiberCue[0]; 
		double SX = fiberCue[1];
		double SY = fiberCue[2];	
		final double[] persistence = getPersistence(woundabmspace, x, y);
		double PM = persistence[0];
		double PX = persistence[1]; 
		double PY = persistence[2];
		
		// Calculate guidance cues strength using maximum value normalization and weighting factors
		final double Sc = (Wc*CM)/Mc;
		final double Sm = (Wm*MM)/Mm;
		final double Ss = (Ws*SM)/Ms;
		final double Sp = (Wp*PM)/Mp;

		// Calculate the resultant guidance cues from the weighted guidance cues unit vectors
		double[] RX = new double[4];
		double[] RY = new double[4];
		double[] RM = new double[4];
		CX = Sc*CX;
		CY = Sc*CY;
		MX = Sm*MX;
		MY = Sm*MY;
		SX = Ss*SX;
		SY = Ss*SY;
		PX = Sp*PX;
		PY = Sp*PY;
		RX[0] = CX+PX+SX+MX;
		RY[0] = CY+PY+SY+MY;
		RM[0] = Math.sqrt(Math.pow(RX[0], 2)+Math.pow(RY[0], 2));
		RX[1] = CX+PX-SX+MX;
		RY[1] = CY+PY-SY+MY;
		RM[1] = Math.sqrt(Math.pow(RX[1], 2)+Math.pow(RY[1], 2));
		RX[2] = CX+PX+SX-MX;
		RY[2] = CY+PY+SY-MY;
		RM[2] = Math.sqrt(Math.pow(RX[2], 2)+Math.pow(RY[2], 2));
		RX[3] = CX+PX-SX-MX;
		RY[3] = CY+PY-SY-MY;
		RM[3] = Math.sqrt(Math.pow(RX[3], 2)+Math.pow(RY[3], 2));
		
		// Find maximum resultant magnitude/s
		double maxRM = 0;
		ArrayList<Integer> RMlist = new ArrayList<Integer>();
		int maxResultantIndex = 99;
		for (int i = 0; i < 4; i++) {
			if (RM[i] > maxRM) {
				maxRM = RM[i];
				maxResultantIndex = i;
			} else if (RM[i] == maxRM) {
				if (RMlist.size() < 1) {
					RMlist.add(maxResultantIndex);
				}
				RMlist.add(i);
			}
		}

		// If two resultant magnitudes are equal, pick randomly out of two
		if (RMlist.size() < 1) {
			RMlist.add(maxResultantIndex);
		}
		int maxRMindex = RMlist.get(RandomHelper.nextIntFromTo(0,(RMlist.size()-1)));
		
		// Calculate the resultant angle "theta" and normalized resultant magnitude "rho"
		final double theta = Math.atan2(RY[maxRMindex], RX[maxRMindex]);
		final double rho = RM[maxRMindex]/(rhoTuning+Ss+Sm+Sc+Sp);

		// Select a direction from the probability distribution, exceptions for rho >=0.999 and rho <=0.001
		double[] binDir360 = new double[72];
		for (int i = 0; i < binNum360; i++) {
			binDir360[i] = binSize*i-177.5;
		}
		double guidedAngle; 
		if (rho >= 0.999) {
			guidedAngle = Math.ceil(Math.toDegrees(theta));
		} else if (rho <= 0.001) {
			guidedAngle = binDir360[(int) RandomHelper.nextIntFromTo(0,(binNum360-1))];
		} else {
			Double[] angleSet = getAngleSet(theta, rho);
			if (angleSet.length > 0) {
				int randomAngIdx = (int) RandomHelper.nextIntFromTo(0,angleSet.length-1);
				guidedAngle = angleSet[randomAngIdx];
			} else {
				guidedAngle = binDir360[RandomHelper.nextIntFromTo(0,binDir360.length-1)];
			}
		}
		
		return guidedAngle;
	}

	// Migration method
	@SuppressWarnings({ "rawtypes", "unchecked" })
	private void moveTowards(Grid grid, int x, int y, int extent) {

		// Search for a migration site starting from farthest to nearest
		while (extent > 0) {

			// Get list of sites within extent and the angle window that are available for migration
			double angle = this.angleSelection;
			if (angle < 0) {
				angle = angle+360.0;
			}
			final List<GridPoint> emptySites = getMigrationSite(grid, x, y, extent, angle);
			if (!emptySites.isEmpty()) {

				// Select a random migration site from the list of those available to move to
				final int migrationSiteIdx = RandomHelper.nextIntFromTo(0,(emptySites.size()-1));
				final GridPoint migrationSite = emptySites.get(migrationSiteIdx);
				final int desiredX = migrationSite.getX();
				final int desiredY = migrationSite.getY();

				// Move to new location
				grid.moveTo(this, desiredX, desiredY);
				this.adjCellSpeed = 1;
				extent = 0;

			}// Abort move and reset migration distance
			else if (extent == 1) {
				this.adjCellSpeed = 0;
				this.migrationDistance = 0;
			}

			// Reduce search distance
			extent = extent-1;
		}
	}	

	// Returns a list containing all the grids who's center is within the cell boarder
	private List<GridPoint> cellCoverage(int x, int y) {

		// Iterate through neighbors
		final List<GridPoint> coveredSites = new ArrayList<GridPoint>();
		for (int i = -gridRadius; i <= gridRadius; i++) {
			int xCoor = x+i;
			for (int j = -gridRadius; j <= gridRadius; j++) {
				int yCoor = y+j;

				// Check radial distance
				final double dist = Math.sqrt(Math.pow(xCoor-x,2)+Math.pow(yCoor-y,2))*gridUnitSize;
				if (dist <= cellRadius) {

					// Wrapped space
					if (xCoor < 0) {
						xCoor = xCoor+gridWidth;
					} else if (xCoor >= gridWidth) {
						xCoor = xCoor-gridWidth;
					}
					if (yCoor < 0) {
						yCoor = yCoor+gridHeight;
					} else if (yCoor >= gridHeight) {
						yCoor = yCoor-gridHeight;
					}

					coveredSites.add(new GridPoint(xCoor,yCoor));
				}
			}
		}
		
		return coveredSites;
	}

	// Rotate collagen
	private void collagenRotation(int x, int y, GridValueLayer chemokine, ArrayList<GridValueLayer> gridValueLayerList) {

		// Calculate collagen fiber rotation constant based on chemokine concentration
		final double kColRot = ((chemokine.get(x,y)-concMin)/(concMax-concMin)*(kColRotMax-kColRotMin))+kColRotMin;
		final double Q = Math.ceil(kColRot*timeStep/2);

		// Get all 36 bins of collagen
		double [] colContent = new double[binNum180];
		ArrayList <Double> tempColContent = new ArrayList<Double>();
		for (int j = 0; j < binNum180; j++) {
			colContent[j] = gridValueLayerList.get(j).get(x, y);
			tempColContent.add((double) 0.0);
		}

		// Converts angleSelection range of [-180 180] to range of [-90 90]
		final double radAngle = Math.toRadians(this.angleSelection);
		final double cellAngle = Math.toDegrees(Math.atan(Math.sin(radAngle)/Math.cos(radAngle)));
		for (int j = 0; j < binNum180; j++) {	

			// Get fiber offset from cell angle in radians
			double fiberBin = binSize*j-87.5;
			double newAngle = 360;
			final double fiberOffset = Math.toRadians(cellAngle-fiberBin);
			final double tanFiberOff = Math.tan(fiberOffset);
			final double sinFiberOff = Math.abs(Math.sin(fiberOffset));
			if ((tanFiberOff > 0) && (sinFiberOff != 1)) {// fiberOffset < 90
				newAngle = fiberBin+(kColRot*timeStep/Q*(sinFiberOff));
			} else if ((tanFiberOff <= 0) && (sinFiberOff != 1)) {// fiberOffset > 90
				newAngle = fiberBin-(kColRot*timeStep/Q*(sinFiberOff));
			} else if (Math.abs(Math.sin(fiberOffset)) == 1) {// fiberOffset ~= 90
				double flip = 1;
				final int dice = RandomHelper.nextIntFromTo(0,1);
				if (dice == 1) {
					flip = -1;
				}
				newAngle = fiberBin+(flip*kColRot*timeStep/Q*(sinFiberOff));
			}

			// Detect and correct angles outside of fiber bin range
			if (newAngle <= -90) {
				newAngle = newAngle+180;
			} else if (newAngle > 90) {
				newAngle = newAngle-180;
			}
			
			// Assign newAngle to a fiber bin
			final int newIndex = (int) Math.ceil(newAngle/binSize)+17;
			final double rotBin = binSize*newIndex-87.5;
			
			// Rotate collagen
			int rotIndex = 99;
			final double rotFrac = (newAngle-rotBin)/binSize;
			double colKept;
			double colRotated;
			if (rotFrac >= 0) { //kept
				colKept = tempColContent.get(newIndex)+colContent[j]*(1-rotFrac);
				if (newIndex == 35) {// rotated
					rotIndex = 0;
				} else {
					rotIndex = newIndex+1;
				}
				colRotated = tempColContent.get(rotIndex)+colContent[j]*rotFrac;
			} else { // kept
				colKept = tempColContent.get(newIndex)+colContent[j]*(1+rotFrac);
				if (newIndex == 0) {
					rotIndex = 35;
				} else {
					rotIndex = newIndex-1;
				}
				colRotated = tempColContent.get(rotIndex)-colContent[j]*rotFrac;
			}
			tempColContent.set(newIndex, colKept);
			tempColContent.set(rotIndex, colRotated);	
		}

		// Update all 36 collagen bins
		for (int i = 0; i < binNum180; i++) {
			gridValueLayerList.get(i).set(tempColContent.get(i), x, y);
		}
	}

	// Fiber degradation method
	private void matrixDegradation(int x, int y, GridValueLayer chemokine, ArrayList<GridValueLayer> gridValueLayerList, ArrayList<GridValueLayer> fibringridValueLayerList, double degrationInterval) {

		// Calculate fiber degradation constants with respect to chemokine concentration
		final double kColDeg = ((chemokine.get(x, y)-concMin)/(concMax-concMin)*(kColDegMax-kColDegMin)+kColDegMin);
		final double kFibDeg = kColDeg;

		// Degrade fibers uniformly across all bins
		for (int i = 0; i < binNum180; i++) {

			// Collagen degradation 
			final double colContent = gridValueLayerList.get(i).get(x, y);
			final double newColContent = colContent*(1-kColDeg*degrationInterval);
			gridValueLayerList.get(i).set(newColContent, x, y);

			// Fibrin degradation (only within wound)
			final double distToCenter = Math.sqrt(Math.pow((x-gridWidth/2), 2)+Math.pow((y-gridHeight/2), 2));
			final double woundSize = (double) woundRadius/gridUnitSize;
			if (distToCenter < woundSize) {
				final double fibrinContent = fibringridValueLayerList.get(i).get(x, y);
				final double newFibrinContent = fibrinContent*(1-kFibDeg*degrationInterval);
				fibringridValueLayerList.get(i).set(newFibrinContent, x, y);
			}
		}
	}

	// Aligned collagen deposition method
	private void collagenDeposition(int x, int y, GridValueLayer chemokine, ArrayList<GridValueLayer> gridValueLayerList, double depositionInterval) {

		// Calculate collagen fiber generation rate based on chemokine concentration
		final double kColGen = ((chemokine.get(x, y))-concMin)/(concMax-concMin)*(kColGenMax-kColGenMin)+kColGenMin;

		// Set collagen deposition angle based on deposition type
		double depoAngle;
		final double cellAngle = this.angleSelection;
		if (depoType.equals("Aligned")) {
			
			if (cellAngle <= -90.0) {
				depoAngle = cellAngle+180.0;
			} else if (cellAngle > 90.0) {
				depoAngle = cellAngle-180.0;
			} else {
				depoAngle = cellAngle;
			}
			
		} else {
			depoAngle = RandomHelper.nextDoubleFromTo(-90, 90);
		}

		// Assign to the closest bin
		final int bin = (int) Math.ceil(depoAngle/binSize)+17;

		// Deposit the aligned collagen
		final double deltaCol = depositionInterval*kColGen/cellGridArea;
		final double newColContent = gridValueLayerList.get(bin).get(x, y)+(deltaCol);
		gridValueLayerList.get(bin).set(newColContent, x, y);
	}

	// Returns a random clock time
	private double randomClock() {
		final Random r = new Random();
		final double time = Math.round(0.1*mitosisTime*r.nextGaussian()/timeStep)*timeStep;
		return time;
	}

	/* Tertiary methods */
	// Returns a vector of chemokine cues
	private double[] getChemoCue(woundABMContextSim woundabmspace, int x, int y) {

		// Load chemokine cue value layers
		final GridValueLayer CM = (GridValueLayer) woundabmspace.getValueLayer("CM");
		final GridValueLayer CX = (GridValueLayer) woundabmspace.getValueLayer("CX");
		final GridValueLayer CY = (GridValueLayer) woundabmspace.getValueLayer("CY");

		// Get chemokine cue at cell location
		final double CMabm = CM.get(x, y);
		final double CXabm = CX.get(x, y); 
		final double CYabm = CY.get(x, y);

		return new double[] {CMabm, CXabm, CYabm};
	}

	// Returns a vector of strain cues
	private double[] getStrainCue(woundABMContextSim woundabmspace, int x, int y) {

		// Load strain cue value layers
		final GridValueLayer MM = (GridValueLayer) woundabmspace.getValueLayer("MM");
		final GridValueLayer MX = (GridValueLayer) woundabmspace.getValueLayer("MX");
		final GridValueLayer MY = (GridValueLayer) woundabmspace.getValueLayer("MY");

		// Get strain cue at cell location
		final double MMabm = MM.get(x, y);
		final double MXabm = MX.get(x, y);
		final double MYabm = MY.get(x, y);

		return new double[] {MMabm, MXabm, MYabm};
	}

	// Returns a vector of structural cues from collagen and fibrin content
	private double[] getColFiberCue(woundABMContextSim woundabmspace, int x, int y) {

		// Get grids covered by the cell
		final List<GridPoint> coveredSites = cellCoverage(x, y);

		// Get grid value layers
		final ArrayList<GridValueLayer> collagenLayers= (ArrayList<GridValueLayer>)woundabmspace.GridValueLayerList();
		final ArrayList<GridValueLayer> fibrinLayers = (ArrayList<GridValueLayer>)woundabmspace.fibrinGridValueLayerList();

		// Calculate structural cue magnitude, x, and y components
		double fiberDist = 0;
		double fiberSumX = 0;
		double fiberSumY = 0;
		double fiberSum = 0;
		double fiberAngle;
		for (GridPoint site : coveredSites) {	// Iterate through covered grids
			final int xCoor = site.getX();
			final int yCoor = site.getY();
			for (int i = 0; i < binNum180; i++) {	// Iterate through fiber bins
				fiberAngle = Math.toRadians(2*(binSize*i-87.5));
				fiberDist = collagenLayers.get(i).get(xCoor, yCoor)+fibrinLayers.get(i).get(xCoor, yCoor);
				fiberSumX = fiberSumX+fiberDist*Math.cos(fiberAngle);
				fiberSumY = fiberSumY+fiberDist*Math.sin(fiberAngle);
				fiberSum = fiberSum+fiberDist;
			}
		}
		final double meanFiberSumX = fiberSumX/fiberSum;
		final double meanFiberSumY = fiberSumY/fiberSum;
		final double meanFiberAngle = 0.5*Math.atan2(meanFiberSumY, meanFiberSumX);
		
		// Set structural cue magnitude, x, and y components
		final double SM = Math.sqrt(Math.pow(meanFiberSumX, 2)+Math.pow(meanFiberSumY, 2));
		final double SX = Math.cos(meanFiberAngle);
		final double SY = Math.sin(meanFiberAngle);

		return new double[] {SM, SX, SY, meanFiberAngle};
	}

	// Returns the direction and magnitude of cell persistence
	private double[] getPersistence(woundABMContextSim woundabmspace, int x, int y) {

		// Calculate cell direction magnitude um/hr
		final GridValueLayer chemokine = (GridValueLayer) woundabmspace.getValueLayer("Chemokine");
		final double PM = adjCellSpeed*((chemokine.get(x, y)-concMin)/(concMax-concMin)*
				(cellSpeedMax-cellSpeedMin))+cellSpeedMin;
		final double PX = Math.cos(Math.toRadians(this.angleSelection));
		final double PY = Math.sin(Math.toRadians(this.angleSelection));

		return new double[] {PM, PX, PY};
	}

	// Returns a set of angles
	private Double[] getAngleSet(double theta, double rho) {

		// Get probability distribution
		double maxBinProb = 0.0;
		final double[] binProb = WND(theta, rho);
		
		// Determine maximum bin probability
		final int wndBinNum = binProb.length;
		for (int i = 0; i < wndBinNum; i++) {
			if (binProb[i] >= maxBinProb)
				maxBinProb = binProb[i];
		}
		
		// Sum probability
		int[] binProb100 = new int[binNum360];
		double[] binDir360 = new double[binNum360];
		int totalWeight = 0;
		for (int i = 0; i < binNum360; i++) {
			binDir360[i] = binSize*i-177.5;
			binProb100[i] = (int) Math.round(100*(binProb[i]/maxBinProb));
			totalWeight = totalWeight+binProb100[i];
		}
		
		// Create weighted angle set
		Double[] angleSet = new Double[totalWeight];
		int k = 0;
		for (int i = 0; i < binNum360; i++) {
			if (binProb100[i] > 0) {
				final int binProbNum = binProb100[i];
				for (int j = 0; j < binProbNum; j++) {
					angleSet[k] = (double) binDir360[i];
					k++;
				}
			}
		}
		Collections.shuffle(Arrays.asList(angleSet));
		return angleSet;
	}

	// Returns a wrapped normal distribution
	private double[] WND(double theta, double rho) {

		// Calculate WND bin (-180 to 180)
		final double sigma = -2*Math.log(rho);
		final double sigma2 = 2*sigma;
		double[] wndArray = new double[361];
		for (int i = 0; i < 361; i++) {
			final double dir = Math.toRadians(i-180);
			final double dtDiff = dir-theta;
			double wnd = Math.exp(-Math.pow((dtDiff),2)/(sigma2));
			for (int p = 1; p <= 10; p++) {
				final double PIp2 = 2*Math.PI*p;
				wnd = wnd+Math.exp(-Math.pow((dtDiff+PIp2),2)/(sigma2))+
						Math.exp(-Math.pow((dtDiff-PIp2),2)/(sigma2));
			}
			wndArray[i] = (double) wnd/Math.sqrt(2*Math.PI*sigma);
		}
		
		// Calculate cumulative probability of WND
		double[] cumProb = new double[361];
		cumProb[0] = Math.toRadians(-180)*wndArray[0];
		for (int i = 0; i < 360; i++) {
			cumProb[i+1] = cumProb[i]+0.5*Math.toRadians(1)*(wndArray[i+1]+wndArray[i]);
		}
		
		// Calculate probability function of WND
		double[] binProb = new double[72];
		for (int i = 0; i < binNum360; i++) {
			binProb[i] = cumProb[binSize*(i+1)]-cumProb[binSize*i];
		}
		
		return binProb;
	}

	// Returns a list of sites which are available for a potential daughter cell
	private List<GridPoint> getMigrationSite(Grid<?> grid, int x, int y, int extent, double angle) {

		final int[] angleWindow = new int[] {45,30,18,15,13,11,10,9,8,7};

		double smallestAngle = 360;
		int siteX = -1;
		int siteY = -1;

		// Iterate through potential sites
		List<GridPoint> emptySites = new ArrayList<GridPoint>();
		for (int i = -extent; i <= extent; i++) {
			int xCoor = x+i;
			for (int j = -extent; j <= extent; j++) {
				int yCoor = y+j;

				// Check if the site is within the radial distance
				final double dist = Math.sqrt(Math.pow(xCoor-x,2)+Math.pow(yCoor-y,2));
				if ((dist <= extent && dist > extent-1) || extent == 1) {

					// Correct out of range site angles
					double siteAngle = Math.toDegrees(Math.atan2(yCoor-y, xCoor-x));
					if (siteAngle < 0) {
						siteAngle = siteAngle+360.0;
					}
					if (angle > 270 && siteAngle == 0) {
						siteAngle = 360;
					}

					// Check if the site is the closest site within the angle window
					final double siteAngleDiff = Math.abs(siteAngle-angle);
					if (siteAngleDiff <= angleWindow[extent-1] && siteAngleDiff <= smallestAngle) {

						// Wrapped space
						if (xCoor < 0) {
							xCoor = xCoor+gridWidth;
						} else if (xCoor >= gridWidth) {
							xCoor = xCoor-gridWidth;
						}
						if (yCoor < 0) {
							yCoor = yCoor+gridHeight;
						} else if (yCoor >= gridHeight) {
							yCoor = yCoor-gridHeight;
						}

						// Record site
						smallestAngle = siteAngleDiff;
						siteX = xCoor;
						siteY = yCoor;
					}
				}
			}
		}

		// Check if and appropriate site has been identified
		if (siteX == -1 || siteY == -1) {

		}
		// Check if the appropriate site's surrounding area is free
		else if (!grid.getObjectsAt(siteX, siteY).iterator().hasNext()) {

			if (extent < gridDiameter) {
				if (findNumNeighbor(grid, siteX, siteY) <= 1) {
					emptySites.add(new GridPoint(siteX, siteY));
				}
			} else {
				if (findNumNeighbor(grid, siteX, siteY) < 1) {
					emptySites.add(new GridPoint(siteX, siteY));	
				}
			}
		}

		return emptySites;
	}

	// Returns a list of sites which are available for a potential daughter cell
	private List<GridPoint> getEmptySites(Grid<?> grid, int x, int y, int extent) {

		// Iterate through potential sites
		List<GridPoint> emptySites = new ArrayList<GridPoint>();
		for (int i = -extent; i <= extent; i++) {
			int xCoor = x+i;
			for (int j = -extent; j <= extent; j++) {
				int yCoor = y+j;

				// Check if the site is within the radial distance
				final double dist = Math.sqrt(Math.pow(xCoor-x,2)+Math.pow(yCoor-y,2));
				if (dist <= extent+0.25*gridUnitSize || extent == 1) {

					// Wrapped space
					if (xCoor < 0) {
						xCoor = xCoor+gridWidth;
					} else if (xCoor >= gridWidth) {
						xCoor = xCoor-gridWidth;
					}
					if (yCoor < 0) {
						yCoor = yCoor+gridHeight;
					} else if (yCoor >= gridHeight) {
						yCoor = yCoor-gridHeight;
					}

					// Check if site is free
					if (!grid.getObjectsAt(xCoor, yCoor).iterator().hasNext()) {

						// Check if area around site is free
						if (extent < gridDiameter) {
							if (findNumNeighbor(grid, xCoor, yCoor) <= 1) {
								emptySites.add(new GridPoint(xCoor,yCoor));
							}
						}
						else {
							if (findNumNeighbor(grid, xCoor, yCoor) < 1) {
								emptySites.add(new GridPoint(xCoor,yCoor));	
							}
						}
					}
				}
			}
		}

		Collections.shuffle(emptySites);	
		return emptySites;
	}

	// Returns the number of agents in the exclusion area of a grid point
	public int findNumNeighbor(Grid<?> grid, int x, int y) {

		// Iterate through neighboring sites
		int numNeighbors = 0;
		for (int i = -gridDiameter; i <= gridDiameter; i++) {
			int xCoor = x+i;
			for (int j = -gridDiameter; j <= gridDiameter; j++) {
				int yCoor = y+j;

				// Check if the site is within the radial distance
				final double dist = Math.sqrt(Math.pow(xCoor-x,2)+Math.pow(yCoor-y,2));
				if (dist < gridDiameter) {

					// Wrapped space
					if (xCoor < 0) {
						xCoor = xCoor+gridWidth;
					} else if (xCoor >= gridWidth) {
						xCoor = xCoor-gridWidth;
					}
					if (yCoor < 0) {
						yCoor = yCoor+gridHeight;
					} else if (yCoor >= gridHeight) {
						yCoor = yCoor-gridHeight;
					}

					// Check if site is occupied
					if (grid.getObjectsAt(xCoor, yCoor).iterator().hasNext()) {
						numNeighbors = numNeighbors+1;
					}
				}
			}
		}

		return numNeighbors;
	}
}