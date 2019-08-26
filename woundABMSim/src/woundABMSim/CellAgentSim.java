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

	// Pull parameters
	Parameters p = RunEnvironment.getInstance().getParameters();

	// Geometry
	double sampleWidth = (Double) p.getValue("sampleWidth");	// um
	double sampleHeight = (Double) p.getValue("sampleHeight");	// um
	double gridUnitSize = (Double) p.getValue("gridUnitSize");	// um
	int gridWidth = (int) Math.floor(sampleWidth/gridUnitSize);	// grids
	int gridHeight = (int) Math.floor(sampleHeight/gridUnitSize);	// grids
	private double cellRadius = 5.0;	// um
	private double cellArea = Math.pow(cellRadius, 2)*Math.PI;	// um^2
	private int gridDiameter = (int) Math.round(cellRadius*2/gridUnitSize);	// Cell diameter in grids
	private double woundRadius = 113.0;	// um
	private int cellGridArea;	// Cell area in grids

	// Tracking
	private double apoptosisTime;	// hr
	private double apoptosisAge;	// hr
	private double mitosisTime;	// hr
	private double mitosisAge;	// hr
	private double migrationTime;	// hr
	private double moveCounter;	// hr
	private double depositionTime;	// hr
	private double depositionCounter;	// hr
	private double degradationTime;	// hr
	private double degradationCounter;	// hr

	// Weighting factors
	private double Wp = (Double) p.getValue("Wp");	// Persistance cue weight
	private double Ws = (Double) p.getValue("Ws"); // Structural cue weight
	private double Wm = (Double) p.getValue("Wm");	// Mechanics cue weight
	private double Wc = (Double) p.getValue("Wc");	// Chemokine cue weight

	// Correction factors
	private double alphaGeometric = woundRadius/2828.0;	// Cryo infarct radius = 2828 um;
	private double alphaInfiltration = woundRadius/500.0; // Cryo infarct thickness = 500 um;

	// Migration
	private double cellSpeedMin = alphaInfiltration*1.0; // um/hr	(doesn't need to be defined this way here)
	private double cellSpeedMax = alphaInfiltration*10.0;	// um/hr
	public double angleSelection;	// deg
	private boolean moveboolean;

	// Chemokine
	private double concMin = 0.0;	//
	private double concMax = 1.0;	//

	// Boolean
	private boolean cellRemoval = (Boolean) p.getValue("cellRemoval");
	private boolean colRotation = (Boolean) p.getValue("colRotation");

	// Collagen
	private int fiberBinSIZE = 5;	// 5 degrees per collagen bin
	private double kColRotMax = 3.0;	// deg/hr
	private double kColRotMin = 0.3;	// deg/hr
	private double adjkColRot = 0.5;	// 
	private double kColGenMax = 17.9874;	// fiber/cell/hr
	private double kColGenMin = 0.1799;	// fiber/cell/hr
	private double kColGen;	// fiber/cell/hr
	private double adjColGen;    	// adjust collagen generation to cell grid area
	private double alignedDepFrac = 0.8;      	// 
	private double kColDegMax = 0.0038;   // unitless scaling factor
	private double kColDegMin = 3.7676e-04;	// unitless scaling factor
	private double kColDeg;	// %
	private double kFibDeg;	// %
	private String depoType = (String) p.getValue("depoType"); // "Aligned" or "Random"

	// Cell cycle timing
	private double timeStep =0.5;	// hr
	private double gMitosisTime = (Double) p.getValue("mitosisTime");
	private double gApoptosisTime = (Double) p.getValue("apoptosisTime");
	private double gDepositionTime = (Double) p.getValue("depositionTime");
	private double gDegradationTime = (Double) p.getValue("degradationTime");

	/* Constructors */
	public CellAgentSim() {
		// Assign parameters
		this.mitosisTime = gMitosisTime;
		this.mitosisAge = RandomHelper.nextDoubleFromTo(0, gMitosisTime);
		this.apoptosisTime = gApoptosisTime;
		this.apoptosisAge = mitosisAge;
		this.depositionTime = gDepositionTime;
		this.depositionCounter = RandomHelper.nextDoubleFromTo(0, gDepositionTime);
		this.degradationTime = gDegradationTime;
		this.degradationCounter = RandomHelper.nextDoubleFromTo(0, gDegradationTime);
		this.angleSelection = RandomHelper.nextDoubleFromTo(-180, 180);
		this.moveboolean = true;

		List<GridPoint> coveredSites = cellCoverage(0, 0);
		cellGridArea = coveredSites.size();
		adjColGen = Math.pow(gridUnitSize, 2)/cellArea;	// Math.pow(gridUnitSize, 2)/cellArea //1.0/cellGridArea
	}

	public CellAgentSim(double apoptosisTime, double mitosisTime, double depositionTime,
			double depositionCounter, double degradationTime, double degradationCounter, 
			double apoptosisAge, double mitosisAge, double angleSelection, boolean moveboolean) {

		// Assign parameters
		this.apoptosisTime = apoptosisTime;
		this.apoptosisAge = apoptosisAge;
		this.mitosisTime = mitosisTime;
		this.mitosisAge = mitosisAge;
		this.depositionTime = depositionTime;
		this.depositionCounter = depositionCounter;
		this.degradationTime = degradationTime;
		this.degradationCounter = degradationCounter;
		this.angleSelection = angleSelection;
		this.moveboolean = moveboolean;

		List<GridPoint> coveredSites = cellCoverage(0, 0);
		cellGridArea = coveredSites.size();
		adjColGen = Math.pow(gridUnitSize, 2)/cellArea;	// Math.pow(gridUnitSize, 2)/cellArea //1.0/cellGridArea
	}

	/* Primary scheduled methods in order of execution */
	// Infarction method
	@ScheduledMethod(start = 0, priority = 3)
	public void cellRemoval() {

		// Remove cells in the wound area
		if (cellRemoval == true) {
			woundABMContextSim woundabmspace = (woundABMContextSim) ContextUtils.getContext(this);
			Grid<?> grid = (Grid<?>) woundabmspace.getProjection("Cell Grid");
			GridPoint pt = grid.getLocation(this);	

			double removalDist = (woundRadius/gridUnitSize);
			double dist = Math.sqrt(Math.pow((pt.getX()-gridWidth/2), 2)+Math.pow((pt.getY()-gridHeight/2), 2));
			if (dist <= (double) removalDist) {
				woundabmspace.remove(this);
			}
		}
	}

	// Initialize chemokine-dependent apoptosis, mitosis, and migration time (add to step method) 
	@ScheduledMethod(start = 0.2, interval = 0.5, priority = 2)
	public void setChemokineProcess() {
		woundABMContextSim woundabmspace = (woundABMContextSim) ContextUtils.getContext(this);
		Grid<?> grid = (Grid<?>) woundabmspace.getProjection("Cell Grid");
		GridValueLayer chemokine = (GridValueLayer) woundabmspace.getValueLayer("Chemokine");
		GridPoint pt = grid.getLocation(this);

		double effectiveChemokine = (chemokine.get(pt.getX(), pt.getY())-concMin)/(concMax-concMin);

		// Define chemokine dependent mitosis to ensure cells are always coming from the outside of wound
		double mitosisMin = gMitosisTime;
		double mitosisMax = gApoptosisTime;
		this.mitosisTime = effectiveChemokine*(mitosisMin-mitosisMax)+mitosisMax;

		// Define chemokine dependent migration to speed up wound closure
		double cellSpeed = effectiveChemokine*(cellSpeedMax-cellSpeedMin)+cellSpeedMin;
		this.migrationTime = gridUnitSize/cellSpeed;
	}

	// Agent decision tree
	@ScheduledMethod(start = 0.5, interval = 0.5, priority = 1)
	public void step() {
		woundABMContextSim woundABMSpace = (woundABMContextSim) ContextUtils.getContext(this);
		Grid<?> grid = (Grid<?>) woundABMSpace.getProjection("Cell Grid");
		GridPoint pt = grid.getLocation(this);

		// Get location
		int x = pt.getX();
		int y = pt.getY();

		// Check for mitosis
		if (mitosisAge >= mitosisTime) {	
			mitosisAge = randomClock();	// Quiescent cell sets mitosis age to random number
			mitose(woundABMSpace, grid, x, y);

		} else {

			// Update cell orientation
			this.angleSelection = guidanceCue(woundABMSpace, x, y);

			// Check for migration
			if (moveCounter >= migrationTime) {	
				moveTowards(grid, x, y);
				moveCounter = 0.0;
			}

			// Get grids covered by the cell
			List<GridPoint> coveredSites = cellCoverage(x, y);

			// Deposit, degrade, and rotate collagen in the grids beneath the cell
			for (GridPoint site : coveredSites) {

				// Get site x and y coordinates
				int siteX = site.getX();
				int siteY = site.getY();

				// Retrieve value layers
				GridValueLayer chemokine = (GridValueLayer) woundABMSpace.getValueLayer("Chemokine");
				ArrayList<GridValueLayer> collagenLayers = (ArrayList<GridValueLayer>)woundABMSpace.GridValueLayerList();

				// Check for rotation
				if (colRotation == true) {
					collagenRotation(siteX, siteY, chemokine, collagenLayers);
				}

				// Check for degradation
				if (degradationCounter >= degradationTime) {
					ArrayList<GridValueLayer> fibrinLayers = (ArrayList<GridValueLayer>)woundABMSpace.fibrinGridValueLayerList();	
					matrixDegradation(siteX, siteY, chemokine, collagenLayers, fibrinLayers, degradationTime);
					degradationCounter = timeStep;
				} else {
					degradationCounter += timeStep;
				}

				// Check for deposition
				if (depositionCounter >= depositionTime) {	
					collagenDeposition(siteX, siteY, chemokine, collagenLayers, depositionTime);
					depositionCounter = timeStep;
				} else {
					depositionCounter += timeStep;
				}
			}
		}

		// Check for apoptosis
		if (apoptosisAge >= apoptosisTime) {
			woundABMSpace.remove(this);
		}

		// Increment the counters
		apoptosisAge = apoptosisAge+timeStep;
		mitosisAge = mitosisAge+timeStep;
		moveCounter = moveCounter+timeStep;	
	}

	/* Secondary methods */
	// Mitosis method (fix get(site)) (move constants to class loading or static method)
	@SuppressWarnings({ "rawtypes", "unchecked" })
	public void mitose(woundABMContextSim woundabmspace, Grid grid, int x, int y) {

		// Create list of vacant neighboring sites
		List<GridPoint> emptySites = getEmptySites(grid, x, y, gridDiameter);
		if (emptySites.size() > 0) {

			// Select a random site from list of empty sites
			GridPoint site = emptySites.get(RandomHelper.nextIntFromTo(0,(emptySites.size()-1)));
			int xCoor = site.getX();
			int yCoor = site.getY();

			// Add daughter (angle from mother's location to new site, random ages)
			double dAngle = Math.toDegrees(Math.atan2(yCoor-y, xCoor-x));
			double dMitosisClock = randomClock();
			CellAgentSim cellagent = new CellAgentSim(apoptosisTime, mitosisTime, depositionTime,
					0.0, degradationTime, 0.0, dMitosisClock, dMitosisClock, dAngle, true);
			woundabmspace.add(cellagent);
			grid.moveTo(cellagent, xCoor, yCoor);
		}
	}

	// Returns an angle based on guidance cues (set colMVL and MVA?) (correct rho?) (move cellSpeed Allocation)
	private double guidanceCue(woundABMContextSim woundabmspace, int x, int y) {

		// Maximum chemokine gradient based on diffusion, generation, and degradation parameters
		double kcdeg = 0.001;         // 1/s
		double Dc = 100;              	// um2/s
		double L = 2828;              	// um
		double lam = Math.sqrt(kcdeg/Dc);  	// 1/um
		double maxConcGrad = (lam*L*Math.sinh(lam*L))/(Math.cosh(lam*L)+Math.sinh(lam*L))/(alphaGeometric*L);    // 1/um

		// Maximum strain anisotropy based on cryoinfarct data
		double eii = 0.05;
		double ejj = 0;
		double eij = 0;

		// Maximum cell speed
		cellSpeedMax = alphaInfiltration*10; // um/tick

		// Calculate guidance cue normalization factors based on maximum values
		double Mc = 0.5*cellRadius*maxConcGrad;
		double Mm = Math.sqrt(Math.pow(eii, 2)-2*eii*ejj+Math.pow(ejj, 2)+4*Math.pow(eij,2))/4; // Ms = 0.0125	 
		double Ms = 1;
		double Mp = cellSpeedMax;

		// Acquire chemokine, strain, fiber, and persistence cues
		double[] chemoCue = getChemoCue(woundabmspace, x, y);
		double CM = chemoCue[0]; 
		double CX = chemoCue[1];
		double CY = chemoCue[2];
		double[] strainCue = getStrainCue(woundabmspace, x, y); 
		double MM = strainCue[0]; 
		double MX = strainCue[1];
		double MY = strainCue[2];
		double[] fiberCue = getColFiberCue(woundabmspace, x, y);
		double SM = fiberCue[0]; 
		double SX = fiberCue[1];
		double SY = fiberCue[2];	
		double[] persistence = getPersistence(woundabmspace, x, y);
		double PM = persistence[0];
		double PX = persistence[1]; 
		double PY = persistence[2];

		// Write GridValueLayer of colMVA and colMVL
		GridValueLayer colMVA = (GridValueLayer) woundabmspace.getValueLayer("Collagen MVA");
		GridValueLayer colMVL = (GridValueLayer) woundabmspace.getValueLayer("Collagen MVL");
		colMVL.set(fiberCue[0], x, y);
		colMVA.set(Math.toDegrees(fiberCue[3]), x, y);

		// Calculate guidance cues strength using maximum value normalization an weighting factors
		double Sc = (Wc*CM)/Mc;
		double Sm = (Wm*MM)/Mm;
		double Ss = (Ws*SM)/Ms;
		double Sp = (Wp*PM)/Mp;

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

		// Find maximum resultant magnitude
		double maxRM = 0;
		ArrayList<Integer> RMlist = new ArrayList<Integer>();
		int maxTempRMidx = 99;
		for (int i = 0; i < 4; i++) {
			if (RM[i] > maxRM) {
				maxRM = RM[i];
				maxTempRMidx = i;
			} else if (RM[i] == maxRM) {
				if (RMlist.size() < 1) {
					RMlist.add(maxTempRMidx);
				}
				RMlist.add(i);
			}
		}

		// If RM are equal pick randomly out of two
		if (RMlist.size() < 1) {
			RMlist.add(maxTempRMidx);
		}
		int maxRMindex = RMlist.get(RandomHelper.nextIntFromTo(0,(RMlist.size()-1)));

		// Persistence tuning factor
		double pMin = 0.25;     // min persistence time (hr)
		double pNoCue = 1.2;    // persistence time in the absence of external directional cues
		double rhoTuning = 2*Wp*pMin/(pNoCue-pMin); // rhoTuning = 0.1754 or = 0.35;

		// Calculate the resultant angle "theta" and normalized resultant magnitude "rho"
		double theta = Math.atan2(RY[maxRMindex], RX[maxRMindex]);
		double rho = RM[maxRMindex]/(rhoTuning+Ss+Sm+Sc+Sp);

		// Correct rho outside of allowed range
		if (rho > 0.9999) {
			rho = 0.9999;
		}
		if (rho < 0.0001) {
			rho = 0.0001;
		}

		// Select a direction from the probability distribution, exceptions for rho >=0.999 and rho <=0.001
		int binSize = 5;
		int binNum = 360/binSize;
		double[] binDir360 = new double[72];
		for (int i = 0; i < binNum; i++) {
			binDir360[i] = binSize*i-177.5;
		}
		double guidedAngle; 
		if (rho >= 0.999) {
			guidedAngle = Math.ceil(Math.toDegrees(theta));
		} else if (rho <= 0.001) {
			guidedAngle = binDir360[(int) RandomHelper.nextIntFromTo(0,(binNum-1))];
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
	public void moveTowards(Grid grid, int x, int y) {

		// Search for a migration site starting from farthest to nearest
		int extent = gridDiameter;
		while (extent > 0) {

			// Get list of sites within extent and the angle window that are available for migration
			double angle = this.angleSelection;
			if (angle < 0) {
				angle = angle+360.0;
			}
			List<GridPoint> emptySites = getMigrationSite(grid, x, y, extent, angle);
			if (!emptySites.isEmpty()) {

				// Select a random migration site from the list of those available to move to
				int migrationSiteIdx = RandomHelper.nextIntFromTo(0,(emptySites.size()-1));
				GridPoint migrationSite = emptySites.get(migrationSiteIdx);
				int desiredX = migrationSite.getX();
				int desiredY = migrationSite.getY();

				// Move to new location
				grid.moveTo(this, desiredX, desiredY);
				this.moveboolean = true;
				extent = 0;

				// Cancel migration
			} else if (extent == 1) {
				this.moveboolean = false;
			}

			// Reduce search distance
			extent = extent-1;
		}
	}	

	// Returns a list containing all the grids who's center is within 2.5um of the cell boarder
	private List<GridPoint> cellCoverage(int x, int y) {

		// Iterate through neighbors
		int extent = (int) Math.floor(cellRadius/gridUnitSize)+1;
		List<GridPoint> coveredSites = new ArrayList<GridPoint>();
		for (int i = -extent; i <= extent; i++) {
			for (int j = -extent; j <= extent; j++) {
				int xCoor = x+i;
				int yCoor = y+j;

				// Check radial distance
				double dist = Math.sqrt(Math.pow(xCoor-x,2)+Math.pow(yCoor-y,2))*gridUnitSize;
				if (dist <= cellRadius+1) {

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

	// Rotate collagen (shortened bin assignment)
	public void collagenRotation(int x, int y, GridValueLayer chemokine, ArrayList<GridValueLayer> gridValueLayerList) {

		// Calculate collagen fiber rotation constant based on chemokine concentration
		double kColRot = ((chemokine.get(x,y)-concMin)/(concMax-concMin)*(kColRotMax-kColRotMin))+kColRotMin;
		double Q = Math.ceil(kColRot*timeStep/2);

		// Get all 36 bins of collagen
		int binNum= gridValueLayerList.size();
		double [] colContent = new double[binNum];
		ArrayList <Double> tempColContent = new ArrayList<Double>();
		for (int j = 0; j < binNum; j++) {
			colContent[j] = gridValueLayerList.get(j).get(x, y);
			tempColContent.add((double) 0.0);
		}

		// Converts angleSelection range of [-180 180] to range of [-90 90]
		double radAngle = Math.toRadians(this.angleSelection);
		double cellAngle = Math.toDegrees(Math.atan(Math.sin(radAngle)/Math.cos(radAngle)));
		double fiberBin;
		for (int j = 0; j < binNum; j++) {	

			// Get fiber offset from cell angle in radians
			fiberBin = fiberBinSIZE*j-87.5;
			double newAngle = 360;
			double fiberOffset = Math.toRadians(cellAngle-fiberBin);
			double tanFiberOff = Math.tan(fiberOffset);
			double sinFiberOff = Math.abs(Math.sin(fiberOffset));
			if ((tanFiberOff > 0) && (sinFiberOff != 1)) {// fiberOffset < 90
				// (M) newAngle = fiberBin+(kColRot*adjkColRot*timeStep*(sinFiberOff));
				newAngle = fiberBin+(kColRot*timeStep/Q*(sinFiberOff));
			} else if ((tanFiberOff <= 0) && (sinFiberOff != 1)) {// fiberOffset > 90
				// (M) newAngle = fiberBin-(kColRot*adjkColRot*timeStep*(sinFiberOff));
				newAngle = fiberBin-(kColRot*timeStep/Q*(sinFiberOff));
			} else if (Math.abs(Math.sin(fiberOffset)) == 1) {// fiberOffset ~= 90
				double flip = 1;
				int dice = RandomHelper.nextIntFromTo(0,1);
				if (dice == 1) {
					flip = -1;
				}
				// (M) newAngle = fiberBin+(flip*kColRot*adjkColRot*timeStep*(sinFiberOff));
				newAngle = fiberBin+(flip*kColRot*timeStep/Q*(sinFiberOff));
			}

			// Detect and correct angles outside of fiber bin range
			if (newAngle <= -90) {
				newAngle = newAngle+180;
			} else if (newAngle > 90) {
				newAngle = newAngle-180;
			}

			// Map newAngle to the nearest bin (can be partially substituted by code in Deposition method)
			double min = 360;
			int minIndex = 99;
			double rotBins = 360;
			for (int i = 0; i < binNum; i++) {
				double fiberDirBins = fiberBinSIZE*i-87.5;
				double angleDiff = Math.abs(fiberDirBins-newAngle);
				if (angleDiff < min) {
					min = Math.abs(fiberDirBins-newAngle);
					minIndex = i;
					rotBins = fiberDirBins;
				} else if (angleDiff == min) {
					int dice = RandomHelper.nextIntFromTo(0,1);
					if (dice == 1) {
						min = Math.abs(fiberDirBins-newAngle);
						minIndex = i;
						rotBins = fiberDirBins;
					}
				}
			}

			// Rotate collagen
			int count2indx = 99;
			double rotFrac = (newAngle-rotBins)/5;
			double temp1;
			double temp2;
			if (rotFrac >= 0) { //kept
				temp1 = tempColContent.get(minIndex)+colContent[j]*(1-rotFrac);
				if (minIndex == 35) {// rotated
					count2indx = 0;
				} else {
					count2indx = minIndex+1;
				}
				temp2 = tempColContent.get(count2indx)+colContent[j]*rotFrac;	
			} else { // kept
				temp1 = tempColContent.get(minIndex)+colContent[j]*(1+rotFrac);
				if (minIndex == 0) {
					count2indx = 35;
				} else {
					count2indx = minIndex-1;
				}
				temp2 = tempColContent.get(count2indx)-colContent[j]*rotFrac;
			}

			tempColContent.set(minIndex, temp1);
			tempColContent.set(count2indx, temp2);	
		}

		// Update all 36 collagen bins
		for (int i = 0; i < binNum; i++) {
			gridValueLayerList.get(i).set(tempColContent.get(i), x, y);
		}
	}

	// Fiber degradation method
	public void matrixDegradation(int x, int y, GridValueLayer chemokine, ArrayList<GridValueLayer> gridValueLayerList, ArrayList<GridValueLayer> fibringridValueLayerList, double degrationInterval) {

		// Calculate fiber degradation constants with respect to chemokine concentration
		kColDeg = ((chemokine.get(x, y)-concMin)/(concMax-concMin)*(kColDegMax-kColDegMin)+kColDegMin);
		kFibDeg = kColDeg;

		// Degrade fibers uniformly across all bins
		int binNum = gridValueLayerList.size();
		for (int i = 0; i < binNum; i++) {

			// Collagen degradation 
			double colContent = gridValueLayerList.get(i).get(x, y);
			double newColContent = colContent*(1-kColDeg*degrationInterval);
			gridValueLayerList.get(i).set(newColContent, x, y);

			// Fibrin degradation (only within wound)
			double distToCenter = Math.sqrt(Math.pow((x-gridWidth/2), 2)+Math.pow((y-gridHeight/2), 2));
			double woundSize = (double) woundRadius/gridUnitSize;
			if (distToCenter < woundSize) {
				double fibrinContent = fibringridValueLayerList.get(i).get(x, y);
				double newFibrinContent = fibrinContent*(1-kFibDeg*degrationInterval);
				fibringridValueLayerList.get(i).set(newFibrinContent, x, y);
			}
		}
	}

	// Aligned collagen deposition method (MATLAB deltaCol is divided by number of local collagen fibers) (alignedDepFrac is not in MATLAB model)
	public void collagenDeposition(int x, int y, GridValueLayer chemokine, ArrayList<GridValueLayer> gridValueLayerList, double depositionInterval) {

		// Calculate collagen fiber generation rate based on chemokine concentration
		kColGen = ((chemokine.get(x, y))-concMin)/(concMax-concMin)*(kColGenMax-kColGenMin)+kColGenMin;
		//this.chemolevel = chemokine.get(x, y);

		// Set collagen deposition angle based on deposition type
		double depoAngle;
		double cellAngle = this.angleSelection;
		if (depoType.equals("Aligned")) {
			depoAngle = cellAngle;
		} else {
			depoAngle = RandomHelper.nextDoubleFromTo(-180, 180);
		}

		// Correct angles outside of the collagen bin range
		if (cellAngle <= -90.0) {
			depoAngle = cellAngle+180.0;
		} else if (cellAngle > 90.0) {
			depoAngle = cellAngle-180.0;
		}
		depoAngle = Math.toDegrees(Math.atan(Math.sin(Math.toRadians(depoAngle))/
				Math.cos(Math.toRadians(depoAngle))));

		// Assign to the closest bin
		int bin = (int) Math.ceil(depoAngle/fiberBinSIZE)+17;

		// Deposit the aligned collagen
		double deltaCol = depositionInterval*kColGen*adjColGen;
		double newColContent = gridValueLayerList.get(bin).get(x, y)+(deltaCol*alignedDepFrac);
		gridValueLayerList.get(bin).set(newColContent, x, y);

		// Deposit a portion of the collagen randomly
		int randBin = RandomHelper.nextIntFromTo(0,35);
		double newColContent2 = gridValueLayerList.get(randBin).get(x, y)+deltaCol*(1-alignedDepFrac);
		gridValueLayerList.get(randBin).set(newColContent2, x, y);

	}

	// Returns a random clock time
	public double randomClock() {
		Random r = new Random();
		double time = Math.round(0.1*mitosisTime*r.nextGaussian()/timeStep)*timeStep;
		return time;
	}

	/* Tertiary methods */
	// Returns a vector of chemokine cues
	private double[] getChemoCue(woundABMContextSim woundabmspace, int x, int y) {

		// Load strain cue value layers
		GridValueLayer CM = (GridValueLayer) woundabmspace.getValueLayer("CM");
		GridValueLayer CX = (GridValueLayer) woundabmspace.getValueLayer("CX");
		GridValueLayer CY = (GridValueLayer) woundabmspace.getValueLayer("CY");

		// Get strain cue at cell location
		double CMabm = CM.get(x, y);
		double CXabm = CX.get(x, y); 
		double CYabm = CY.get(x, y);

		return new double[] {CMabm, CXabm, CYabm};
	}

	// Returns a vector of strain cues
	private double[] getStrainCue(woundABMContextSim woundabmspace, int x, int y) {

		// Load strain cue value layers
		GridValueLayer MM = (GridValueLayer) woundabmspace.getValueLayer("MM");
		GridValueLayer MX = (GridValueLayer) woundabmspace.getValueLayer("MX");
		GridValueLayer MY = (GridValueLayer) woundabmspace.getValueLayer("MY");

		// Get strain cue at cell location
		double MMabm = MM.get(x, y);
		double MXabm = MX.get(x, y);
		double MYabm = MY.get(x, y);

		return new double[] {MMabm, MXabm, MYabm};
	}

	// Returns a vector of structural cues from collagen and fibrin content
	private double[] getColFiberCue(woundABMContextSim woundabmspace, int x, int y) {

		// Get grids covered by the cell
		List<GridPoint> coveredSites = cellCoverage(x, y);

		// Get grid value layers
		ArrayList<GridValueLayer> collagenLayers= (ArrayList<GridValueLayer>)woundabmspace.GridValueLayerList();
		ArrayList<GridValueLayer> fibrinLayers = (ArrayList<GridValueLayer>)woundabmspace.fibrinGridValueLayerList();
		int binNum = collagenLayers.size();

		// Calculate structural cue magnitude, x, and y components
		double fiberDist = 0;
		double fiberSumX = 0;
		double fiberSumY = 0;
		double fiberSum = 0;
		double fiberAngle;
		for (GridPoint site : coveredSites) {	// Iterate through covered grids
			int xCoor = site.getX();
			int yCoor = site.getY();
			for (int i = 0; i < binNum; i++) {	// Iterate through fiber bins
				fiberAngle = Math.toRadians(2*(fiberBinSIZE*i-87.5));
				fiberDist = collagenLayers.get(i).get(xCoor, yCoor)+fibrinLayers.get(i).get(xCoor, yCoor);
				fiberSumX = fiberSumX+fiberDist*Math.cos(fiberAngle);
				fiberSumY = fiberSumY+fiberDist*Math.sin(fiberAngle);
				fiberSum = fiberSum+fiberDist;
			}
		}
		double meanFiberSumX = fiberSumX/fiberSum;
		double meanFiberSumY = fiberSumY/fiberSum;
		double meanFiberAngle = 0.5*Math.atan2(meanFiberSumY, meanFiberSumX);

		// Set structural cue magnitude, x, and y components
		double SM = Math.sqrt(Math.pow(meanFiberSumX, 2)+Math.pow(meanFiberSumY, 2));
		double SX = Math.cos(meanFiberAngle);
		double SY = Math.sin(meanFiberAngle);

		return new double[] {SM, SX, SY, meanFiberAngle};
	}

	// Returns the direction and magnitude of cell persistence
	private double[] getPersistence(woundABMContextSim woundabmspace, int x, int y) {

		// Set adjacent cell speed based on whether the cell is going to migrate
		double adjCellSpeed;
		if (!moveboolean) {
			adjCellSpeed =RandomHelper.nextDoubleFromTo(0, 1);
		} else {
			adjCellSpeed = 1.0;
		}

		// Calculate cell direction magnitude um/hr
		GridValueLayer chemokine = (GridValueLayer) woundabmspace.getValueLayer("Chemokine");
		double PM = adjCellSpeed*((chemokine.get(x, y)-concMin)/(concMax-concMin)*
				(cellSpeedMax-cellSpeedMin))+cellSpeedMin;
		double PX = Math.cos(Math.toRadians(this.angleSelection));
		double PY = Math.sin(Math.toRadians(this.angleSelection));

		return new double[] {PM, PX, PY};
	}

	// Returns a set of angles (slow?)
	private Double[] getAngleSet(double theta, double rho) {

		// Get probability distribution
		String typeModal = "Single"; // 'Single' or 'Bimodal'
		double sumBinProb = 0.0;
		double maxBinProb = 0.0;
		double theta2 = theta+Math.PI;
		double[] wndBinProb1 = WND(theta, rho);
		double[] wndBinProb2 = WND(theta2, rho);
		int wndBinNum = wndBinProb1.length;
		double[] bimodalBinProb = new double[wndBinNum];
		for (int i = 0; i < wndBinNum; i++) {
			if (typeModal.equals("Bimodal")) {
				bimodalBinProb[i] = 0.5*(wndBinProb1[i]+wndBinProb2[i]);
			} else {
				bimodalBinProb[i] = wndBinProb1[i];
			}
			sumBinProb = sumBinProb+bimodalBinProb[i];
			if (bimodalBinProb[i] >= maxBinProb)
				maxBinProb = bimodalBinProb[i];
		}

		// Sum probability
		int binSize = 5;
		int binNum = 360/binSize;
		int[] binProb100 = new int[binNum];
		double[] binDir360 = new double[binNum];	
		int totalWeight = 0;
		for (int i = 0; i < binNum; i++) {
			binDir360[i] = binSize*i-177.5;
			binProb100[i] = (int) Math.round(100*(bimodalBinProb[i]/maxBinProb));
			totalWeight = totalWeight+binProb100[i];
		}

		// Create weighted angle set 
		Double[] angleset = new Double[totalWeight];
		int k = 0;
		for (int i = 0; i < binNum; i++) {
			if (binProb100[i] > 0) {
				int binProbNum = binProb100[i];
				for (int j = 0; j < binProbNum; j++) {
					angleset[k] = (double) binDir360[i];
					k++;
				}
			}
		}

		Collections.shuffle(Arrays.asList(angleset));
		return angleset;
	}

	// Returns a wrapped normal distribution
	private double[] WND(double theta, double rho) {

		// Calculate WND bin (-180 to 180)
		double sigma = -2*Math.log(rho);
		double dir;
		double wnd;
		double[] wndArray = new double[361];
		for (int i = 0; i < 361; i++) {
			dir = Math.toRadians(i-180);
			wnd = Math.exp(-Math.pow((dir-theta),2)/(2*sigma));
			for (int p = 1; p <= 10; p++) {
				wnd = wnd+Math.exp(-Math.pow((dir-theta+2*Math.PI*p),2)/(2*sigma))+
						Math.exp(-Math.pow((dir-theta-2*Math.PI*p),2)/(2*sigma));
			}
			wnd = (double) wnd/Math.sqrt(2*Math.PI*sigma);
			wndArray[i] = wnd;
		}

		// Calculate cumulative probability of WND
		double[] cumProb = new double[361];
		cumProb[0] = Math.toRadians(-180)*wndArray[0];
		for (int i = 0; i < 360; i++) {
			cumProb[i+1] = cumProb[i]+0.5*Math.toRadians(1)*(wndArray[i+1]+wndArray[i]);
		}

		// Calculate probability function of WND
		int binNum = 360/5;
		double[] binProb = new double[72];
		for (int i = 0; i < binNum; i++) {
			binProb[i] = cumProb[5*(i+1)]-cumProb[5*i];
		}

		return binProb;
	}

	// Returns a list of sites which are available for a potential daughter cell
	private List<GridPoint> getMigrationSite(Grid<?> grid, int x, int y, int extent, double angle) {

		int[] angleWindow = new int[] {45,30,18,15,13,11,10,9,8,7};

		double smallestAngle = 360;
		int siteX = -1;
		int siteY = -1;

		// Iterate through potential sites
		List<GridPoint> emptySites = new ArrayList<GridPoint>();
		for (int i = -extent; i <= extent; i++) {
			for (int j = -extent; j <= extent; j++) {
				int xCoor = x+i;
				int yCoor = y+j;

				// Check if the site is within the radial distance
				double dist = Math.sqrt(Math.pow(xCoor-x,2)+Math.pow(yCoor-y,2));
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
					double siteAngleDiff = Math.abs(siteAngle-angle);
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
			}
			else {
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
			for (int j = -extent; j <= extent; j++) {
				int xCoor = x+i;
				int yCoor = y+j;

				// Check if the site is within the radial distance
				double dist = Math.sqrt(Math.pow(xCoor-x,2)+Math.pow(yCoor-y,2));
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
	int findNumNeighbor(Grid<?> grid, int x, int y) {

		// Iterate through neighboring sites
		int extent = gridDiameter;
		int numNeighbors = 0;
		for (int i = -extent; i <= extent; i++) {
			for (int j = -extent; j <= extent; j++) {
				int xCoor = x+i;
				int yCoor = y+j;

				// Check if the site is within the radial distance
				double dist = Math.sqrt(Math.pow(xCoor-x,2)+Math.pow(yCoor-y,2));
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