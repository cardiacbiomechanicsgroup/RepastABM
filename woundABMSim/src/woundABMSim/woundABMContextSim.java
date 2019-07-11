package woundABMSim;

import au.com.bytecode.opencsv.CSVReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import org.apache.commons.math3.special.Erf;
import repast.simphony.context.DefaultContext;
import repast.simphony.context.space.grid.GridFactoryFinder;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.schedule.ISchedule;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridBuilderParameters;
import repast.simphony.space.grid.GridPoint;
import repast.simphony.space.grid.RandomGridAdder;
import repast.simphony.space.grid.WrapAroundBorders;
import repast.simphony.valueLayer.GridValueLayer;

public class woundABMContextSim extends DefaultContext<Object> {
	
	// Geometry
	static int sampleWidth = 480;											// um
	static int sampleHeight = 480;											// um
	static double gridUnitSize = 2.5;										// um (10, 5, or 2.5)
	static int gridWidth = (int) Math.floor(sampleWidth/gridUnitSize);		// grids
	static int gridHeight = (int) Math.floor(sampleHeight/gridUnitSize);	// grids
	private double cellRadius = 5.0;										// um
	private double woundRadius = 113.0; 									// um
	
	// Value layer information
	ArrayList<ArrayList<double[]>> initialFiberArrayList;
	private ArrayList<GridValueLayer> gridvaluelayerlist;
	private ArrayList<GridValueLayer> fibringridvaluelayerlist;
	private String[] basicValueLayer = {"Collagen MVA", "Collagen MVL", "Collagen Sum", "Fibrin Sum"};
	private String[] cueValueLayer = {"CGM", "CGX", "CGY", "SM", "SX", "SY"};
	private String mechs = "Biax";								// "UniaxC", "UniaxL", or "Biax"
	private double fiberDensity = 178.25;							// fibers/um^2/7um thickness
	private int fiberBinSize = 5;									// degrees
	
	private String fiberDist = (String) "Circumferential";	// "Circumferential" "Longitudinal" or "Uniform"
	private boolean includeFib = (boolean) true;			// true or false
	
	// Output file names
	private Date date = new Date();
	private DateFormat df = new SimpleDateFormat("yyyyMMdd_HHmmss");
	private String modelTag = df.format(date);
	private String woundParamFilename = "output\\"+mechs+"_"+gridUnitSize+"_"+modelTag+"_WoundParameters.csv";
	private String woundStatFilename = "output\\"+mechs+"_"+gridUnitSize+"_"+modelTag+"_WoundStat.csv";
	private String woundColAnDistFilename = "output\\"+mechs+"_"+gridUnitSize+"_"+modelTag+"_WoundColFiberAngDist.csv";
	
/* Constructor */ 
	public woundABMContextSim() {
		// Create projection (imposes structures on agents)
		super("woundABMContextSim");
		GridFactoryFinder.createGridFactory(null)
				.createGrid("Cell Grid",this,new GridBuilderParameters<Object>(new WrapAroundBorders(),
						new RandomGridAdder<Object>(), true, gridWidth, gridHeight));
		
		// Create collagen grid value layer
		GridValueLayer collagen = new GridValueLayer("Collagen", 1.0, true, gridWidth, gridHeight);
		this.addValueLayer(collagen);
		
		// Create chemokine grid value layer
		GridValueLayer chemokine = new GridValueLayer("Chemokine", true, new WrapAroundBorders(), gridWidth, gridHeight);
		this.addValueLayer(chemokine);
		
		// Create collagen and cell mean vector angle, mean vector length, and sum grid value layers
		GridValueLayer basiclayer;
		int basicLayerNum = basicValueLayer.length;
		for (int i = 0; i < basicLayerNum; i++) {
			basiclayer = new GridValueLayer(basicValueLayer[i], 1.0, true,gridWidth, gridHeight);
			this.addValueLayer(basiclayer);
		}
		
		// Create chemokine and strain cue magnitude, x component, and y component grid value layers
		GridValueLayer cuelayer;
		int cueLayerNum = cueValueLayer.length;
		for (int i = 0; i < cueLayerNum; i++) {
			cuelayer = new GridValueLayer(cueValueLayer[i], 1.0, true,gridWidth, gridHeight);
			this.addValueLayer(cuelayer);
		}
		
		// Create 36 collagen bin grid value layers
		String[] collagenLayerName = new String[36];
		GridValueLayer layer;
		int fiberLayerNum = 36;
		for (int i = 0; i < fiberLayerNum; i++) {
			collagenLayerName[i] = "CollagenLayer"+i;
			layer = new GridValueLayer(collagenLayerName[i], 1.0, true,gridWidth, gridHeight);
			this.addValueLayer(layer);
		}
		gridvaluelayerlist = new ArrayList<GridValueLayer>();

		// Create 36 fibrin bin grid value layers
		String[] fibrinLayerName = new String[36];
		GridValueLayer fibrinlayer;
		for (int i = 0; i < fiberLayerNum; i++) {
			fibrinLayerName[i] = "FibrinLayer"+i;
			fibrinlayer = new GridValueLayer(fibrinLayerName[i], 1.0, true,gridWidth, gridHeight);
			this.addValueLayer(fibrinlayer);
		}
		fibringridvaluelayerlist = new ArrayList<GridValueLayer>();
	}

/* Primary scheduled methods in order of execution */
	// Initialize and set collagen grid value layers (move to data collection)
	@ScheduledMethod(start = 0, priority = 0.5, interval = 1)
	public void CollagenValueLayer() {
		GridValueLayer colSum = (GridValueLayer) getValueLayer("Collagen Sum");
		GridValueLayer colMVA = (GridValueLayer) getValueLayer("Collagen MVA");
		GridValueLayer colMVL = (GridValueLayer) getValueLayer("Collagen MVL");
		GridValueLayer fibSum = (GridValueLayer) getValueLayer("Fibrin Sum");
		
		double fiberBin = 0;
		for (int y = 0; y < gridHeight; y++) {
			for (int x = 0; x < gridWidth; x++) {
				double tempColSum = 0;
				double tempFibSum = 0;
				double localFiberSumX = 0;
				double localFiberSumY = 0;
				int binNum = gridvaluelayerlist.size();
				for (int i = 0; i < binNum; i++) { 
					fiberBin = fiberBinSize*(i)-87.5;
					double value = gridvaluelayerlist.get(i).get(x, y);
					localFiberSumX = localFiberSumX+value*Math.cos(Math.toRadians(2*fiberBin));
					localFiberSumY = localFiberSumY+value*Math.sin(Math.toRadians(2*fiberBin));
					tempColSum = tempColSum+value;
					tempFibSum = tempFibSum+fibringridvaluelayerlist.get(i).get(x, y);
				}
				double localMeanColSumX = localFiberSumX/tempColSum;
				double localMeanColSumY = localFiberSumY/tempColSum;
				double localColMVA = 0.5*Math.toDegrees(Math.atan2(localMeanColSumY, localMeanColSumX));
				double localColMVL = Math.sqrt(Math.pow(localMeanColSumX,2)+Math.pow(localMeanColSumY,2));		
				
				// Update grid value layers			
				colSum.set(tempColSum, x, y);
				colMVA.set(localColMVA, x, y);
				colMVL.set(localColMVL, x, y);	
				fibSum.set(tempFibSum, x, y);
			}
		}
	}
	
	// Initialize collagen angle bin grid value layers
	@ScheduledMethod(start = 0, priority = 1)
	public void initializeAngleCollagenLayer() {
		
		// Set distribution across all bins based on fiber distribution type
		double[] distrib;
		if (fiberDist.equals("Circumferential")) {	// (rho = 0.75; Theta = 0)
			distrib = new double[] {.000043,.000079,.000174,.000377,.000782,.001542,.002885,.005124,.008638,
					.013821, .020988, .03025, .041382, .053729, .06621, .07744,.085965,.090574,.090574,
					.085965,.07744,.06621,.053729,.041382,.03025,.020988,.013821,.008638,.005124,.002885,
					.001542,.000782,.000377,.000174,.000079,.000043};
		} else if (fiberDist.equals("Longitudinal")) {	// (rho = 0.75; Theta = 90) 
			distrib = new double[] {.0912,.0865,.0778,.0664,.0537,.0412,.0300,.0207,.0136,.0084,.0050,.0028,
					.0015,.0007,.0004,.0002,.0001,.0000,.0000,.0001,.0002,.0004,.0007,.0015,.0028,.0050,
					.0084,.0136,.0207,.0300,.0412,.0537,.0664,.0778,.0865,.0912};
		} else {	// Uniform
			distrib = new double [] {.027778,.027778,.027778,.027778,.027778,.027778,.027778,.027778,
					.027778,.027778,.027778,.027778,.027778,.027778,.027778,.027778,.027778,.027778,.027778,
					.027778,.027778,.027778,.027778,.027778,.027778,.027778,.027778,.027778,.027778,.027778,
					.027778,.027778,.027778,.027778,.027778,.027778};
		}
		
		// Scale fiber content (from 2.5) to grid size
		double scalingFactor = Math.pow(gridUnitSize/2.5,2);
		
		// Set each grid value layer
		ArrayList<GridValueLayer> gridValueLayerList = (ArrayList<GridValueLayer>) GridValueLayerList();
		ArrayList<GridValueLayer> fibringridValueLayerList = (ArrayList<GridValueLayer>) fibrinGridValueLayerList();
		for (int i = 0; i < gridWidth; i++) {
			for (int j = 0; j < gridHeight; j++) {
				int binNum = gridValueLayerList.size();
				for (int k = 0; k < binNum; k++) {
					gridValueLayerList.get(k).set(36*scalingFactor*distrib[k], i, j);
					if (includeFib == true) {
						fibringridValueLayerList.get(k).set(25*36*scalingFactor*distrib[k], i, j);
					} else {
						fibringridValueLayerList.get(k).set(0, i, j);
					}
				}
			}
		}
	}
	
	// Initialize chemokine gradient grid value layer
	@ScheduledMethod(start = 0, priority = 3)
	public void initializeChemoLayer() throws Exception {

		// If the file with the correct gridWidth does not exist, interpolate one
		String csvFilename = "input\\ConcMean_"+gridWidth+".csv";
		File file = new File(csvFilename);
		if (!file.exists()) {
			GridInterpolator gi = new GridInterpolator();
			gi.interpolateCSV("input\\ConcMean_48.csv", gridWidth);
		}
		
		// Set chemokine value layer to CSV values
		CSVReader csvConc = new CSVReader(new FileReader(csvFilename.toString()));
		GridValueLayer chemokine = (GridValueLayer) getValueLayer("Chemokine");
		setValuelayer (csvConc, chemokine);
	}

	// Initialize grid value layers of cues
	@ScheduledMethod(start = 0, priority = 2)
	public void initializeCueLayers() throws Exception {
		String[] cueFileName = {"input\\SM_"+mechs+"_"+gridWidth+".csv", "input\\SX_"+mechs+"_"+gridWidth+".csv", 
				"input\\SY_"+mechs+"_"+gridWidth+".csv", "input\\CGM_"+gridWidth+".csv", "input\\CGX_"+gridWidth+".csv", 
				"input\\CGY_"+gridWidth+".csv"};
		GridValueLayer SM = (GridValueLayer) getValueLayer("SM");
		GridValueLayer SX = (GridValueLayer) getValueLayer("SX");	
		GridValueLayer SY = (GridValueLayer) getValueLayer("SY");
		GridValueLayer CGM = (GridValueLayer) getValueLayer("CGM");
		GridValueLayer CGX = (GridValueLayer) getValueLayer("CGX");
		GridValueLayer CGY = (GridValueLayer) getValueLayer("CGY");

		// Check that each file exists by type
		int cueNum = cueFileName.length;
		for (int i = 0; i < cueNum; i++) {
			File file = new File(cueFileName[i]);

			// If a strain files does not exist recreate all strain GVLs
			if (!file.exists() && i < 3) {
				for (int j = 0; j < gridWidth; j++) {	
					for (int k = 0; k < gridHeight; k++) {
						double [] initialStrainCue = calculateStrainCue(j, k);
						SM.set(initialStrainCue[0], k, j);
						SX.set(initialStrainCue[1], k, j);
						SY.set(initialStrainCue[2], k, j);
					}
				}
				
				// Write the new strain grid value layers to csv files for future model runs
				FileWriter SMfile = new FileWriter(cueFileName[0]);
				FileWriter SXfile = new FileWriter(cueFileName[1]);
				FileWriter SYfile = new FileWriter(cueFileName[2]);
				BufferedWriter SMwriter = new BufferedWriter(SMfile);
				BufferedWriter SXwriter = new BufferedWriter(SXfile);
				BufferedWriter SYwriter = new BufferedWriter(SYfile);
				writeValLayertoCSV(SM, SMwriter);
				writeValLayertoCSV(SX, SXwriter);
				writeValLayertoCSV(SY, SYwriter);
				i = 2; // Skip testing any remaining strain files

				// If a chemo files does not exist recreate all chemo GVLs
			} else if (!file.exists() && i >=  3) {
				for (int j = 0; j < gridWidth; j++) {
					for (int k = 0; k < gridHeight; k++) {
						double [] initialChemoCue = calculateChemoCue(j, k);
						CGM.set(initialChemoCue[0], j, k);
						CGX.set(initialChemoCue[1], j, k);
						CGY.set(initialChemoCue[2], j, k);
					}
				}

				// Write the new chemo grid value layers to csv files for future model runs
				FileWriter CGMfile = new FileWriter(cueFileName[3]);
				FileWriter CGXfile = new FileWriter(cueFileName[4]);
				FileWriter CGYfile = new FileWriter(cueFileName[5]);
				BufferedWriter CGMwriter = new BufferedWriter(CGMfile);
				BufferedWriter CGXwriter = new BufferedWriter(CGXfile);
				BufferedWriter CGYwriter = new BufferedWriter(CGYfile);
				writeValLayertoCSV(CGM, CGMwriter);
				writeValLayertoCSV(CGX, CGXwriter);
				writeValLayertoCSV(CGY, CGYwriter);
				i = 5; // Skip testing any remaining chemo files
			}
		}

		// Set value layers to CSV
		CSVReader csvSM = new CSVReader(new FileReader(cueFileName[0].toString()));
		setValuelayer (csvSM, SM);
		CSVReader csvSX = new CSVReader(new FileReader(cueFileName[1].toString()));
		setValuelayer (csvSX, SX);
		CSVReader csvSY = new CSVReader(new FileReader(cueFileName[2].toString()));
		setValuelayer (csvSY, SY);
		CSVReader csvCGM = new CSVReader(new FileReader(cueFileName[3].toString()));
		setValuelayer (csvCGM, CGM);
		CSVReader csvCGX = new CSVReader(new FileReader(cueFileName[4].toString()));
		setValuelayer (csvCGX, CGX);
		CSVReader csvCGY = new CSVReader(new FileReader(cueFileName[5].toString()));
		setValuelayer (csvCGY, CGY);
	}

	/* Secondary methods */
	// Create an array list that holds 36 grid value layers of collagen bins
	public ArrayList<GridValueLayer> GridValueLayerList() {

		int fiberLayerNum = 36;
		String [] layerName = new String[fiberLayerNum];
		for (int i = 0; i < fiberLayerNum; i++) {
			layerName[i] = "CollagenLayer"+i;
			GridValueLayer variableLayer = (GridValueLayer) getValueLayer(layerName[i]);
			if (gridvaluelayerlist.size() < fiberLayerNum) {
				gridvaluelayerlist.add(i, variableLayer);
			} else {
				gridvaluelayerlist.set(i, variableLayer);
			}
		}	
		return gridvaluelayerlist;
	}

	// An ArryList holds 36 GridValueLayer of fibrin bins (combine with gridvaluelayerlist method)
	public ArrayList<GridValueLayer> fibrinGridValueLayerList() {

		int fiberLayerNum = 36;
		String [] layerName = new String[fiberLayerNum];
		for (int i = 0; i < fiberLayerNum; i++) {
			layerName[i] = "FibrinLayer"+i;
			GridValueLayer fibrinvariableLayer = (GridValueLayer) getValueLayer(layerName[i]);
			if (fibringridvaluelayerlist.size() < fiberLayerNum) {
				fibringridvaluelayerlist.add(i, fibrinvariableLayer);
			} else {
				fibringridvaluelayerlist.set(i, fibrinvariableLayer);
			}
		}	
		return fibringridvaluelayerlist;
	}

	// Read in a CSV file and set it as grid value layer
	public void setValuelayer (CSVReader csvfile, GridValueLayer valueLayer) throws IOException {
		int j = 0;
		String[] row = null;
		while (j < gridHeight) {
			while ((row = csvfile.readNext()) !=  null) {
				for (int i = 0; i < gridWidth; i++) {
					double doubleRow = Double.parseDouble(row[i]);
					valueLayer.set(doubleRow, i, j);
				}
				j++;
			}
		}
		csvfile.close();
	}

	// Write a grid value layer to a CSV file
	public void writeValLayertoCSV(GridValueLayer gridvaluelayer, BufferedWriter writer) {
		double[][] dataArray = new double[gridHeight][gridWidth];
		StringBuilder dataStringBuilder = new StringBuilder();
		double counter = 0.0;
		for (int row = 0; row < gridHeight; row++) {
			for (int col = 0; col < gridWidth; col++) {
				counter = (double) gridvaluelayer.get(col, row);
				dataArray[row][col] = counter;
				dataStringBuilder.append(dataArray[row][col]+"");
				if (col < gridWidth)
					dataStringBuilder.append(",");
			}
			dataStringBuilder.append("\n");
		}
		try {
			writer.write(dataStringBuilder.toString());
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	// Calculate chemokine guidance cue vectors (remove /1.5)
	private double[] calculateChemoCue(int x, int y) {
		GridValueLayer chemokine = (GridValueLayer) getValueLayer("Chemokine");

		// Integrate chemokine concentration over cell perimeter
		double gradientX = 0.0;
		double gradientY = 0.0;
		List<GridPoint> edgeSites = cellPerimeter(x, y, cellRadius);
		edgeSites.add(new GridPoint(edgeSites.get(0).getX(),edgeSites.get(0).getY()));
		int numSites = edgeSites.size();
		for (int i = 1; i < numSites; i++) {
			GridPoint site1 = edgeSites.get(i-1);
			GridPoint site2 = edgeSites.get(i);
			int x1 = site1.getX();
			int y1 = site1.getY();
			int x2 = site2.getX();
			int y2 = site2.getY();
			double ang1 = Math.atan2((y-y1),(x1-x));
			double ang2 = Math.atan2((y-y2),(x2-x));
			if (ang2-ang1 < 0) {
				ang1 = ang1*-1;
			}

			// Integrate using the trapezoidal method
			double chemo1 = chemokine.get(x1, y1);
			double chemo2 = chemokine.get(x2, y2);
			double trapPart = 0.5*(ang2-ang1);
			gradientX = gradientX+trapPart*(chemo1*Math.cos(ang1)+chemo2*Math.cos(ang2));
			gradientY = gradientY+trapPart*(chemo1*Math.sin(ang1)+chemo2*Math.sin(ang2));
		}

		// Calculate chemokine gradient cue vectors
		double meanGradientX = 0.25/Math.PI*gradientX;
		double meanGradientY = 0.25/Math.PI*gradientY;
		double CGMre = Math.sqrt(Math.pow(meanGradientX, 2)+Math.pow(meanGradientY, 2)); 
		double CGXre = meanGradientX/CGMre;
		double CGYre = meanGradientY/CGMre;

		// Correct division by zero
		if(Double.isNaN(CGXre)) {
			CGXre = 0.0;
		}
		if(Double.isNaN(CGYre)) {
			CGYre = 0.0;
		}
		return new double[] {CGMre, CGXre, CGYre};
	}

	// Calculate strain guidance cue vectors
	private double[] calculateStrainCue(int x, int y) {

		// Set wound strain based on mechanics (checked)
		double eWoundXX;
		double eWoundYY;
		double eWoundXY = 0.0;
		if (mechs.equals("UniaxL")) {			// Uniaxial logitudinal condition
			eWoundXX = 0.0;
			eWoundYY = 0.05;
		} else if (mechs.equals("UniaxC")) {	// Uniaxial circumferential condition
			eWoundXX = 0.05;
			eWoundYY = 0.0;		
		} else {								// Biax condition
			eWoundXX = 0.05;
			eWoundYY = 0.05;
		}

		// Strains in remote tissue (checked)
		double eTissueXX = -0.05;
		double eTissueYY = -0.05;
		double eTissueXY = 0.0;

		// Functional border zone as fraction of average infarct diameter (checked)
		double FBZFracXX = 0.08;
		double FBZFracYY = 0.08;
		double FBZFracXY = 0.08;
		double woundRadiusX = woundRadius/gridUnitSize;
		double woundRadiusY = woundRadius/gridUnitSize;
		double productRoot = Math.sqrt(4*woundRadiusX*woundRadiusY);
		double stdvFBZXX = FBZFracXX*productRoot;
		double stdvFBZYY = FBZFracYY*productRoot;
		double stdvFBZXY = FBZFracXY*productRoot;

		// Integrate strain over cell perimeter (checked)
		double meanStrainX = 0.0;
		double meanStrainY = 0.0;
		double woundCenterX = (gridWidth-1)*0.5;
		double woundCenterY = (gridHeight-1)*0.5;
		for (int i = -180; i < 180; i++) {

			// X and Y coordinates of cell border (checked)
			double ang1 = i*Math.PI/180;
			double perimX1 = cellRadius/gridUnitSize*Math.cos(ang1)+x-woundCenterX;
			double perimY1 = cellRadius/gridUnitSize*Math.sin(ang1)+y-woundCenterY;
			
			double ang2 = (i+1)*Math.PI/180;
			double perimX2 = cellRadius/gridUnitSize*Math.cos(ang2)+x-woundCenterX;
			double perimY2 = cellRadius/gridUnitSize*Math.sin(ang2)+y-woundCenterY;

			// Compute strain (checked)
			double perimR1 = Math.sqrt(Math.pow(perimX1, 2)+Math.pow(perimY1, 2));
			double perimT1 = Math.atan2(perimY1,perimX1);
			double woundBorderR1 = woundRadiusX*woundRadiusY/(Math.sqrt(Math.pow((woundRadiusX*Math.sin(perimT1)),2)
					+Math.pow((woundRadiusY*Math.cos(perimT1)),2)));
			double errorInput1 = -(perimR1-woundBorderR1)/Math.sqrt(2);
			double strainXX1 = eWoundXX+0.5*(eTissueXX-eWoundXX)*Erf.erfc(errorInput1/stdvFBZXX);
			double strainYY1 = eWoundYY+0.5*(eTissueYY-eWoundYY)*Erf.erfc(errorInput1/stdvFBZYY);
			double strainXY1 = eWoundXY+0.5*(eTissueXY-eWoundXY)*Erf.erfc(errorInput1/stdvFBZXY);
			
			double perimR2 = Math.sqrt(Math.pow(perimX2, 2)+Math.pow(perimY2, 2));
			double perimT2 = Math.atan2(perimY2,perimX2);
			double woundBorderR2 = woundRadiusX*woundRadiusY/(Math.sqrt(Math.pow((woundRadiusX*Math.sin(perimT2)),2)
					+Math.pow((woundRadiusY*Math.cos(perimT2)),2)));
			double errorInput2 = -(perimR2-woundBorderR2)/Math.sqrt(2);
			double strainXX2 = eWoundXX+0.5*(eTissueXX-eWoundXX)*Erf.erfc(errorInput2/stdvFBZXX);
			double strainYY2 = eWoundYY+0.5*(eTissueYY-eWoundYY)*Erf.erfc(errorInput2/stdvFBZYY);
			double strainXY2 = eWoundXY+0.5*(eTissueXY-eWoundXY)*Erf.erfc(errorInput2/stdvFBZXY);

			// Calculate cell's normal strain (checked)
			double normX1 = Math.cos(ang1);
			double normY1 = Math.sin(ang1);
			double normStrain1 = (strainXX1*normX1+strainXY1*normY1)*normX1+(strainXY1*normX1+strainYY1*normY1)*normY1;
			
			double normX2 = Math.cos(ang2);
			double normY2 = Math.sin(ang2);
			double normStrain2 = (strainXX2*normX2+strainXY2*normY2)*normX2+(strainXY2*normX2+strainYY2*normY2)*normY2;



			// Calculate strain components by integration using the trapezoidal method (checked)
			double trapPart = 0.25*(ang2-ang1)/Math.PI;
			meanStrainX = meanStrainX+trapPart*(normStrain2*Math.cos(2*ang2)+normStrain1*Math.cos(2*ang1));
			meanStrainY = meanStrainY+trapPart*(normStrain2*Math.sin(2*ang2)+normStrain1*Math.sin(2*ang1));
			
		}
		// Strain components
		double meanStrainAngle = 0.5*Math.atan2(meanStrainY, meanStrainX);
		double SXre = Math.cos(meanStrainAngle);
		double SYre = Math.sin(meanStrainAngle);
		double SMre = Math.sqrt(Math.pow(meanStrainX, 2)+Math.pow(meanStrainY, 2));

		return new double[] {SMre, SXre, SYre};
	}

	// Returns a list containing all the grids within one half grid from the perimeter
	private List<GridPoint> cellPerimeter(int x, int y, double radius) {

		// Build list of neighbors
		int extent = (int) Math.floor(radius/gridUnitSize);
		if (extent == 0) {
			extent = 1;
		}

		// Build perimeter list from neighbor list by checking the distance
		List<GridPoint> edgeSites = new ArrayList<GridPoint>();
		for (int i = -extent; i <= 0; i++) {	// Quadrant 3
			for (int j = 0; j <= extent; j++) {
				int xCoor = x+i;
				int yCoor = y+j;
				GridPoint point = new GridPoint(xCoor,yCoor);
				double dist = Math.sqrt(Math.pow(xCoor-x,2)+Math.pow(yCoor-y,2))*gridUnitSize;
				if (dist <= radius+gridUnitSize && dist >= radius) {
					edgeSites.add(point);
				}
			}
		}
		for (int i = 0; i <= extent ; i++) {	// Quadrant 4
			for (int j = extent; j >= 0; j--) {
				int xCoor = x+i;
				int yCoor = y+j;
				GridPoint point = new GridPoint(xCoor,yCoor);
				if (!edgeSites.contains(point)) {
					double dist = Math.sqrt(Math.pow(xCoor-x,2)+Math.pow(yCoor-y,2))*gridUnitSize;
					if (dist <= radius+gridUnitSize && dist >= radius) {
						edgeSites.add(point);
					}
				}
			}
		}
		for (int i = extent; i >= 0 ; i--) {	// Quadrant 1
			for (int j = 0; j >= -extent; j--) {
				int xCoor = x+i;
				int yCoor = y+j;
				GridPoint point = new GridPoint(xCoor,yCoor);
				if (!edgeSites.contains(point)) {
					double dist = Math.sqrt(Math.pow(xCoor-x,2)+Math.pow(yCoor-y,2))*gridUnitSize;
					if (dist <= radius+gridUnitSize && dist >= radius) {
						edgeSites.add(point);
					}
				}
			}
		}
		for (int i = 0; i >= -extent ; i--) {	// Quadrant 2
			for (int j = -extent; j <= 0; j++) {
				int xCoor = x+i;
				int yCoor = y+j;
				GridPoint point = new GridPoint(xCoor,yCoor);
				if (!edgeSites.contains(point)) {
					double dist = Math.sqrt(Math.pow(xCoor-x,2)+Math.pow(yCoor-y,2))*gridUnitSize;
					if (dist <= radius+gridUnitSize && dist >= radius) {
						edgeSites.add(point);
					}
				}
			}
		}

		return edgeSites;
	}

	// Returns the grid size to other classes
	public static double getGridSize() {
		return gridUnitSize;
	}
	
	// Returns the model dimensions
	public static int[] getModelDimension() {
		return new int[] {gridHeight, gridWidth};
	}
	
/* Data collection methods */
	// Writes model parameters to console and to csv file (change scheduling to count fibroblasts before infarction)
	@SuppressWarnings("unchecked")
	@ScheduledMethod(start = 0, priority = 4)
	public void writeParameters() throws IOException {
		Grid<?> grid = (Grid<?>) getProjection("Cell Grid");
		Iterable<CellAgentSim> cells = (Iterable<CellAgentSim>) grid.getObjects();
		int cellSum = (int) cells.spliterator().getExactSizeIfKnown();
		
		// Build data strings
		String heading = "GridSize,Width,Height,IntCells,Mechanics,IntFiberDist,Fibrin";
		String output = Double.toString(gridUnitSize)+","+Integer.toString(sampleWidth)+","+
				Integer.toString(sampleHeight)+","+Integer.toString(cellSum)+","+mechs+","+
				fiberDist+","+Boolean.toString(includeFib);
		String[][] fibroblastParameters = CellAgentSim.getFibroblastParameters();
		
		StringBuilder paramBuilder = new StringBuilder();
		paramBuilder.append(heading+","+fibroblastParameters[0][0]+"\n");
		paramBuilder.append(output+","+fibroblastParameters[0][1]+"\n");
		
		// Write data to console
		System.out.println(paramBuilder);
		
		// Write data to CSV
		FileWriter csvWoundParameters = new FileWriter(woundParamFilename);
		try {
			BufferedWriter woundStatWriter = new BufferedWriter(csvWoundParameters);
			woundStatWriter.write(paramBuilder.toString());
			woundStatWriter.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	// Calculate cell and collagen MVA, MVL, and fraction in the wound (cut iteration through 36 gridvalue layer lists)
	@ScheduledMethod(start = 0.5, interval = 24, priority = 2)
	public void writeCollagenSTAT() throws IOException {
		Grid<?> grid = (Grid<?>) getProjection("Cell Grid");
		GridValueLayer colSum = (GridValueLayer) getValueLayer("Collagen Sum");
		GridValueLayer fibSum = (GridValueLayer) getValueLayer("Fibrin Sum");

		double woundGridNum = (double) (woundRadius/gridUnitSize);
		double wCenter = gridWidth/2;
		double hCenter = gridHeight/2;

		// Collect collagen and fibrin in the wound area
		double woundColSum = 0;
		double woundFibSum = 0;
		int binNum = gridvaluelayerlist.size();
		double [] binColSum = new double[binNum];
		double [] binFibSum = new double[binNum];
		for (int y = 0; y < gridHeight; y++) {
			for (int x = 0; x < gridWidth; x++) {
				if (Math.sqrt(Math.pow((x-wCenter), 2)+Math.pow((y-hCenter), 2)) <= woundGridNum) {
					woundColSum = woundColSum+colSum.get(x, y);
					woundFibSum = woundFibSum+fibSum.get(x, y);
					for (int i = 0; i < binNum; i++) {
						binColSum[i] = binColSum[i]+gridvaluelayerlist.get(i).get(x, y);	
						binFibSum[i] = binFibSum[i]+fibringridvaluelayerlist.get(i).get(x, y);
					}
				}
			}
		}

		// Calculate collagen content, MVA, and MVL
		double woundArea = Math.PI*Math.pow(woundRadius, 2);
		double colFrac = woundColSum/(woundArea*fiberDensity);
		double colSumX = 0;
		double colSumY = 0;
		for (int i = 0; i < binNum; i++) {	
			double fiberBin = fiberBinSize*(i)-87.5;
			double fiberBinRad = Math.toRadians(2*fiberBin);
			colSumX = colSumX+binColSum[i]*Math.cos(fiberBinRad);
			colSumY = colSumY+binColSum[i]*Math.sin(fiberBinRad);
		}
		double meanColSumX = colSumX/woundColSum;
		double meanColSumY = colSumY/woundColSum;
		double colMVA = 0.5*Math.toDegrees(Math.atan2(meanColSumY, meanColSumX));
		double colMVL = Math.sqrt(Math.pow(meanColSumX, 2)+Math.pow(meanColSumY, 2));

		// Calculate fibrin content, MVA, and MVL
		double fibFrac = woundFibSum/(woundArea*fiberDensity);
		double fibSumX = 0;
		double fibSumY = 0;
		for (int i = 0; i < binNum; i++) {	
			double fiberBin = fiberBinSize*(i)-87.5;
			double fiberBinRad = Math.toRadians(2*fiberBin);
			fibSumX = fibSumX+binFibSum[i]*Math.cos(fiberBinRad);
			fibSumY = fibSumY+binFibSum[i]*Math.sin(fiberBinRad);
		}
		double meanFibSumX = fibSumX/woundFibSum;
		double meanFibSumY = fibSumY/woundFibSum;
		double fibMVA = 0.5*Math.toDegrees(Math.atan2(meanFibSumY, meanFibSumX));
		double fibMVL = Math.sqrt(Math.pow(meanFibSumX, 2)+Math.pow(meanFibSumY, 2));

		// Collect cells in the wound area
		int cellSum = 0;
		double cellSumX = 0;
		double cellSumY = 0;
		@SuppressWarnings("unchecked")
		Iterable<CellAgentSim> cells = (Iterable<CellAgentSim>) grid.getObjects();
		for (CellAgentSim cell : cells) {
			GridPoint pt = grid.getLocation(cell);
			if (Math.sqrt(Math.pow((pt.getX()-wCenter), 2)+Math.pow((pt.getY()-hCenter), 2)) <= woundGridNum) {
				double cellAngle = cell.angleSelection;
				double cellAngleRad = Math.toRadians(2*cellAngle);
				cellSumX = cellSumX+Math.cos(cellAngleRad);
				cellSumY = cellSumY+Math.sin(cellAngleRad);
				cellSum = cellSum+1;
			}
		}

		// Calculate cell content, MVA, and MVL
		double cellArea=Math.pow(cellRadius, 2)*Math.PI;
		double cellFrac = cellSum*cellArea/woundArea;
		double meanCellSumX = cellSumX/cellSum;
		double meanCellSumY = cellSumY/cellSum; 
		double cellMVA = 0.5*Math.toDegrees(Math.atan2(meanCellSumY, meanCellSumX));
		double cellMVL = Math.sqrt(Math.pow(meanCellSumX, 2)+Math.pow(meanCellSumY, 2));

		// Retrieve tick time
		ISchedule schedule = RunEnvironment.getInstance().getCurrentSchedule();
		double time = schedule.getTickCount();

		// Output collagen and cell content, MVA, and MVL
		String[] heading = {"Time", "ColFRC","ColMVA", "ColMVL", "CellFRC", "CellMVA", "CellMVL", "FibFRC", "FibMVA", "FibMVL"};
		double[] output = {time, colFrac, colMVA, colMVL, cellFrac, cellMVA, cellMVL, fibFrac, fibMVA, fibMVL};
		int outputNum = output.length;

		// Initiate CSV file writer
		FileWriter csvWoundStat;
		FileWriter csvWoundColAngDist;
		if (time < 1) {
			csvWoundStat = new FileWriter(woundStatFilename);
			csvWoundColAngDist = new FileWriter(woundColAnDistFilename);
		} else {
			csvWoundStat = new FileWriter(woundStatFilename, true);
			csvWoundColAngDist = new FileWriter(woundColAnDistFilename, true);
		}

		// Build data strings
		StringBuilder woundStatBuilder = new StringBuilder();
		if (time < 1) {
			for (int i = 0; i < outputNum; i++) {
				if (i == outputNum-1) {
					woundStatBuilder.append(heading[i]);
				}
				else {
					woundStatBuilder.append(heading[i]+",");
				}
			}
			woundStatBuilder.append("\n");
		}
		for (int i = 0; i < outputNum; i++) {
			if (i == outputNum-1) {
				woundStatBuilder.append(output[i]);
			}
			else {
				woundStatBuilder.append(output[i]+",");
			}
		}
		woundStatBuilder.append("\n");

		StringBuilder woundColAngDistBuilder = new StringBuilder();
		woundColAngDistBuilder.append(time+",");
		for (int i = 0; i < binNum; i++) {
			if (i == binNum-1) {
				woundColAngDistBuilder.append(binColSum[i]);
			}
			else {
				woundColAngDistBuilder.append(binColSum[i]+",");
			}
		}
		woundColAngDistBuilder.append("\n");	

		// Write data to CSV
		try {
			BufferedWriter woundStatWriter = new BufferedWriter(csvWoundStat);
			woundStatWriter.write(woundStatBuilder.toString());
			woundStatWriter.close();
			BufferedWriter woundColAngDistWriter = new BufferedWriter(csvWoundColAngDist);
			woundColAngDistWriter.write(woundColAngDistBuilder.toString());
			woundColAngDistWriter.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}