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
import repast.simphony.parameter.Parameters;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridBuilderParameters;
import repast.simphony.space.grid.GridPoint;
import repast.simphony.space.grid.RandomGridAdder;
import repast.simphony.space.grid.WrapAroundBorders;
import repast.simphony.valueLayer.GridValueLayer;

public class woundABMContextSim extends DefaultContext<Object> {
	
	// Pull parameters
	Parameters p = RunEnvironment.getInstance().getParameters();
	
	// Fibroblast parameters
	private double Wp = (Double) p.getValue("Wp");	// Persistance cue weight
	private double Ws = (Double) p.getValue("Ws"); 	// Structural cue weight
	private double Wm = (Double) p.getValue("Wm");	// Mechanics cue weight
	private double Wc = (Double) p.getValue("Wc");	// Chemokine cue weight
	private String depoType = (String) p.getValue("depoType"); 	// "Aligned" or "Random"
	private double gMitosisTime = (Double) p.getValue("mitosisTime");
	private double gApoptosisTime = (Double) p.getValue("apoptosisTime");
	private double gDepositionTime = (Double) p.getValue("depositionTime");
	private double gDegradationTime = (Double) p.getValue("degradationTime");
	private boolean cellRemoval = (Boolean) p.getValue("cellRemoval");
	private boolean colRotation = (Boolean) p.getValue("colRotation");
	
	// Geometry
	double sampleWidth = (Double) p.getValue("sampleWidth");			// um
	double sampleHeight = (Double) p.getValue("sampleHeight");			// um
	double gridUnitSize = (Double) p.getValue("gridUnitSize");		// um
	int gridWidth = (int) Math.floor(sampleWidth/gridUnitSize);		// grids
	int gridHeight = (int) Math.floor(sampleHeight/gridUnitSize);	// grids
	private double cellRadius = 5.0;										// um
	private double woundRadius = 113.0; 									// um
	
	// Value layer information
	ArrayList<ArrayList<double[]>> initialFiberArrayList;
	private ArrayList<GridValueLayer> gridvaluelayerlist;
	private ArrayList<GridValueLayer> fibringridvaluelayerlist;
	private String[] basicValueLayer = {"Collagen Sum", "Fibrin Sum"};
	private String[] cueValueLayer = {"CM", "CX", "CY", "MM", "MX", "MY"};
	private String mechs = (String) p.getValue("mechs");			// "UniaxC" "UniaxL" or "Biax"
	private double fiberDensity = 178.25;							// fibers/um^2/7um thickness
	private int binSize = 5;									// degrees
	private String initialFiberDist = (String) p.getValue("initialFiberDist");	// "Circumferential" "Longitudinal" or "Uniform"
	private boolean includeFibrin = (Boolean) p.getValue("includeFibrin");		// true or false
	private double scalePercent = (double) p.getValue("initialColPercent")/3.0;		// 3% normally
	
	// Output file specification
	private String outputDir = "output";
	private Date date = new Date();
	private DateFormat df = new SimpleDateFormat("yyyyMMdd_HHmmss");
	private String modelTag = df.format(date);
	private String woundParamFilename = outputDir+"\\"+mechs+"_"+gridUnitSize+"_"+modelTag+"_WoundParameters.csv";
	private String woundStatFilename = outputDir+"\\"+mechs+"_"+gridUnitSize+"_"+modelTag+"_WoundStat.csv";
	private String woundColAnDistFilename = outputDir+"\\"+mechs+"_"+gridUnitSize+"_"+modelTag+"_WoundColFiberAngDist.csv";
	
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
	// Initialize and set collagen grid value layers
	@ScheduledMethod(start = 0, priority = 0.5, interval = 0.5)
	public void CollagenValueLayer() {
		GridValueLayer colSum = (GridValueLayer) getValueLayer("Collagen Sum");
		GridValueLayer fibSum = (GridValueLayer) getValueLayer("Fibrin Sum");
		
		for (int y = 0; y < gridHeight; y++) {
			for (int x = 0; x < gridWidth; x++) {
				double tempColSum = 0;
				double tempFibSum = 0;
				int binNum = gridvaluelayerlist.size();
				for (int i = 0; i < binNum; i++) {
					double value = gridvaluelayerlist.get(i).get(x, y);
					tempColSum = tempColSum+value;
					tempFibSum = tempFibSum+fibringridvaluelayerlist.get(i).get(x, y);
				}
				
				// Update grid value layers			
				colSum.set(tempColSum, x, y);	
				fibSum.set(tempFibSum, x, y);
			}
		}
	}
	
	// Initialize collagen angle bin grid value layers
	@ScheduledMethod(start = 0, priority = 1)
	public void initializeAngleCollagenLayer() {
		
		// Set distribution across all bins based on fiber distribution type
		double[] distrib;
		if (initialFiberDist.equals("Circumferential")) {	// (rho = 0.75; Theta = 0)
			distrib = new double[] {.000043,.000079,.000174,.000377,.000782,.001542,.002885,.005124,.008638,
					.013821, .020988, .03025, .041382, .053729, .06621, .07744,.085965,.090574,.090574,
					.085965,.07744,.06621,.053729,.041382,.03025,.020988,.013821,.008638,.005124,.002885,
					.001542,.000782,.000377,.000174,.000079,.000043};
		} else if (initialFiberDist.equals("Longitudinal")) {	// (rho = 0.75; Theta = 90) 
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
					gridValueLayerList.get(k).set(scalePercent*36*scalingFactor*distrib[k], i, j);
					if (includeFibrin == true) {
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
		GridValueLayer chemokine = (GridValueLayer) getValueLayer("Chemokine");
		
		// Check if input file is in current directory (added for batch runs)
		String path = "C:\\Users\\acb4ft\\Desktop\\WoundABM\\woundABMSim\\input";
		
		// If the file with the correct gridWidth does not exist, create one
		String csvFilename = path+"\\ConcMean_"+gridWidth+".csv";
		File file = new File(csvFilename);
		if (!file.exists()) {
			System.out.println("ConcMean file not found");
		} else {	// Set chemokine value layer to CSV values
			CSVReader csvConc = new CSVReader(new FileReader(csvFilename.toString()));
			setValuelayer(csvConc, chemokine);
		}
	}

	// Initialize grid value layers of cues
	@ScheduledMethod(start = 0, priority = 2)
	public void initializeCueLayers() throws Exception {
		
		// Check if input file is in current directory (added for batch runs)
		String path;
		File f = new File("input");
		if (!f.exists()) {	// change path to your directory
			path = "C:\\Users\\acb4ft\\Desktop\\WoundABM\\woundABMSim\\input";
		} 
		else {
			path = "input";
		}
		
		String[] cueFileName = {path+"\\MM_"+mechs+"_"+gridWidth+".csv", path+"\\MX_"+mechs+"_"+
				gridWidth+".csv", path+"\\MY_"+mechs+"_"+gridWidth+".csv", path+"\\CM_"+gridWidth+".csv", 
				path+"\\CX_"+gridWidth+".csv", path+"\\CY_"+gridWidth+".csv"};
		GridValueLayer MM = (GridValueLayer) getValueLayer("MM");
		GridValueLayer MX = (GridValueLayer) getValueLayer("MX");	
		GridValueLayer MY = (GridValueLayer) getValueLayer("MY");
		GridValueLayer CM = (GridValueLayer) getValueLayer("CM");
		GridValueLayer CX = (GridValueLayer) getValueLayer("CX");
		GridValueLayer CY = (GridValueLayer) getValueLayer("CY");

		// Check that each file exists by type
		int cueNum = cueFileName.length;
		for (int i = 0; i < cueNum; i++) {
			File file = new File(cueFileName[i]);

			// If a strain files does not exist recreate all strain GVLs
			if (!file.exists() && i < 3) {
				for (int j = 0; j < gridWidth; j++) {	
					for (int k = 0; k < gridHeight; k++) {
						double [] initialStrainCue = calculateStrainCue(j, k);
						MM.set(initialStrainCue[0], k, j);
						MX.set(initialStrainCue[1], k, j);
						MY.set(initialStrainCue[2], k, j);
					}
				}
				
				// Write the new strain grid value layers to csv files for future model runs
				FileWriter MMfile = new FileWriter(cueFileName[0]);
				FileWriter MXfile = new FileWriter(cueFileName[1]);
				FileWriter MYfile = new FileWriter(cueFileName[2]);
				BufferedWriter MMwriter = new BufferedWriter(MMfile);
				BufferedWriter MXwriter = new BufferedWriter(MXfile);
				BufferedWriter MYwriter = new BufferedWriter(MYfile);
				writeValLayertoCSV(MM, MMwriter);
				writeValLayertoCSV(MX, MXwriter);
				writeValLayertoCSV(MY, MYwriter);
				i = 2; // Skip testing any remaining strain files

				// If a chemo files does not exist recreate all chemo GVLs
			} else if (!file.exists() && i >=  3) {
				for (int j = 0; j < gridWidth; j++) {
					for (int k = 0; k < gridHeight; k++) {
						double [] initialChemoCue = calculateChemoCue(j, k);
						CM.set(initialChemoCue[0], j, k);
						CX.set(initialChemoCue[1], j, k);
						CY.set(initialChemoCue[2], j, k);
					}
				}

				// Write the new chemo grid value layers to csv files for future model runs
				FileWriter CMfile = new FileWriter(cueFileName[3]);
				FileWriter CXfile = new FileWriter(cueFileName[4]);
				FileWriter CYfile = new FileWriter(cueFileName[5]);
				BufferedWriter CMwriter = new BufferedWriter(CMfile);
				BufferedWriter CXwriter = new BufferedWriter(CXfile);
				BufferedWriter CYwriter = new BufferedWriter(CYfile);
				writeValLayertoCSV(CM, CMwriter);
				writeValLayertoCSV(CX, CXwriter);
				writeValLayertoCSV(CY, CYwriter);
				i = 5; // Skip testing any remaining chemo files
			}
		}

		// Set value layers to CSV
		CSVReader csvMM = new CSVReader(new FileReader(cueFileName[0].toString()));
		setValuelayer (csvMM, MM);
		CSVReader csvMX = new CSVReader(new FileReader(cueFileName[1].toString()));
		setValuelayer (csvMX, MX);
		CSVReader csvMY = new CSVReader(new FileReader(cueFileName[2].toString()));
		setValuelayer (csvMY, MY);
		CSVReader csvCM = new CSVReader(new FileReader(cueFileName[3].toString()));
		setValuelayer (csvCM, CM);
		CSVReader csvCX = new CSVReader(new FileReader(cueFileName[4].toString()));
		setValuelayer (csvCX, CX);
		CSVReader csvCY = new CSVReader(new FileReader(cueFileName[5].toString()));
		setValuelayer (csvCY, CY);
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

	// An ArryList holds 36 GridValueLayer of fibrin bins
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

		String[] row = null;
		for (int j = 0; j < gridHeight; j++) {
			row = csvfile.readNext();
			for (int i = 0; i < gridWidth; i++) {
				double doubleRow = Double.parseDouble(row[i]);
				valueLayer.set(doubleRow, i, j);
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

	// Calculate chemokine guidance cue vectors
	private double[] calculateChemoCue(int x, int y) {
		GridValueLayer chemokine = (GridValueLayer) getValueLayer("Chemokine");

		// Integrate chemokine concentration over cell perimeter
		double gradientX = 0.0;
		double gradientY = 0.0;
		double cellGridRadius = cellRadius/gridUnitSize;
		List<double[]> data = perimeterData(x, y, cellGridRadius, chemokine);
		int num = data.size();
		for (int i = 1; i < num; i++) {
			double[] data1 = data.get(i-1);
			double[] data2 = data.get(i);
			double ang1 = data1[0];
			double ang2 = data2[0];

			// Integrate using the trapezoidal method
			double chemo1 = data1[3];
			double chemo2 = data2[3];
			double trapPart = 0.5*(ang2-ang1);
			gradientX = gradientX+trapPart*(chemo1*Math.cos(ang1)+chemo2*Math.cos(ang2));
			gradientY = gradientY+trapPart*(chemo1*Math.sin(ang1)+chemo2*Math.sin(ang2));
		}

		// Calculate chemokine gradient cue vectors
		double meanGradientX = 0.5/Math.PI*gradientX;
		double meanGradientY = 0.5/Math.PI*gradientY;
		double CMre = Math.sqrt(Math.pow(meanGradientX, 2)+Math.pow(meanGradientY, 2)); 
		double CXre = meanGradientX/CMre;
		double CYre = meanGradientY/CMre;

		// Correct division by zero
		if(Double.isNaN(CXre)) {
			CXre = 0.0;
		}
		if(Double.isNaN(CYre)) {
			CYre = 0.0;
		}
		return new double[] {CMre, CXre, CYre};
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
		double MXre = Math.cos(meanStrainAngle);
		double MYre = Math.sin(meanStrainAngle);
		double MMre = Math.sqrt(Math.pow(meanStrainX, 2)+Math.pow(meanStrainY, 2));

		return new double[] {MMre, MXre, MYre};
	}

	// Returns data at the perimeter of the cell for the specified grid value layer 
	private List<double[]> perimeterData(int x, int y, double radius, GridValueLayer layer) {

		// Build perimeter list from neighbor list
		List<double[]> data = new ArrayList<double[]>();
		for (int angle = 0; angle <= 360; angle++) {
			double angleRad = Math.toRadians(angle);
			double xCoor = x+radius*Math.cos(angleRad);
			double yCoor = y+radius*Math.sin(angleRad);

			// Determine surrounding coordinates
			int x1 = (int) Math.floor(xCoor);
			int x2 = (int) Math.ceil(xCoor);
			int y1 = (int) Math.floor(yCoor);
			int y2 = (int) Math.ceil(yCoor);
			
			// Wrapped space (only applies to layer calls)
			int x1w = wrapCoordinate(x1, gridWidth);
			int y1w = wrapCoordinate(y1, gridHeight);
			
			// Bilinearly interpolate the data for the current x and y coordinates
			if (x1 == xCoor && y1 == yCoor) {
				double fxy = layer.get(x1w,y1w);
				double[] value = new double[] {angleRad,xCoor,yCoor,fxy};
				data.add(value);
			} else if (x1 == xCoor) {
				int y2w = wrapCoordinate(y2, gridHeight);
				double fxy = (yCoor-y1)*(layer.get(x1w, y2w)-layer.get(x1w, y1w))/(y2-y1)+layer.get(x1w, y1w);
				double[] value = new double[] {angleRad,xCoor,yCoor,fxy};
				data.add(value);
			} else if (y1 == yCoor) {
				int x2w = wrapCoordinate(x2, gridWidth);
				double fxy = (xCoor-x1)*(layer.get(x2w, y1w)-layer.get(x1w, y1w))/(x2-x1)+layer.get(x1w, y1w);
				double[] value = new double[] {angleRad,xCoor,yCoor,fxy};
				data.add(value);
			} else {
				int x2w = wrapCoordinate(x2, gridWidth);
				int y2w = wrapCoordinate(y2, gridHeight);
				double xDif = (x2-x1);
				double yDif = (y2-y1);
				double xRatio1 = (x2-xCoor)/xDif;
				double xRatio2 = (xCoor-x1)/xDif;
				double fxy1 = xRatio1*layer.get(x1w, y1w)+xRatio2*layer.get(x2w, y1w);
				double fxy2 = xRatio1*layer.get(x1w, y2w)+xRatio2*layer.get(x2w, y2w);
				double fxy = (y2-yCoor)/yDif*fxy1+(yCoor-y1)/yDif*fxy2;
				double[] value = new double[] {angleRad,xCoor,yCoor,fxy};
				data.add(value);
			}
		}

		return data;
	}
	
	// Checks if a coordinate falls outside of the model space and wraps it if it does
	public int wrapCoordinate(int oldCoor, int dimLength) {
		
		// Check and wrap coordinate
		int newCoor = oldCoor;
		if (oldCoor < 0) {
			newCoor = oldCoor+dimLength;
		} else if (oldCoor > dimLength) {
			newCoor = oldCoor-dimLength;
		}
		
		return newCoor;
	}
	
/* Data collection methods */
	// Writes model parameters to console and to csv file
	@SuppressWarnings("unchecked")
	@ScheduledMethod(start = 0, priority = 4)
	public void writeParameters() throws IOException {
		Grid<?> grid = (Grid<?>) getProjection("Cell Grid");
		Iterable<CellAgentSim> cells = (Iterable<CellAgentSim>) grid.getObjects();
		int cellSum = (int) cells.spliterator().getExactSizeIfKnown();
		
		// Make output file if it does not already exist
		File directory = new File(outputDir);
		if (!directory.exists()) {
			directory.mkdir();
		}
		
		// Build data strings
		String heading = "GridSize,Width,Height,IntCells,Mechanics,IntFiberDist,Fibrin,DepoType,CollagenRot,..."
				+ "Infarction,MitosisTime,ApoptosisTime,DepoTime,DegTime,Wp,Wc,Ws,Wm";
		String output = Double.toString(gridUnitSize)+","+Double.toString(sampleWidth)+","+
				Double.toString(sampleHeight)+","+Integer.toString(cellSum)+","+mechs+","+
				initialFiberDist+","+Boolean.toString(includeFibrin)+","+depoType+","+
				Boolean.toString(colRotation)+","+Boolean.toString(cellRemoval)+","+
				Double.toString(gMitosisTime)+","+Double.toString(gApoptosisTime)+","+
				Double.toString(gDepositionTime)+","+Double.toString(gDegradationTime)+","+
				Double.toString(Wp)+","+Double.toString(Wc)+","+Double.toString(Ws)+","+
				Double.toString(Wm);
		
		StringBuilder paramBuilder = new StringBuilder();
		paramBuilder.append(heading+"\n");
		paramBuilder.append(output+"\n");
		
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
	
	// Calculate cell and collagen MVA, MVL, and fraction in the wound
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
			double fiberBin = 2*Math.toRadians(binSize*(i)-87.5);
			colSumX = colSumX+binColSum[i]*Math.cos(fiberBin);
			colSumY = colSumY+binColSum[i]*Math.sin(fiberBin);
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
			double fiberBin = 2*Math.toRadians(binSize*(i)-87.5);
			fibSumX = fibSumX+binFibSum[i]*Math.cos(fiberBin);
			fibSumY = fibSumY+binFibSum[i]*Math.sin(fiberBin);
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
				double cellAngle = 2*Math.toRadians(cell.angleSelection);
				cellSumX = cellSumX+Math.cos(cellAngle);
				cellSumY = cellSumY+Math.sin(cellAngle);
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