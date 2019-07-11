package woundABMSim;

import java.awt.Color;
import java.util.HashMap;
import java.util.Map;

import repast.simphony.valueLayer.ValueLayer;
import repast.simphony.visualizationOGL2D.ValueLayerStyleOGL;

public class ChemokineStyleSim implements ValueLayerStyleOGL {

	private ValueLayer chemokine;
	Map<Integer, Color> colorMap = new HashMap<Integer,Color>();
	
	// Create color gradient
	public ChemokineStyleSim() {
		for (int i = 0; i < 100; i++) {
			colorMap.put(i, new Color(128, (int) (128* i /99), 128));
		}
	}

	// Translate value layer to colors for display
	public Color getColor(double... coordinates) {
		double val = (double) chemokine.get(coordinates)*100;
		Color color = colorMap.get((int) val);
		
		if (val <= 100 && val > 10) {
			color = Color.white;
		} else if (color == null) {
			color = Color.blue;
		} else {
			color = Color.gray;
		}
		return color;
	}
	
	// Specify grid cell size
	public float getCellSize() {
		return 15.0f;
	}
	
	// Get chemokine value layer ("Chemokine" specified in xml)
	public void init(ValueLayer chemokine) {
		this.chemokine = chemokine;
	}
}