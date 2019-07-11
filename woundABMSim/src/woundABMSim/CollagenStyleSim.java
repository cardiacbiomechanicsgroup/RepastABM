package woundABMSim;

import java.awt.Color;
import java.util.HashMap;
import java.util.Map;

import repast.simphony.valueLayer.ValueLayer;
import repast.simphony.visualizationOGL2D.ValueLayerStyleOGL;

public class CollagenStyleSim implements ValueLayerStyleOGL {
	double gridUnitSize = woundABMContextSim.getGridSize();
	private ValueLayer collagen;
	private Map<Integer, Color> colorMap = new HashMap<Integer, Color>();
	
	// Create color gradient
	public CollagenStyleSim() {
		for (int i = 0; i < 10; i++) {
			colorMap.put(i, new Color(i / 10.0f, 0f, 0f));
		}
	}
	
	// Translate value layer to colors for display
	public Color getColor(double... coordinates) {
		int val = (int) Math.round(collagen.get(coordinates)*Math.pow((10/gridUnitSize),2)/700);
		Color color = colorMap.get(val);
		
		if (val <= 10) {
			return color;
		} else if (color == null) {
			color = Color.red;
		} else {
			color = Color.black;
		}
		return color;
	}

	// Specify grid cell size
	public float getCellSize() {
		return 15.0f;
	}
	
	// Get collagen value layer ("Collagen Sum" specified in xml)
	public void init(ValueLayer colSum) {
		this.collagen = colSum;
	}
}