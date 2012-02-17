/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

package xmipp.viewer;

import java.util.ArrayList;

import ij.ImagePlus;
import xmipp.ij.XmippImageConverter;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.utils.DEBUG;

public class MetadataGallery extends ImageGallery {
	private static final long serialVersionUID = 1L;

	// Label to be rendered
	protected int renderLabel;
	boolean renderLabels;
	protected int displayLabel;
	
	// Also store the visible ones to fast access
	ArrayList<ColumnInfo> visibleLabels;
	
//	int[] labels;
//	String[] labelStr;

	// Ids of objects stored in metadata
	long[] ids;

	public MetadataGallery(GalleryData data) throws Exception {
		super(data);
	}
	
	/** Update the columns display information */
	public void updateColumnInfo(ArrayList<ColumnInfo> newInfo){
		int n = newInfo.size();
		boolean changed = false;
		renderLabels = false;
		visibleLabels.clear();
		for (int i = 0; i < n; ++i){
			ColumnInfo ci1 = data.labels.get(i);
			ColumnInfo ci2 = newInfo.get(i);
			if (ci1.visible != ci2.visible){
				ci1.visible = ci2.visible;
				changed = true;
			}
			if (ci1.render != ci2.render){
				ci1.render = ci2.render;
				changed = true;
			}
			if (ci1.visible) 
				visibleLabels.add(ci1);
			if (ci1.render)
				renderLabels = true;
		}
		if (changed) {
			calculateCellSize();
			cols = visibleLabels.size();
			//	fireTableDataChanged();
			fireTableStructureChanged();
		}
	}

	// Load initial dimensions
	protected ImageDimension loadDimension() throws Exception {
		//md = new MetaData(filename);
		ids = data.md.findObjects();
		renderLabel = MDLabel.MDL_IMAGE;
		displayLabel = MDLabel.MDL_IMAGE;
		ImageGeneric image = getImage(0, renderLabel);
		ImageDimension dim = new ImageDimension(image);
		dim.setZDim(ids.length);
		//Set information about columns
		visibleLabels = new ArrayList<ColumnInfo>();
		renderLabels = false;
		for (ColumnInfo ci: data.labels) {
			if (ci.visible)
				visibleLabels.add(ci);
			if (ci.render)
				renderLabels = true;
		}
		image.destroy();
		return dim;
	}
	
	@Override
	protected ImageItem createItem(int index, String key) throws Exception {
		return createImageItem(index, renderLabel, displayLabel, key);
	}
	
	/** Function to create an image item 
	 * @throws Exception */
	protected ImageItem createImageItem(int index, int renderLabel, int displayLabel, String key) throws Exception{
		ImageGeneric image = getImage(index, renderLabel);
		image.readApplyGeo(data.md, ids[index], thumb_width, thumb_height);
		ImagePlus imp = XmippImageConverter.convertImageGenericToImageJ(image);
		String labelStr = data.md.getValueString(displayLabel, ids[index]);
		return new ImageItem(key, labelStr, imp);		
	}

	@Override
	protected String getItemKey(int index) throws Exception {
		return getItemKey(index, renderLabel);
	}
	
	/** Return a key string using label 
	 * @throws Exception */
	protected String getItemKey(int index, int label) throws Exception{
		String filename = data.md.getValueString(label, ids[index]);
		return String.format("%s(%d,%d)", filename, thumb_width, thumb_height);
	}
	
	@Override
	public String getTitle() {
		return String.format("Metadata: %s (%d)", filename, n);
	}

	/**
	 * Create and returns the image generic from a given an index
	 * 
	 * @param index
	 *            Index of the image
	 * @return ImageGeneric created
	 * @throws Exception
	 *             if can not load image
	 */
	protected ImageGeneric getImage(int index, int label) throws Exception {
		String imgFn = data.md.getValueString(label, ids[index]);
		return new ImageGeneric(imgFn);
	}

	@Override
	protected double[] getMinAndMax() {
		try {
			return data.md.getStatistics(false);
		} catch (Exception ex) {
			DEBUG.printException(ex);
		}
		return null;
	}
}
