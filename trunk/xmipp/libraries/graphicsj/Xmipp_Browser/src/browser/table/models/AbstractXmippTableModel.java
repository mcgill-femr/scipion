/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.models;

import browser.Cache;
import browser.DEBUG;
import browser.ICONS_MANAGER;
import browser.imageitems.ImageConverter;
import browser.imageitems.tableitems.AbstractTableImageItem;
import ij.IJ;
import ij.ImagePlus;
import java.io.File;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import javax.swing.table.AbstractTableModel;
import xmipp.MDLabel;

/**
 *
 * @author Juanjo Vega
 */
public abstract class AbstractXmippTableModel extends AbstractTableModel {

    protected String filename;  // Filename (specified by user).
    protected String path;  // Normalized path.
    protected ArrayList<AbstractTableImageItem> data = new ArrayList<AbstractTableImageItem>();
    protected LinkedList<AbstractTableImageItem> selectedItems = new LinkedList<AbstractTableImageItem>();
    protected int rows, cols;
    protected Cache<String, ImagePlus> cache = new Cache<String, ImagePlus>();
    protected boolean normalize = false;
    protected double min = Double.POSITIVE_INFINITY, max = Double.NEGATIVE_INFINITY;
    protected boolean showLabels = false;
    protected int selectedLabel = MDLabel.MDL_IMAGE;
    protected int labelsValues[];
    protected boolean sort = false;
    protected int index[];

    public AbstractXmippTableModel() {
        super();
    }

    public AbstractXmippTableModel(String filename) {
        super();

        this.filename = filename;
        path = new File(new File(filename).toURI().normalize()).getAbsolutePath();

        String message = null;

        File f = new File(path);
        if (f.exists()) {
            message = populateTable(path);
        } else {
            message = "File not found: " + filename;
        }

        if (message != null) {
            IJ.error(message);
        }

//        printMD(md);
    }

    public void setShowLabels(boolean showLabels) {
        this.showLabels = showLabels;
    }

    public boolean isShowingLabels() {
        return showLabels;
    }

    public abstract String[] getLabels();

    public int getSelectedLabel() {
        return selectedLabel;
    }

    public void setSelectedLabel(int selectedLabel) {
        this.selectedLabel = labelsValues[selectedLabel];
    }

    protected abstract String populateTable(String filename);

//    private static void printMD(MetaData md) {
//        int labels[] = md.getActiveLabels();
//
//        for (int i : labels) {
//            String label = MetaData.label2Str(i);
//            System.out.println(i + ": " + label);
//        }
//
//        long ids[] = md.findObjects();
//        for (long id : ids) {
//            System.out.println(id + ": "
//                    + md.getValueString(MDLabel.MDL_IMAGE, id, true) + " / "
//                    + md.getValueInt(MDLabel.MDL_ENABLED, id));
//        }
//
//        System.out.println("++++++++++++++++++++++++++++++++++++++++");
//    }
    protected abstract void getMinAndMax();

    public abstract String getFilename();

    public ArrayList<AbstractTableImageItem> getAllItems() {
        return data;
    }

    public LinkedList<AbstractTableImageItem> getSelectedItems() {
        return selectedItems;
    }

    @Override
    public String getColumnName(int column) {
        return "" + (column + 1);
    }

    @Override
    public Class getColumnClass(int c) {
        Object object = getValueAt(0, c);

        return object != null ? object.getClass() : null;
    }

    public int getRowCount() {
        return rows;
    }

    public int getColumnCount() {
        return cols;
    }

    public boolean isSorting() {
        return sort;
    }

    public void setSorting(boolean sort) {
        this.sort = sort;

        DEBUG.printMessage(sort ? " *** Sorting" : " *** UN sorting");

        sort();
    }

    public void updateSort() {
        if (isSorting()) {
            DEBUG.printMessage(" *** Updating sort...");
            sort();
        }
    }

    private void sort() {
        // Creates index and...
        index = new int[getSize()];
        for (int i = 0; i < index.length; i++) {
            index[i] = i;
        }

        // ...sorts.
        mergeSort(data, index, selectedLabel);
    }

    /**
     * Mergesort algorithm.
     * @param data an array of Comparable items.
     * @param index index used to avoid moving data.
     */
    protected static void mergeSort(ArrayList data, int index[], int label) {
        int tmpArray[] = new int[data.size()];
        mergeSort(data, tmpArray, 0, data.size() - 1, index, label);
    }

    /**
     * Internal method that makes recursive calls.
     * @param data an array of Comparable items.
     * @param tmpArray an array to place the merged result.
     * @param left the left-most index of the subarray.
     * @param right the right-most index of the subarray.
     * @param index index used to avoid moving data.
     */
    private static void mergeSort(ArrayList<AbstractTableImageItem> data, int tmpArray[],
            int left, int right, int index[], int label) {
        if (left < right) {
            int center = (left + right) / 2;
            mergeSort(data, tmpArray, left, center, index, label);
            mergeSort(data, tmpArray, center + 1, right, index, label);
            merge(data, tmpArray, left, center + 1, right, index, label);
        }
    }

    /**
     * Internal method that merges two sorted halves of a subarray.
     * @param data an array of Comparable items.
     * @param tmpArray an array to place the merged result.
     * @param leftPos the left-most index of the subarray.
     * @param rightPos the index of the start of the second half.
     * @param rightEnd the right-most index of the subarray.
     * @param index index used to avoid moving data.
     */
    private static void merge(ArrayList<AbstractTableImageItem> data, int tmpArray[],
            int leftPos, int rightPos, int rightEnd, int index[], int label) {
        int leftEnd = rightPos - 1;
        int tmpPos = leftPos;
        int numElements = rightEnd - leftPos + 1;

        // Main loop
        while (leftPos <= leftEnd && rightPos <= rightEnd) {
            if (data.get(index[leftPos]).compareToByLabel(
                    data.get(index[rightPos]), label) <= 0) {
                tmpArray[tmpPos++] = index[leftPos++];//data[index[leftPos++]];
            } else {
                tmpArray[tmpPos++] = index[rightPos++];//data[index[rightPos++]];
            }
        }

        while (leftPos <= leftEnd) { // Copy rest of first half
            tmpArray[tmpPos++] = index[leftPos++];//data[index[leftPos++]];
        }

        while (rightPos <= rightEnd) { // Copy rest of right half
            tmpArray[tmpPos++] = index[rightPos++];//data[index[rightPos++]];
        }

        // Copy tmpArray back
        for (int i = 0; i < numElements; i++, rightEnd--) {
            //data[rightEnd] = tmpArray[rightEnd];
            index[rightEnd] = tmpArray[rightEnd];
        }
    }

    public Object getValueAt(int rowIndex, int columnIndex) {
        int i = getDataIndex(rowIndex, columnIndex);

        // If sorting, data will be accessed through the index.
        if (isSorting() && i < data.size()) {
            i = index[i];
        }

        return i < data.size() ? data.get(i) : null;
    }

    private int getDataIndex(int rowIndex, int columnIndex) {
        return rowIndex * getColumnCount() + columnIndex;
    }

    public int[] getRowColForIndex(int index) {
        int row = index / getColumnCount();
        int col = index % getColumnCount();

        return new int[]{row, col};
    }

    public int getSize() {
        return data.size();
    }

    public abstract String getTitle();

    public void setEnabledFrom(int row, int col, boolean enable) {
        int from = getDataIndex(row, col);

        ArrayList<AbstractTableImageItem> items = new ArrayList<AbstractTableImageItem>();

        for (int i = from; i < data.size(); i++) {
            items.add(data.get(i));
        }

        enableItems(items, enable);
    }

    public void setEnabledTo(int row, int col, boolean enable) {
        int to = getDataIndex(row, col);

        ArrayList<AbstractTableImageItem> items = new ArrayList<AbstractTableImageItem>();

        for (int i = 0; i < getSize() && i <= to; i++) {
            items.add(data.get(i));
        }

        enableItems(items, enable);
    }

    public void selectRange(int initial_row, int initial_col, int final_row, int final_col) {
        int row0 = Math.min(initial_row, final_row);
        int row1 = Math.max(initial_row, final_row);
        int col0 = Math.min(initial_col, final_col);
        int col1 = Math.max(initial_col, final_col);

        for (int i = row0; i <= row1; i++) {
            for (int j = col0; j <= col1; j++) {
                setSelected(i, j, true);
            }
        }
    }

    public void setSelected(int row, int col, boolean selected) {
        AbstractTableImageItem item = (AbstractTableImageItem) getValueAt(row, col);

        if (item != null) {
            item.setSelected(selected);

            if (selected) {
                if (!selectedItems.contains(item)) {
                    selectedItems.addLast(item);
                }
            } else {
                if (selectedItems.contains(item)) {
                    selectedItems.remove(item);
                }
            }
        }
        //fireTableCellUpdated(row, col);
    }

    public void toggleSelected(int row, int col) {
        AbstractTableImageItem item = (AbstractTableImageItem) getValueAt(row, col);

        if (item != null) {
            setSelected(row, col, !item.isSelected());
        }
    }

    // To select image directly
    public void setSelected(int index) {
        clearSelection();

        AbstractTableImageItem item = data.get(index);
        item.setSelected(true);
        selectedItems.addLast(item);
    }

    public void selectAll() {
        clearSelection();

        for (int i = 0; i < getSize(); i++) {
            AbstractTableImageItem item = data.get(i);

            item.setSelected(true);
            selectedItems.addLast(item);
        }
    }

    public void invertSelection() {
        for (int i = 0; i < getRowCount(); i++) {
            for (int j = 0; j < getColumnCount(); j++) {
                toggleSelected(i, j);
            }
        }
    }

    public void clearSelection() {
        int n = selectedItems.size();   // Size will change over iterations while index will be increased.

        for (int i = 0; i < n; i++) {
            selectedItems.getLast().setSelected(false);
            selectedItems.removeLast();
        }
    }

    private void enableItems(List<AbstractTableImageItem> items, boolean enable) {
        if (items != null) {
            for (int i = 0; i < items.size(); i++) {
                items.get(i).setEnabled(enable);
            }

            updateSort();   // Just in case it's sorting by field "enabled".
            fireTableDataChanged();
        }
    }

    public void enableAllItems() {
        enableItems(data, true);
    }

    public void disableAllItems() {
        enableItems(data, false);
    }

    public void enableSelectedItems() {
        enableItems(selectedItems, true);
    }

    public void disableSelectedItems() {
        enableItems(selectedItems, false);
    }

    public void setZoomScale(double zoomScale) {
        for (int i = 0; i < getSize(); i++) {
            data.get(i).setZoomScale(zoomScale);
        }
    }

    public void setRows(int rows) {
        if (rows > getSize()) {
            rows = getSize();
        }

        if (rows > 0) {
            this.rows = rows;

            // Calculates necessary cols for desired number of rows.
            cols = (int) Math.ceil((double) data.size() / (double) rows);

//            DEBUG.printMessage("R:" + rows + " / C:" + cols + " > S: " + (rows * cols));

            fireTableStructureChanged();
        }
    }

    public void setColumns(int cols) {
        if (cols > getSize()) {
            cols = getSize();
        }

        if (cols > 0) {
            this.cols = cols;

            // Calculates necessary rows for the desired number of cols.
            rows = (int) Math.ceil((double) data.size() / (double) cols);

            fireTableStructureChanged();
        }
    }

    public double getInitialZoomScale(int width, int intercellWidth) {
        double scale = 1.0;
        int i = 0;

        for (i = 0; i < getSize(); i++) {
            AbstractTableImageItem item = data.get(i);

            if (item.exists()) {
                scale = (double) width / (double) (item.getWidth() - 2 * intercellWidth);
                break;
            }
        }

        return scale > 1.0 ? 1.0 : scale;
    }

    public void autoAdjustColumns(int width, int intercellWidth) {
        int displayableColumns = (int) Math.floor(
                (double) width / (double) (getCellWidth() - 2 * intercellWidth));

//        double ddisplayableColumns =
//                (double) width / (double) (getCellWidth() - 2 * intercellWidth);
//
//        DEBUG.printMessage(" *** Displayable  columns: " + displayableColumns);
//        DEBUG.printMessage(" *** DDisplayable columns: " + ddisplayableColumns);

        if (getColumnCount() != displayableColumns) {
            setColumns(displayableColumns);
        }
    }

    public int getCellWidth() {
        int i = 0;

        for (i = 0; i < getSize(); i++) {
            if (data.get(i).exists()) {
                return getAllItems().get(i).getThumbnailWidth();
            }
        }

        return ICONS_MANAGER.DEFAULT_PREVIEW_WIDTH;
    }

    public int getCellHeight() {
        int i = 0;

        for (i = 0; i < getSize(); i++) {
            if (data.get(i).exists()) {
//                System.out.println("H: " + getAllItems().get(i).getThumbnailHeight());

                return getAllItems().get(i).getThumbnailHeight();
            }
        }

        return ICONS_MANAGER.DEFAULT_PREVIEW_HEIGHT;
    }

    public void setNormalized(boolean normalize) {
        this.normalize = normalize;

        if (normalize) {
            if (min == Double.POSITIVE_INFINITY && max == Double.NEGATIVE_INFINITY) {
                DEBUG.printMessage(" *** Retrieving minAndMax.");
                getMinAndMax();
            }
        }
    }

    public boolean isNormalizing() {
        return min > Double.NEGATIVE_INFINITY && max < Double.POSITIVE_INFINITY && normalize;
    }

    public double getNormalizeMin() {
        return min;
    }

    public double getNormalizeMax() {
        return max;
    }

    public abstract boolean isStack();

    public abstract boolean isVolume();

    public abstract boolean saveAsMetadata(String path);

    public boolean saveAsStack(String path) {
        try {
            ImagePlus imp = ImageConverter.convertToImagePlus(data);
            IJ.run(imp, "Xmipp writer", "save=" + path);

            return true;
        } catch (Exception ex) {
        }

        return false;
    }
}
