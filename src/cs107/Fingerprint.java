package cs107;

import java.util.List;
import java.util.ArrayList;

/**
 * Provides tools to compare fingerprint.
 */
public class Fingerprint {

    /**
     * The number of pixels to consider in each direction when doing the linear
     * regression to compute the orientation.
     */
    public static final int ORIENTATION_DISTANCE = 16;

    /**
     * The maximum distance between two minutiae to be considered matching.
     */
    public static final int DISTANCE_THRESHOLD = 5;

    /**
     * The number of matching minutiae needed for two fingerprints to be considered
     * identical.
     */
    public static final int FOUND_THRESHOLD = 20;

    /**
     * The distance between two angle to be considered identical.
     */
    public static final int ORIENTATION_THRESHOLD = 20;

    /**
     * The offset in each direction for the rotation to test when doing the
     * matching.
     */
    public static final int MATCH_ANGLE_OFFSET = 2;

    /**
     * Returns an array containing the value of the 8 neighbours of the pixel at
     * coordinates <code>(row, col)</code>.
     * <p>
     * The pixels are returned such that their indices corresponds to the following
     * diagram:<br>
     * ------------- <br>
     * | 7 | 0 | 1 | <br>
     * ------------- <br>
     * | 6 | _ | 2 | <br>
     * ------------- <br>
     * | 5 | 4 | 3 | <br>
     * ------------- <br>
     * <p>
     * If a neighbours is out of bounds of the image, it is considered white.
     * <p>
     * If the <code>row</code> or the <code>col</code> is out of bounds of the
     * image, the returned value should be <code>null</code>.
     *
     * @param image array containing each pixel's boolean value.
     * @param row   the row of the pixel of interest, must be between
     *              <code>0</code>(included) and
     *              <code>image.length</code>(excluded).
     * @param col   the column of the pixel of interest, must be between
     *              <code>0</code>(included) and
     *              <code>image[row].length</code>(excluded).
     * @return An array containing each neighbours' value.
     */

    public static boolean[] getNeighbours(boolean[][] image, int row, int col) {
        assert (image != null); // special case that is not expected (the image is supposed to have been earlier)

        // get dimensions of the rectangular array image
        int imageHeight = image.length;
        boolean[] rowOne = image[0];
        int imageWidth = rowOne.length;

        // neighbours[] is the array to be returned by getNeighbours()
        boolean neighbours[] = new boolean[8];


        // filling the array for the general case with 8 neighbours
        if (row != 0 && col != 0 && row != imageHeight - 1 && col != imageWidth - 1) {
            neighbours[0] = image[row - 1][col];
            neighbours[1] = image[row - 1][col + 1];
            neighbours[2] = image[row][col + 1];
            neighbours[3] = image[row + 1][col + 1];
            neighbours[4] = image[row + 1][col];
            neighbours[5] = image[row + 1][col - 1];
            neighbours[6] = image[row][col - 1];
            neighbours[7] = image[row - 1][col - 1];
        }
        // filling the array but checking for the special bordering cases beforehand
        else {
            neighbours[0] = (row == 0) ? false : image[row - 1][col];
            neighbours[1] = (row == 0 || col == imageWidth - 1) ? false : image[row - 1][col + 1];
            neighbours[2] = (col == imageWidth - 1) ? false : image[row][col + 1];
            neighbours[3] = (row == imageHeight - 1 || col == imageWidth - 1) ? false : image[row + 1][col + 1];
            neighbours[4] = (row == imageHeight - 1) ? false : image[row + 1][col];
            neighbours[5] = (row == imageHeight - 1 || col == 0) ? false : image[row + 1][col - 1];
            neighbours[6] = (col == 0) ? false : image[row][col - 1];
            neighbours[7] = (row == 0 || col == 0) ? false : image[row - 1][col - 1];
        }

        return neighbours;
    }


    /**
     * Computes the number of black (<code>true</code>) pixels among the neighbours
     * of a pixel.
     *
     * @param neighbours array containing each pixel value. The array must respect
     *                   the convention described in
     *                   {@link #getNeighbours(boolean[][], int, int)}.
     * @return the number of black neighbours.
     */
    public static int blackNeighbours(boolean[] neighbours) {

        int blackNeighbours = 0;
        for (int i = 0; i < neighbours.length; i++) {
            if (neighbours[i] == true) {
                blackNeighbours++;
            }
        }
        return blackNeighbours;
    }

    /**
     * Computes the number of white to black transitions among the neighbours of
     * pixel.
     *
     * @param neighbours array containing each pixel value. The array must respect
     *                   the convention described in
     *                   {@link #getNeighbours(boolean[][], int, int)}.
     * @return the number of white to black transitions.
     */
    public static int transitions(boolean[] neighbours) {

        int transitions = 0;
        for (int i = 0; i <= 6; i++) {                                                          // only check until P_6
            if (neighbours[i] == false && (neighbours[i + 1] == true)) {
                transitions++;
            }
        }
        if (neighbours[neighbours.length - 1] == false && neighbours[0] == true) {           // check special last case P_7
            transitions++;
        }
        return transitions;
    }

    /**
     * Returns <code>true</code> if the images are identical and false otherwise.
     *
     * @param image1 array containing each pixel's boolean value.
     * @param image2 array containing each pixel's boolean value.
     * @return <code>True</code> if they are identical, <code>false</code>
     * otherwise.
     */


    public static boolean identical(boolean[][] image1, boolean[][] image2) {
        // iterates through the 2D array of image1 comparing every element of image1to its corresponding element in image2
        // they have the same dimensions
        for (int i = 0; i < image1.length; i++) {
            for (int j = 0; j < image1[i].length; j++) {
                boolean b1 = image1[i][j];
                boolean b2 = image2[i][j];
                if (b1 != b2) {
                    return false;
                }   // hop out of loop immediately if a discrepancy is detected
            }
        }
        return true;
    }

    public static boolean[][] copyArray(boolean[][] image) {
        boolean copy[][] = new boolean[image.length][image[0].length];
        for (int i = 0; i < image.length; i++) {
            for (int j = 0; j < image[i].length; j++) {
                copy[i][j] = image[i][j];
            }
        }
        return copy;
    }

    public static boolean step(boolean[][] image, int i, int j, int step) {

        //checking conditions for unnecessary pixels (i,j) that will be set to false if they meet all the conditions
        if (image[i][j] &&
                (2 <= blackNeighbours(getNeighbours(image, i, j)) && blackNeighbours(getNeighbours(image, i, j)) <= 6) &&
                transitions(getNeighbours(image, i, j)) == 1 &&
                (!getNeighbours(image, i, j)[0] || !getNeighbours(image, i, j)[2] || !getNeighbours(image, i, j)[2 * step + 4]) &&
                (!getNeighbours(image, i, j)[-2 * step + 2] || !getNeighbours(image, i, j)[4] || !getNeighbours(image, i, j)[6])) {
            return false;
        } else return image[i][j];
    }

    /**
     * Internal method used by {@link #thin(boolean[][])}.
     *
     * @param image array containing each pixel's boolean value.
     * @param step  the step to apply, Step 0 or Step 1.
     * @return A new array containing each pixel's value after the step.
     */


    public static boolean[][] thinningStep(boolean[][] image, int step) {
        boolean[][] thinningStep = new boolean[image.length][image[0].length];    // create new array with the dimensions of input image
        for (int i = 0; i < image.length; i++) {                                      // fill new array with evaluated (necessary) elements of input image
            for (int j = 0; j < image[i].length; j++) {
                thinningStep[i][j] = step(image, i, j, step);
            }
        }
        return thinningStep;
    }

    /**
     * Compute the skeleton of a boolean image.
     *
     * @param image array containing each pixel's boolean value.
     * @return array containing the boolean value of each pixel of the image after
     * applying the thinning algorithm.
     */
    public static boolean[][] thin(boolean[][] image) {
        boolean[][] initial = copyArray(image);

        boolean[][] second = thinningStep(initial, 0);      // first step applied on initial
        boolean[][] third = thinningStep(second, 1);        // second step applied on modified initial

        while (!identical(second, third)) {
            initial = copyArray(third);
            second = thinningStep(initial, 0);
            third = thinningStep(second, 1);
        }
        return third;
    }

    /**
     * Computes all pixels that are connected to the pixel at coordinate
     * <code>(row, col)</code> and within the given distance of the pixel.
     *
     * @param image    array containing each pixel's boolean value.
     * @param row      the first coordinate of the pixel of interest.
     * @param col      the second coordinate of the pixel of interest.
     * @param distance the maximum distance at which a pixel is considered.
     * @return An array where <code>true</code> means that the pixel is within
     * <code>distance</code> and connected to the pixel at
     * <code>(row, col)</code>.
     */
    public static boolean[][] connectedPixels(boolean[][] image, int row, int col, int distance) {
        boolean[][] connectedPixels = new boolean[image.length][image[0].length]; //Filled by default of false
        connectedPixels[row][col] = true; //light up the minutia


        boolean stillModifying = true;
        boolean isBlack;
        boolean isConnected;
        //**sub states of isConnected**
        boolean skip_top, top_left, top, top_right, right, bottom_right, bottom, bottom_left, left;
        skip_top = top_left = top = top_right = right = bottom_right = bottom = bottom_left = left = false;

        //*********
        boolean isInDistance;
        while (stillModifying) {
            stillModifying = false;
            for (int i = 0; i < image.length; i++) {
                for (int j = 0; j < image[0].length; j++) {
                    //testing if the cell at i,j is connected to the minutia

                    //Condition: the pixel at this position is black
                    isBlack = image[i][j]; //same as image[i][j] == true

                    // Condition: the cell has at least a black neighbour
                    int nb_of_black_neighbors = blackNeighbours(getNeighbours(connectedPixels, i, j));
                    isConnected = nb_of_black_neighbors > 0;

                    //Condition: The pixel is in the square of size 2*distance + 1 centered on the minutia
                    isInDistance = i <= row + distance && i >= row - distance && j <= col + distance && j >= col - distance;

                    if (isBlack && isConnected && isInDistance) {
                        if (!connectedPixels[i][j]) {
                            connectedPixels[i][j] = true;
                            stillModifying = true; // is used to spot when to terminate
                        }
                    }

                }
            }
        }
        return connectedPixels;
    }

    /**
     * Computes the slope of a minutia using linear regression.
     *
     * @param connectedPixels the result of
     *                        {@link #connectedPixels(boolean[][], int, int, int)}.
     * @param row             the row of the minutia.
     * @param col             the col of the minutia.
     * @return the slope.
     */
    public static double computeSlope(boolean[][] connectedPixels, int row, int col) {
        // compute somme of x*y
        double sumXTimesY = 0;

        // compute somme of x^2
        double sumXSquared = 0;
        // compute somme of y^2
        double sumYSquared = 0;
        //slope
        double slope;

        //useful placeholders
        double x, y;

        // compute sumXTimesY, sumXSquared, sumYSquared
        for (int i = 0; i < connectedPixels.length; ++i) {
            for (int j = 0; j < connectedPixels[0].length; ++j) {
                if (connectedPixels[i][j] && (i != row || j != col)) { // if black and is not the minutia
                    // cf. formulas p.17
                    x = j - col;
                    y = row - i;

                    sumXTimesY += x * y;
                    sumXSquared += Math.pow(x, 2);
                    sumYSquared += Math.pow(y, 2);
                }
            }
        }

        //takes care of the case where the slope is infinite(the line is vertical)
        if (sumXSquared == 0) {
            return Double.POSITIVE_INFINITY;
        }
        if (sumXSquared >= sumYSquared) {
            slope = sumXTimesY / sumXSquared;
        } else {
            slope = sumYSquared / sumXTimesY;
        }
        return slope;
    }

    /**
     * Computes the orientation of a minutia in radians.
     *
     * @param connectedPixels the result of
     *                        {@link #connectedPixels(boolean[][], int, int, int)}.
     * @param row             the row of the minutia.
     * @param col             the col of the minutia.
     * @param slope           the slope as returned by
     *                        {@link #computeSlope(boolean[][], int, int)}.
     * @return the orientation of the minutia in radians.
     */
    public static double computeAngle(boolean[][] connectedPixels, int row, int col, double slope) {
        double defaultAngle;
        // treats the case where slope was infinity
        if (slope == Double.POSITIVE_INFINITY) {
            defaultAngle = Math.PI / 2;
        } else {
            defaultAngle = Math.atan(slope);
        }

        double correctAngle;
        int nbPixelsAbove = 0;
        int nbPixelsBelow = 0;
        // useful variables
        double x, y;

        // compute nbPixelsAbove and nbPixelsBelow
        for (int i = 0; i < connectedPixels.length; ++i) {
            for (int j = 0; j < connectedPixels[0].length; ++j) {
                if (connectedPixels[i][j] && (i != row || j != col)) { // if black and is not the minutia
                    // cf. formulas p.17
                    x = j - col;
                    y = row - i;

                    if (y >= (-1 / slope) * x) {
                        ++nbPixelsAbove;
                    } else {
                        ++nbPixelsBelow;
                    }
                }
            }
        }

        // gives the defaultAngle the proper direction -> correctAngle.
        // I could put the if and the else if together but less readable.
        if (defaultAngle > 0 && nbPixelsBelow > nbPixelsAbove) { // test if it should be in quadrant 3
            correctAngle = defaultAngle + Math.PI;
        } else if (defaultAngle < 0 && nbPixelsAbove > nbPixelsBelow) { // test if it should be in quadrant 2
            correctAngle = defaultAngle + Math.PI;
        } else { // else it means that it is correctly in quadrant 1 or 4.
            correctAngle = defaultAngle;
        }
        return correctAngle;
    }

    /**
     * Computes the orientation of the minutia that the coordinate <code>(row,
     * col)</code>.
     *
     * @param image    array containing each pixel's boolean value.
     * @param row      the first coordinate of the pixel of interest.
     * @param col      the second coordinate of the pixel of interest.
     * @param distance the distance to be considered in each direction to compute
     *                 the orientation.
     * @return The orientation in degrees.
     */
    public static int computeOrientation(boolean[][] image, int row, int col, int distance) {

        boolean[][] connectedPixels = connectedPixels(image, row, col, distance);
        double slope = computeSlope(connectedPixels, row, col);
        double angle = computeAngle(connectedPixels, row, col, slope);
        int angleDegree = (int) Math.round(Math.toDegrees(angle));

        if (angleDegree < 0) { // handles negative angles
            angleDegree += 360;
        }

        return angleDegree;
    }

    /**
     * Extracts the minutiae from a thinned image.
     *
     * @param image array containing each pixel's boolean value.
     * @return The list of all minutiae. A minutia is represented by an array where
     * the first element is the row, the second is column, and the third is
     * the angle in degrees.
     * @see #thin(boolean[][])
     */
    public static List<int[]> extract(boolean[][] image) {
        //** useful variables **
        int orientation;
        List<int[]> minutiae = new ArrayList<int[]>(); // the int[] characterising the minutia is [row, col, orientation(degrees)]
        //**********************

        // Go through the whole image, get the no. of transitions, check if it is a minutia, if so add it to minutiae array with the direction data.
        for (int i = 1; i < image.length - 1; ++i) { // we don't loop through the edges because we want minutiae with 8 neigh.
            for (int j = 1; j < image[0].length - 1; ++j) {
                int nbTransition = transitions(getNeighbours(image, i, j));
                if (image[i][j] && (nbTransition == 1 || nbTransition == 3)) { // this tests if it is a minutia
                    orientation = computeOrientation(image, i, j, ORIENTATION_DISTANCE);
                    minutiae.add(new int[]{i, j, orientation});  // [row, col, orientation(degrees)]
                }
            }

        }

        return minutiae;
    }

    /**
     * Applies the specified rotation to the minutia.
     *
     * @param minutia   the original minutia.
     * @param centerRow the row of the center of rotation.
     * @param centerCol the col of the center of rotation.
     * @param rotation  the rotation in degrees.
     * @return the minutia rotated around the given center.
     */
    public static int[] applyRotation(int[] minutia, int centerRow, int centerCol, int rotation) {

        int x = minutia[1] - centerCol;
        int y = centerRow - minutia[0];
        double rotationRadians = Math.toRadians(rotation);                                   // rotation from degrees to radians
        double newX_notRounded = x * Math.cos(rotationRadians) - y * Math.sin(rotationRadians);
        double newY_notRounded = x * Math.sin(rotationRadians) + y * Math.cos(rotationRadians);
        int newX = (int) Math.round(newX_notRounded);                                        // cast from double to int
        int newY = (int) Math.round(newY_notRounded);                                        // cast from double to int
        int newRow = centerRow - newY;
        int newCol = newX + centerCol;
        double newOrientation_notRounded = (minutia[2] + rotation) % 360;
        int newOrientation = (int) newOrientation_notRounded;
        int[] rotatedMinutia = {newRow, newCol, newOrientation};                              // fill rotated results into new array
        return rotatedMinutia;
    }

    /**
     * Applies the specified translation to the minutia.
     *
     * @param minutia        the original minutia.
     * @param rowTranslation the translation along the rows.
     * @param colTranslation the translation along the columns.
     * @return the translated minutia.
     */
    public static int[] applyTranslation(int[] minutia, int rowTranslation, int colTranslation) {
        int newRow = minutia[0] - rowTranslation;
        int newCol = minutia[1] - colTranslation;
        int newOrientation = minutia[2];
        int[] translatedMinutia = {newRow, newCol, newOrientation};    // fill rotated results into new array
        return translatedMinutia;
    }

    /**
     * Computes the row, column, and angle after applying a transformation
     * (translation and rotation).
     *
     * @param minutia        the original minutia.
     * @param centerCol      the column around which the point is rotated.
     * @param centerRow      the row around which the point is rotated.
     * @param rowTranslation the vertical translation.
     * @param colTranslation the horizontal translation.
     * @param rotation       the rotation.
     * @return the transformed minutia.
     */
    public static int[] applyTransformation(int[] minutia, int centerRow, int centerCol, int rowTranslation,
                                            int colTranslation, int rotation) {
        int[] rotatedMinutia = applyRotation(minutia, centerRow, centerCol, rotation);
        int[] rotated_translatedMinutia = applyTranslation(rotatedMinutia, rowTranslation, colTranslation);
        return rotated_translatedMinutia;
    }

    /**
     * Computes the row, column, and angle after applying a transformation
     * (translation and rotation) for each minutia in the given list.
     *
     * @param minutiae       the list of minutiae.
     * @param centerCol      the column around which the point is rotated.
     * @param centerRow      the row around which the point is rotated.
     * @param rowTranslation the vertical translation.
     * @param colTranslation the horizontal translation.
     * @param rotation       the rotation.
     * @return the list of transformed minutiae.
     */
    public static List<int[]> applyTransformation(List<int[]> minutiae, int centerRow, int centerCol, int rowTranslation,
                                                  int colTranslation, int rotation) {

        ArrayList<int[]> transformedMinutiae = new ArrayList<int[]>();

        for (int i = 0; i < minutiae.size(); i++) {
            transformedMinutiae.add(i, applyTransformation(minutiae.get(i), centerRow, centerCol, rowTranslation, colTranslation, rotation));
        }
        return transformedMinutiae;
    }

    public static boolean checkSuperImposed(int[] list1Minutia, int[] list2Minutia, int maxDistance, int maxOrientation) {
        double euclideanDistance = Math.sqrt((list1Minutia[0] - list2Minutia[0]) * (list1Minutia[0] - list2Minutia[0])
                + (list1Minutia[1] - list2Minutia[1]) * (list1Minutia[1] - list2Minutia[1]));

        int differenceOrientation = Math.abs(list1Minutia[2] - list2Minutia[2]);

        if (euclideanDistance <= maxDistance && differenceOrientation <= maxOrientation) {
            return true;                                                                // the two compared minutia match
        }

        return false;
    }

    /**
     * Counts the number of overlapping minutiae.
     *
     * @param minutiae1      the first set of minutiae.
     * @param minutiae2      the second set of minutiae.
     * @param maxDistance    the maximum distance between two minutiae to consider
     *                       them as overlapping.
     * @param maxOrientation the maximum difference of orientation between two
     *                       minutiae to consider them as overlapping.
     * @return the number of overlapping minutiae.
     */


    public static int matchingMinutiaeCount(List<int[]> minutiae1, List<int[]> minutiae2, int maxDistance, int maxOrientation) {
        int count = 0;

        for (int i = 0; i < minutiae1.size(); i++) {
            for (int j = 0; j < minutiae2.size(); j++) {
                if (checkSuperImposed(minutiae1.get(i), minutiae2.get(j), maxDistance, maxOrientation)) {
                    count++;
                }
            }
        }

        return count;
    }

    /**
     * Compares the minutiae from two fingerprints.
     *
     * @param minutiae1 the list of minutiae of the first fingerprint.
     * @param minutiae2 the list of minutiae of the second fingerprint.
     * @return Returns <code>true</code> if they match and <code>false</code>
     * otherwise.
     */
    public static boolean match(List<int[]> minutiae1, List<int[]> minutiae2) {

        for (int i = 0; i < minutiae1.size(); i++) {
            int[] m_1 = minutiae1.get(i);
            for (int j = 0; j < minutiae2.size(); j++) {
                int[] m_2 = minutiae2.get(j);
                int rotation = m_2[2] - m_1[2];
                for (int k = rotation - MATCH_ANGLE_OFFSET; k < rotation + MATCH_ANGLE_OFFSET; k++) {
                    List<int[]> transformedMinutiae =
                            applyTransformation(minutiae2, minutiae1.get(i)[0], minutiae1.get(i)[1],
                                    minutiae2.get(j)[0] - minutiae1.get(i)[0], minutiae2.get(j)[1] - minutiae1.get(i)[1], rotation);

                    int resultCount = matchingMinutiaeCount(minutiae1, transformedMinutiae, DISTANCE_THRESHOLD, ORIENTATION_THRESHOLD);
                    if (resultCount >= FOUND_THRESHOLD) return true;

                }

            }

        }

        return false;

    }
}
