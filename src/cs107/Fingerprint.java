package cs107;

import jdk.swing.interop.DropTargetContextWrapper;

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
   *
   */

    public static boolean[] getNeighbours(boolean[][] image, int row, int col) {
	  assert (image != null); // special case that is not expected (the image is supposed to have been checked
                              // earlier)

      // get dimensions of the rectangular array image
      int imageHeight = image.length;
      boolean [] rowOne = image[0];
      int imageWidth = rowOne.length;

      // neighbours[] is the array to be returned by getNeighbours()
      boolean neighbours [] = new boolean[8]; 

      // filling the array and checking for special cases beforehand
      neighbours[0] = (row==0)? false : image[row-1][col];
      neighbours[1] = (row==0 || col==imageWidth-1)? false : image[row-1][col+1];
      neighbours[2] = (col==imageWidth-1)? false : image[row][col+1];
      neighbours[3] = (row==imageHeight-1 || col==imageWidth-1)?false : image[row+1][col+1];
      neighbours[4] = (row==imageHeight-1)? false : image[row+1][col];
      neighbours[5] = (row==imageHeight-1 || col==0)? false : image[row+1][col-1];
      neighbours[6] = (col==0)? false : image[row][col-1];
      neighbours[7] = (row==0 || col==0)?false : image[row-1][col-1];

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
    for (int i=0; i<neighbours.length; i++) {
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
      for (int a=0; a<=6; a++) {                                                          // only check until P_6
          if (neighbours[a] == false &&  (neighbours[a + 1] == true)) {
              transitions++;
          }
      }
      if (neighbours[neighbours.length-1] == false && neighbours[0] == true)  {           // check special last case P_7
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
   *         otherwise.
   */


  public static boolean identical(boolean[][] image1, boolean[][] image2) {
      // iterates through the 2D array of image1 comparing every element of image1to its corresponding element in image2
      // assume they have the same dimensions because elements can only flip their truth value
	  for (int i=0; i<image1.length; i++) {
          for (int j=0; j<image1[i].length; j++) {
              boolean b1 = image1[i][j];
              boolean b2 = image2[i][j];
              if (b1!=b2) return false;
          }
      }
	  return true;
  }

  /**
   * Internal method used by {@link #thin(boolean[][])}.
   *
   * @param image array containing each pixel's boolean value.
   * @param step  the step to apply, Step 0 or Step 1.
   * @return A new array containing each pixel's value after the step.
   */
  public static boolean[][] thinningStep(boolean[][] image, int step) {
	  //TODO implement
	  return null;
  }
  
  /**
   * Compute the skeleton of a boolean image.
   *
   * @param image array containing each pixel's boolean value.
   * @return array containing the boolean value of each pixel of the image after
   *         applying the thinning algorithm.
   */
  public static boolean[][] thin(boolean[][] image) {
	  //TODO implement
	  return null;
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
   *         <code>distance</code> and connected to the pixel at
   *         <code>(row, col)</code>.
   */
  public static boolean[][] connectedPixels(boolean[][] image, int row, int col, int distance) {
      boolean[][] connectedPixels = new boolean[image.length][image[0].length]; //Filled by default of falses
      connectedPixels[row][col] = true; //light up the minutia


      boolean still_modifying = true;
      boolean is_black;
      boolean is_connected;
      //**sub states of is_connected**
      boolean skip_top, top_left, top, top_right, right, bottom_right, bottom, bottom_left, left;
      skip_top = top_left = top = top_right = right = bottom_right = bottom = bottom_left = left = false;

      //*********
      boolean is_in_distance;
      while(still_modifying){
          still_modifying = false;
          for(int i = 0; i< image.length; i++) {
              for(int j = 0; j < image[0].length; j++) {
                //testing if the cell at i,j is connected to the minutia
                //Condition: the pixel at this position is black
                is_black = image[i][j]; //same as image[i][j] == true

                //Condition: the pixel is connected to our minutia therefore connected to a black pixel
                //All the try/catch are to not have any problems with pixels that are on the edge.
                // Didn't find a simpler way to handle it.All others wher more complex.
                top_left = top = top_right = right = bottom_right = bottom = bottom_left = left = false;
                try{
                    top_right = connectedPixels[i-1][j+1];
                }catch(Exception err){
                    ; // instead of "top_right = false;" because top_right is already false by default.
                }
                try{
                    right = connectedPixels[i][j+1];
                }catch(Exception err){
                    ;
                }
                try{
                    bottom_right = connectedPixels[i+1][j+1];
                }catch(Exception err){
                    ;
                }
                try{
                    bottom = connectedPixels[i+1][j];
                }catch(Exception err){
                    ;
                }
                try{
                    bottom_left = connectedPixels[i+1][j-1];
                }catch(Exception err){
                    ;
                }
                try{
                    left = connectedPixels[i][j-1];
                }catch(Exception err){
                    ;
                }
                try{
                    top_left = connectedPixels[i-1][j-1];
                }catch(Exception err){
                    ;
                }
                try{
                    top = connectedPixels[i-1][j];
                }catch(Exception err){
                    ;
                }

                is_connected = right || left || bottom || top || top_left || top_right || bottom_right || bottom_left;

                //Condition: The pixel is in the square of size 2*distance + 1 centered on the minutia
                is_in_distance = i <= row + distance && i >= row - distance && j <= col + distance && j >= col - distance;

                if(is_black && is_connected && is_in_distance){
                  if(!connectedPixels[i][j]) {
                    connectedPixels[i][j] = true;
                    still_modifying = true; // is used to spot when to terminate
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
      double sum_x_times_y = 0;

      // compute somme of x^2
      double sum_x_squared = 0;
      // compute somme of y^2
      double sum_y_squared = 0;
      //slope
      double slope;

      //useful placeholders
      double x, y;

      // compute sum_x_times_y, sum_x_squared, sum_y_squared
      for(int i = 0; i < connectedPixels.length; ++i) {
          for (int j = 0; j < connectedPixels[0].length; ++j) {
            if(connectedPixels[i][j] && (i != row || j != col)){ // if black and is not the minutia
                // cf. formulas p.17
                x = j-col;
                y = row-i;

                sum_x_times_y += x*y;
                sum_x_squared += Math.pow(x,2);
                sum_y_squared += Math.pow(y,2);
            }
          }
      }

      //takes care of the case where the slope is infinite(the line is vertical)
      if(sum_x_squared == 0){
          return Double.POSITIVE_INFINITY;
      }
      if(sum_x_squared >= sum_y_squared){
          slope = sum_x_times_y/sum_x_squared;
      } else{
          slope = sum_y_squared/sum_x_times_y;
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
      double default_angle;
      // treats the case where slope was infinity
      if(slope == Double.POSITIVE_INFINITY){
          default_angle = Math.PI/2;
      }else{
          default_angle = Math.atan(slope);
      }

      double correct_angle;
      int nb_pixels_above = 0;
      int nb_pixels_below = 0;
      // useful variables
      double x, y;

      // compute nb_pixels_above and nb_pixels_below
      for(int i = 0; i < connectedPixels.length; ++i) {
          for (int j = 0; j < connectedPixels[0].length; ++j) {
              if(connectedPixels[i][j] && (i != row || j != col)){ // if black and is not the minutia
                  // cf. formulas p.17
                  x = j-col;
                  y = row-i;

                  if(y >= (-1/slope)*x){
                      ++nb_pixels_above;
                  }else{
                      ++nb_pixels_below;
                  }
              }
          }
      }

      // gives the default_angle the proper direction -> correct_angle.
      // I could put the if and the else if together but less readable.
      if(default_angle > 0 && nb_pixels_below > nb_pixels_above){ // test if it should be in quadrant 3
          correct_angle = default_angle + Math.PI;
      } else if(default_angle < 0 && nb_pixels_above > nb_pixels_below){ // test if it should be in quadrant 2
          correct_angle = default_angle + Math.PI;
      } else{ // else it means that it is correctly in quadrant 1 or 4.
          correct_angle = default_angle;
      }
	  return correct_angle;
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

      boolean [][] connected_pixels = connectedPixels(image, row, col, distance);
      double slope = computeSlope(connected_pixels, row, col);
      double angle = computeAngle(connected_pixels, row, col, slope);
      int angle_degree = (int)Math.round(Math.toDegrees(angle));

      if(angle_degree < 0){ // handles negative angles
          angle_degree += 360;
      }

	  return angle_degree;
  }

  /**
   * Extracts the minutiae from a thinned image.
   *
   * @param image array containing each pixel's boolean value.
   * @return The list of all minutiae. A minutia is represented by an array where
   *         the first element is the row, the second is column, and the third is
   *         the angle in degrees.
   * @see #thin(boolean[][])
   */
  public static List<int[]> extract(boolean[][] image) {
	  //TODO test it once transition() and getNeighbours() is done.
      //** useful variables **
      int orientation;
      List<int[]> minutiae = new ArrayList<int[]>(); // the int[] characterising the minutia is [row, col, orientation(degrees)]
      //**********************

      // Go through the whole image, get the nb of transitions, check if it is a minutia, if so add it to minutiae array with the direction data.
      for(int i = 1; i < image.length-1; ++i) { // we don't loop through the edges because we want minutiae with 8 neigh.
          for (int j = 1; j < image[0].length-1; ++j) {
              int nb_transition = transitions(getNeighbours(image, i, j));
              if(image[i][j] && (nb_transition == 1 || nb_transition == 3)){ // this tests if it is a minutia
                  orientation = computeOrientation(image, i, j, ORIENTATION_DISTANCE);
                  minutiae.add(new int[] {i, j, orientation});  // [row, col, orientation(degrees)]
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
	  //TODO implement
	  return null;
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
	  //TODO implement
	  return null;
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
	  //TODO implement
	  return null;
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
	  //TODO implement
	  return null;
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
  public static int matchingMinutiaeCount(List<int[]> minutiae1, List<int[]> minutiae2, int maxDistance,
      int maxOrientation) {
	  //TODO implement
	  return 0;
  }

  /**
   * Compares the minutiae from two fingerprints.
   *
   * @param minutiae1 the list of minutiae of the first fingerprint.
   * @param minutiae2 the list of minutiae of the second fingerprint.
   * @return Returns <code>true</code> if they match and <code>false</code>
   *         otherwise.
   */
  public static boolean match(List<int[]> minutiae1, List<int[]> minutiae2) {
	  //TODO implement
	  return false;
  }
}
