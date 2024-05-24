import com.opencsv.CSVReader;
import com.opencsv.exceptions.CsvException;
import java.awt.*;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import org.knowm.xchart.*;
import org.knowm.xchart.style.lines.SeriesLines;
import org.knowm.xchart.style.markers.SeriesMarkers;

public class kMeans {

  public static void main(String[] args) {
    double[][] irisData = new double[0][0];
    int numClusters1 = 2;
    int numClusters2 = 3;
    String csvFile = "irisdata.csv";

    try (CSVReader reader = new CSVReader(new FileReader(csvFile))) {
      List<String[]> data = reader.readAll();

      String[] labels = data.remove(0);

      irisData = new double[data.size()][data.get(0).length - 1];

      for (int i = 0; i < data.size(); i++) {
        String[] row = data.get(i);
        for (int j = 0; j < row.length - 1; j++) {
          irisData[i][j] = Double.parseDouble(row[j]);
        }
      }
    } catch (IOException | CsvException | NumberFormatException e) {
      e.printStackTrace();
    }
    // Plotting for Objective Function
    XYChart objectiveChart = new XYChartBuilder()
      .width(800)
      .height(600)
      .title("Objective Function Value Q1b")
      .xAxisTitle("Iteration")
      .yAxisTitle("Objective Function Value")
      .build();
    objectiveChart.getStyler().setMarkerSize(6);

    plotLearningCurve(
      irisData,
      numClusters1,
      objectiveChart,
      new Color(255, 0, 0)
    );
    plotLearningCurve(
      irisData,
      numClusters2,
      objectiveChart,
      new Color(0, 0, 255)
    );

    new SwingWrapper<>(objectiveChart).displayChart();

    // Plotting for Cluster Centers Overlaid on Data
    XYChart clusterCentersChart1 = new XYChartBuilder()
      .width(800)
      .height(600)
      .title("Cluster Centers Overlaid on Data - k=2  Q1c")
      .xAxisTitle("Petal Length")
      .yAxisTitle("Petal Width")
      .build();
    clusterCentersChart1
      .getStyler()
      .setDefaultSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Scatter);
    XYChart clusterCentersChart2 = new XYChartBuilder()
      .width(800)
      .height(600)
      .title("Cluster Centers Overlaid on Data - k=3  Q1c")
      .xAxisTitle("Petal Length")
      .yAxisTitle("Petal Width")
      .build();
    clusterCentersChart2
      .getStyler()
      .setDefaultSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Scatter);

    plotClustersWithData(
      irisData,
      numClusters1,
      clusterCentersChart1,
      new Color(255, 0, 0),
      "k=2"
    ); //k=2
    plotClustersWithData(
      irisData,
      numClusters2,
      clusterCentersChart2,
      new Color(0, 0, 255),
      "k=3"
    ); //k=3

    new SwingWrapper<>(clusterCentersChart1).displayChart();
    new SwingWrapper<>(clusterCentersChart2).displayChart();
    //plotting decision boundaries
    plotDecisionBoundaries(irisData, 2);
    plotDecisionBoundaries(irisData, 3);
  }

  public static void printMatrix(double[][] matrix) {
    int rows = matrix.length;
    int cols = matrix[0].length;

    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        System.out.print(matrix[i][j] + "\t"); // Use "\t" for tab spacing
      }
      System.out.println(); // Move to the next line for the next row
    }
  }

  private static double[][] initializeCentroids(double[][] data, int k) {
    double[][] centroids = new double[k][data[0].length];
    Random random = new Random();

    for (int i = 0; i < k; i++) {
      int randomIndex = random.nextInt(data.length);
      centroids[i] = data[randomIndex].clone();
    }
    return centroids;
  }

  private static int[] assignToClusters(double[][] data, double[][] centroids) {
    int[] clusters = new int[data.length];

    for (int i = 0; i < data.length; i++) {
      double minDistance = Double.MAX_VALUE;
      int assignedCluster = -1;

      for (int j = 0; j < centroids.length; j++) {
        double distance = calculateDistanceSquared(data[i], centroids[j]);
        if (distance < minDistance) {
          minDistance = distance;
          assignedCluster = j;
        }
      }

      clusters[i] = assignedCluster;
    }

    return clusters;
  }

  private static List<double[]> getClusterPoints(
    double[][] data,
    int[] clusters,
    int clusterIndex
  ) {
    List<double[]> clusterPoints = new ArrayList<>();

    for (int i = 0; i < data.length; i++) {
      if (clusters[i] == clusterIndex) {
        clusterPoints.add(data[i].clone());
      }
    }

    return clusterPoints;
  }

  private static double[] calculateNewCentroid(List<double[]> clusterPoints) {
    int dimensions = clusterPoints.get(0).length;
    double[] newCentroid = new double[dimensions];

    for (int d = 0; d < dimensions; d++) {
      double sum = 0;
      for (double[] point : clusterPoints) {
        sum += point[d];
      }
      newCentroid[d] = sum / clusterPoints.size();
    }

    return newCentroid;
  }

  private static double calculateObjectiveFunction(
    double[][] data,
    double[][] centroids,
    int[] clusters
  ) {
    double objective = 0;

    for (int n = 0; n < data.length; n++) {
      for (int k = 0; k < centroids.length; k++) {
        double distanceSquared = calculateDistanceSquared(
          data[n],
          centroids[k]
        );
        objective += clusters[n] == k ? distanceSquared : 0;
      }
    }

    return objective;
  }

  private static double calculateDistanceSquared(
    double[] point1,
    double[] point2
  ) {
    double sumSquared = 0;

    for (int i = 0; i < point1.length; i++) {
      double diff = point1[i] - point2[i];
      sumSquared += diff * diff;
    }

    return sumSquared;
  }

  private static void printArray(int[] array) {
    for (int value : array) {
      System.out.print(value + " ");
    }
    System.out.println();
  }

  private static void plotLearningCurve(
    double[][] irisData,
    int numClusters,
    XYChart chart,
    Color color
  ) {
    // Initialize centroids randomly
    double[][] centroids = initializeCentroids(irisData, numClusters);

    // Maximum number of iterations
    int maxIter = 100;

    // Objective function values for each iteration
    double[] objectiveValues = new double[maxIter];

    // K-means algorithm
    for (int iter = 0; iter < maxIter; iter++) {
      // Assign each point to the nearest centroid
      int[] clusters = assignToClusters(irisData, centroids);

      // Update centroids using the learning rule
      for (int k = 0; k < numClusters; k++) {
        List<double[]> clusterPoints = getClusterPoints(irisData, clusters, k);
        if (!clusterPoints.isEmpty()) {
          centroids[k] = calculateNewCentroid(clusterPoints);
        }
      }

      // Calculate and store the objective function value
      objectiveValues[iter] =
        calculateObjectiveFunction(irisData, centroids, clusters);
    }

    // Print objective function values for each iteration
    /*System.out.println(
      "Objective Function Values for K=" +
      numClusters +
      ": " +
      Arrays.toString(objectiveValues)
    );*/

    // Plot the learning curve
    chart
      .addSeries("K=" + numClusters, null, objectiveValues)
      .setMarker(SeriesMarkers.CIRCLE)
      .setLineColor(color);
  }

  private static void plotClustersWithData(
    double[][] irisData,
    int numClusters,
    XYChart chart,
    Color color,
    String seriesName
  ) {
    double[][] centroids = initializeCentroids(irisData, numClusters);
    int maxIter = 100;
    double convergenceThreshold = 0.0001; // Set your desired threshold here

    List<double[][]> intermediateCentroids = new ArrayList<>();
    intermediateCentroids.add(centroids.clone());

    boolean converged = false;

    for (int iter = 0; iter < maxIter && !converged; iter++) {
      int[] clusters = assignToClusters(irisData, centroids);
      double[][] prevCentroids = centroids.clone();

      for (int k = 0; k < numClusters; k++) {
        List<double[]> clusterPoints = getClusterPoints(irisData, clusters, k);
        if (!clusterPoints.isEmpty()) {
          centroids[k] = calculateNewCentroid(clusterPoints);
        }
      }

      intermediateCentroids.add(centroids.clone());

      // Check for convergence
      converged = true;
      for (int i = 0; i < centroids.length; i++) {
        double[] prev = prevCentroids[i];
        double[] current = centroids[i];
        for (int j = 0; j < prev.length; j++) {
          if (Math.abs(prev[j] - current[j]) > convergenceThreshold) {
            converged = false;
            break;
          }
        }
        if (!converged) {
          break;
        }
      }
    }

    plotData(irisData, chart, seriesName);

    plotCentroids(intermediateCentroids.get(0), "Initial", chart, Color.red);
    plotCentroids(
      intermediateCentroids.get(intermediateCentroids.size() / 2),
      "Intermediate",
      chart,
      Color.yellow
    );
    plotCentroids(
      intermediateCentroids.get(intermediateCentroids.size() - 1),
      "Converged",
      chart,
      Color.green
    );
  }

  private static void plotData(
    double[][] irisData,
    XYChart chart,
    String seriesName
  ) {
    double[] xData = new double[irisData.length];
    double[] yData = new double[irisData.length];

    for (int i = 0; i < irisData.length; i++) {
      xData[i] = irisData[i][2]; // Petal Length
      yData[i] = irisData[i][3]; // Petal Width
    }

    chart.addSeries(seriesName, xData, yData);
  }

  private static void plotCentroids(
    double[][] centroids,
    String label,
    XYChart chart,
    Color pcolor
  ) {
    double[] xCentroids = new double[centroids.length];
    double[] yCentroids = new double[centroids.length];
    //printMatrix(centroids);

    for (int i = 0; i < centroids.length; i++) {
      xCentroids[i] = centroids[i][2]; // Petal Length
      yCentroids[i] = centroids[i][3]; // Petal Width
    }
    chart
      .addSeries(
        label + " " + System.currentTimeMillis(),
        xCentroids,
        yCentroids
      )
      .setMarker(SeriesMarkers.CIRCLE)
      .setMarkerColor(pcolor);
  }

  public static void plotDecisionBoundaries(
    double[][] irisData,
    int numClusters
  ) {
    // Fit k-means and obtain centroids
    double[][] centroids = initializeCentroids(irisData, numClusters);
    int maxIter = 100;
    double convergenceThreshold = 0.0001; // Set your desired threshold here

    boolean converged = false;

    //calculatingCentroids
    for (int iter = 0; iter < maxIter && !converged; iter++) {
      int[] clusters = assignToClusters(irisData, centroids);

      double[][] prevCentroids = centroids.clone();

      for (int k = 0; k < numClusters; k++) {
        List<double[]> clusterPoints = getClusterPoints(irisData, clusters, k);
        if (!clusterPoints.isEmpty()) {
          centroids[k] = calculateNewCentroid(clusterPoints);
        }
      }

      // Check for convergence
      converged = true;
      for (int i = 0; i < centroids.length; i++) {
        double[] prev = prevCentroids[i];
        double[] current = centroids[i];
        for (int j = 0; j < prev.length; j++) {
          if (Math.abs(prev[j] - current[j]) > convergenceThreshold) {
            converged = false;
            break;
          }
        }
        if (!converged) {
          break;
        }
      }
    }

    double[][] projectedCentroids = new double[numClusters][2]; // 2D space for Petal Length and Petal Width

    for (int i = 0; i < numClusters; i++) {
      projectedCentroids[i][0] = centroids[i][2]; // Petal Length
      projectedCentroids[i][1] = centroids[i][3]; // Petal Width
    }
    // Create a mesh grid covering the dataset range
    double[][] meshGrid = createMeshGrid(irisData);

    // Predict clusters for each point on the mesh grid
    int[] meshPredictions = assignToClusters(meshGrid, projectedCentroids);
    //printArray(meshPredictions);

    // Prepare data for plotting decision boundaries
    double[] xData = new double[meshGrid.length];
    double[] yData = new double[meshGrid.length];
    int[] clusterAssignments = new int[meshGrid.length];

    for (int i = 0; i < meshGrid.length; i++) {
      xData[i] = meshGrid[i][0];
      yData[i] = meshGrid[i][1];
      clusterAssignments[i] = meshPredictions[i];
    }

    // Plot decision boundaries
    XYChart chart = new XYChartBuilder()
      .width(800)
      .height(600)
      .title("Decision Boundaries for k = " + numClusters + " Q1d")
      .xAxisTitle("Petal Length")
      .yAxisTitle("Petal Width")
      .build();
    chart
      .getStyler()
      .setDefaultSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Scatter);

    // Plotting each cluster separately
    for (int k = 0; k < numClusters; k++) {
      Color color;
      if (k == 0) {
        color = Color.pink;
      } else if (k == 1) {
        color = Color.cyan;
      } else if (k == 2) {
        color = Color.yellow;
      } else color = Color.ORANGE;
      // Filtering xData and yData based on cluster assignments
      List<Double> clusterXList = new ArrayList<>();
      List<Double> clusterYList = new ArrayList<>();

      for (int i = 0; i < xData.length; i++) {
        if (clusterAssignments[i] == k) {
          clusterXList.add(xData[i]);
          clusterYList.add(yData[i]);
        }
      }

      // Convert ArrayLists to arrays for plotting
      double[] clusterX = clusterXList
        .stream()
        .mapToDouble(Double::doubleValue)
        .toArray();
      double[] clusterY = clusterYList
        .stream()
        .mapToDouble(Double::doubleValue)
        .toArray();

      // Plot the series for the cluster
      chart
        .addSeries("Cluster " + (k + 1), clusterX, clusterY)
        .setMarker(SeriesMarkers.CIRCLE)
        .setMarkerColor(color);
    }
    plotData(irisData, chart, "Iris Data");
    plotCentroids(centroids, "centroids", chart, Color.black);
    new SwingWrapper<>(chart).displayChart();
  }

  private static double[][] createMeshGrid(double[][] data) {
    // Finding the minimum and maximum values for Petal Length and Petal Width
    double minPetalLength = Double.MAX_VALUE;
    double maxPetalLength = Double.MIN_VALUE;
    double minPetalWidth = Double.MAX_VALUE;
    double maxPetalWidth = Double.MIN_VALUE;

    for (double[] row : data) {
      double petalLength = row[2];
      double petalWidth = row[3];

      minPetalLength = Math.min(minPetalLength, petalLength);
      maxPetalLength = Math.max(maxPetalLength, petalLength);
      minPetalWidth = Math.min(minPetalWidth, petalWidth);
      maxPetalWidth = Math.max(maxPetalWidth, petalWidth);
    }

    // Define the number of points on the grid for Petal Length and Petal Width
    int numPoints = 100;
    double[][] meshGrid = new double[numPoints * numPoints][2];

    // Create the mesh grid
    double stepPetalLength =
      (maxPetalLength - minPetalLength) / (numPoints - 1);
    double stepPetalWidth = (maxPetalWidth - minPetalWidth) / (numPoints - 1);

    int index = 0;
    for (int i = 0; i < numPoints; i++) {
      for (int j = 0; j < numPoints; j++) {
        double petalLength = minPetalLength + i * stepPetalLength;
        double petalWidth = minPetalWidth + j * stepPetalWidth;

        meshGrid[index][0] = petalLength;
        meshGrid[index][1] = petalWidth;

        index++;
      }
    }
    return meshGrid;
  }
}
