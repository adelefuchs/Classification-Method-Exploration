import ChartDirector.*;
import com.opencsv.CSVReader;
import com.opencsv.exceptions.CsvException;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import javax.imageio.ImageIO;
import org.knowm.xchart.*;
import org.knowm.xchart.style.markers.SeriesMarkers;

public class linearDecisionBoundaries {

  public static void main(String[] args) {
    double[] weights = { -.2, -.1 }; // Adjust weights
    double bias = 1; // Adjust bias
    double[] weights2 = { -0.5, -0.5 };
    double bias2 = 1;

    double[][] irisData = new double[0][0];
    String[] irisClassification = new String[0];
    double[] petalLen = new double[0];
    double[] petalWid = new double[0];
    int numClusters1 = 2;
    int numClusters2 = 3;
    String csvFile = "irisdata.csv";

    try (CSVReader reader = new CSVReader(new FileReader(csvFile))) {
      List<String[]> data = reader.readAll();

      String[] labels = data.remove(0);

      irisData = new double[100][data.get(0).length - 1];
      irisClassification = new String[100];
      petalLen = new double[100];
      petalWid = new double[100];

      for (int i = 0; i < 100; i++) {
        String[] row = data.get(i + 50);
        irisClassification[i] = row[row.length - 1];
        petalLen[i] = Double.parseDouble(row[row.length - 3]);
        petalWid[i] = Double.parseDouble(row[row.length - 2]);
        for (int j = 0; j < row.length - 1; j++) {
          irisData[i][j] = Double.parseDouble(row[j]);
        }
      }
    } catch (IOException | CsvException | NumberFormatException e) {
      e.printStackTrace();
    }
    //2a
    org.knowm.xchart.XYChart plainChart = new XYChartBuilder()
      .width(800)
      .height(600)
      .title("Iris Data")
      .xAxisTitle("Petal Length")
      .yAxisTitle("Petal Width")
      .build();
    plainChart
      .getStyler()
      .setDefaultSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Scatter);
    plotData(
      irisData,
      irisClassification,
      plainChart,
      new String[] { "versicolor", "virginica" }
    );
    new SwingWrapper<>(plainChart).displayChart();

    //2bc
    org.knowm.xchart.XYChart boundaryChart = new XYChartBuilder()
      .width(800)
      .height(600)
      .title("Decision Boundary")
      .xAxisTitle("Petal Length")
      .yAxisTitle("Petal Width")
      .build();
    boundaryChart
      .getStyler()
      .setDefaultSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Line);

    // Plot the decision boundary on the new chart
    double[] lineBounds = plotDecisionBoundary(
      irisData,
      irisClassification,
      boundaryChart,
      weights,
      bias,
      "Decision Boundary"
    );

    // Display the chart with the decision boundary
    new SwingWrapper<>(boundaryChart).displayChart();

    // Create a 3D surface chart
    SurfaceChart surfaceChart = new SurfaceChart(800, 600);

    double[] functionOut = calculateAllOutputs(
      petalLen,
      petalWid,
      weights,
      bias
    );
    // Set the data for the 3D surface chart - 3 arrays, x -petal len,y - petal width,z - output function
    surfaceChart.setData(petalLen, petalWid, functionOut);

    // Set the x-axis, y-axis, and z-axis labels
    surfaceChart.xAxis().setTitle("Petal Length");
    surfaceChart.yAxis().setTitle("Petal Width");
    surfaceChart.zAxis().setTitle("Output");
    // Set the title of the chart
    surfaceChart.addTitle("Neural Network Output Surface");

    String javaFilePath =
      linearDecisionBoundaries.class.getProtectionDomain()
        .getCodeSource()
        .getLocation()
        .getPath();
    String parentDirectoryPath = new File(javaFilePath).getParent();
    String project2FolderPath = parentDirectoryPath;

    saveChartAsImage(surfaceChart, project2FolderPath, "3dsurface.png");

    int[] classints = new int[irisClassification.length];
    //covert classes to ints
    for (int i = 0; i < irisClassification.length; i++) {
      if (irisClassification[i].equals("setosa")) {
        classints[i] = 0;
      } else if (irisClassification[i].equals("versicolor")) {
        classints[i] = 1;
      } else {
        classints[i] = 2;
      }
    }
    //3B
    // Calculate mean squared error using the provided function
    double mse = calculateMSE(irisData, classints, weights, bias);
    System.out.println("Mean Squared Error: " + mse);
    double mse2 = calculateMSE(irisData, classints, weights2, bias2);
    System.out.println("Mean Squared Error: " + mse2);

    org.knowm.xchart.XYChart chart = new XYChartBuilder()
      .width(800)
      .height(600)
      .title("Iris Data")
      .xAxisTitle("Petal Length")
      .yAxisTitle("Petal Width")
      .build();
    chart
      .getStyler()
      .setDefaultSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Line);

    plotData(
      irisData,
      irisClassification,
      chart,
      new String[] { "versicolor", "virginica" }
    );
    // Define x-values for the decision boundary range
    double[] xValues = { 4, 7 };

    // Calculate y-values for decision boundary 1 using weights1 and bias1
    double[] yValues1 = new double[xValues.length];
    double[] lals = new double[xValues.length];
    for (int i = 0; i < xValues.length; i++) {
      // Calculate y = (-w1 * x - b) / w2 for each x-value
      yValues1[i] = (-weights[0] * xValues[i] - bias) / weights[1];
      lals[i] = computeOutput(xValues, weights, bias);
    }

    // Calculate y-values for decision boundary 2 using weights2 and bias2
    double[] yValues2 = new double[xValues.length];
    for (int i = 0; i < xValues.length; i++) {
      // Calculate y = (-w1 * x - b) / w2 for each x-value
      yValues2[i] = (-weights2[0] * xValues[i] - bias2) / weights2[1];
    }

    // Plot decision boundary 1
    chart
      .addSeries("Decision Boundary 1", xValues, yValues1)
      .setXYSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Line)
      .setLineColor(Color.RED);

    // Plot decision boundary 2
    chart
      .addSeries("Decision Boundary 2", xValues, yValues2)
      .setXYSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Line)
      .setLineColor(Color.RED);
    // Show the plot
    new SwingWrapper<>(chart).displayChart();

    //2e plotting what my model predicted
    org.knowm.xchart.XYChart predChart = new XYChartBuilder()
      .width(800)
      .height(600)
      .title("Predicted Iris Data Classification with Linear Decision Boundary")
      .xAxisTitle("Petal Length")
      .yAxisTitle("Petal Width")
      .build();
    predChart
      .getStyler()
      .setDefaultSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Scatter);

    String[] predictedClasses = new String[150];
    predictedClasses =
      predictPointsWDB(
        lineBounds,
        irisData,
        irisClassification,
        predChart,
        weights,
        bias,
        petalLen,
        petalWid
      );
    plotData(
      irisData,
      predictedClasses,
      predChart,
      new String[] { "Versicolor", "Virginica" }
    );
    new SwingWrapper<>(predChart).displayChart();

    //3e
    double learningRate = 0.001;
    int numIterations = 5;
    org.knowm.xchart.XYChart gradientStepChart = new XYChartBuilder()
      .width(800)
      .height(600)
      .title("Step of Gradient on Iris Data")
      .xAxisTitle("Petal Length")
      .yAxisTitle("Petal Width")
      .build();
    gradientStepChart
      .getStyler()
      .setDefaultSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Line);

    plotData(
      irisData,
      irisClassification,
      chart,
      new String[] { "versicolors", "virginicas" }
    );

    updateAndPlot(
      irisData,
      classints,
      weights,
      bias,
      learningRate,
      numIterations,
      gradientStepChart
    );
    new SwingWrapper<>(gradientStepChart).displayChart();

    //4ab
    double learningRate2 = 0.000001;
    org.knowm.xchart.XYChart dataAndBoundaryChart = new XYChartBuilder()
      .width(800)
      .height(600)
      .title("Iris Data with Decision Boundary")
      .xAxisTitle("Petal Length")
      .yAxisTitle("Petal Width")
      .build();
    dataAndBoundaryChart
      .getStyler()
      .setDefaultSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Scatter);

    org.knowm.xchart.XYChart learningCurveChart = new XYChartBuilder()
      .width(800)
      .height(600)
      .title("Learning Curve")
      .xAxisTitle("Iterations")
      .yAxisTitle("Mean Squared Error")
      .build();
    learningCurveChart
      .getStyler()
      .setDefaultSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Line);

    double convergenceThreshold = .00005;
    updateAndPlotWithGradientDescent(
      irisData,
      classints,
      weights,
      bias,
      learningRate2,
      convergenceThreshold,
      dataAndBoundaryChart,
      learningCurveChart
    );

    new SwingWrapper<>(dataAndBoundaryChart).displayChart();
    new SwingWrapper<>(learningCurveChart).displayChart();
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

  private static void plotData(
    double[][] irisData,
    String[] irisClassification,
    org.knowm.xchart.XYChart chart,
    String[] seriesName
  ) {
    double[] xData1 = new double[irisData.length];
    double[] yData1 = new double[irisData.length];
    double[] xData2 = new double[irisData.length];
    double[] yData2 = new double[irisData.length];

    int index1 = 0;
    int index2 = 0;

    for (int i = 0; i < irisData.length; i++) {
      if (irisData[i][2] == 0 || irisData[i][3] == 0) {} else if (
        irisClassification[i].equals("versicolor")
      ) {
        xData1[index1] = irisData[i][2]; // Petal Length
        yData1[index1] = irisData[i][3]; // Petal Width
        index1++;
      } else if (irisClassification[i].equals("virginica")) {
        xData2[index2] = irisData[i][2]; // Petal Length
        yData2[index2] = irisData[i][3]; // Petal Width
        index2++;
      }
    }
    chart
      .addSeries(seriesName[0], xData1, yData1)
      .setMarkerColor(Color.green)
      .setMarker(SeriesMarkers.CIRCLE); // Set marker for data points
    chart
      .addSeries(seriesName[1], xData2, yData2)
      .setMarkerColor(Color.blue)
      .setMarker(SeriesMarkers.CIRCLE); // Set marker for data points

    chart
      .getStyler()
      .setDefaultSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Scatter); // Set scatter plot style
  }

  public static double sigmoid(double z) {
    return 1.0 / (1.0 + Math.exp(-z));
  }

  public static double computeOutput(
    double[] inputFeatures,
    double[] weights,
    double bias
  ) {
    if (inputFeatures.length != weights.length) {
      throw new IllegalArgumentException(
        "Input features and weights must have the same length"
      );
    }

    double z = bias;
    for (int i = 0; i < inputFeatures.length; i++) {
      z += inputFeatures[i] * weights[i];
    }

    return sigmoid(z);
  }

  public static double[] plotDecisionBoundary(
    double[][] irisData,
    String[] irisClassification,
    org.knowm.xchart.XYChart chart,
    double[] weights,
    double bias,
    String seriesName
  ) {
    // Define x-values for the decision boundary range
    double[] xValues = { 0.1, 7 }; //idk what these should be maybe the range of the data? but of which len or wid idk

    // Calculate y-values for decision boundary 1 using weights1 and bias1
    double[] yValues1 = new double[xValues.length];
    for (int i = 0; i < xValues.length; i++) {
      // Calculate y = (-w1 * x - b) / w2 for each x-value
      yValues1[i] = (-weights[0] * xValues[i] - bias) / weights[1];
      //yValues1[i] = computeOutput(xValues, weights, bias);
    }

    chart
      .addSeries(
        seriesName,
        new double[] { xValues[0], yValues1[0] },
        new double[] { xValues[1], yValues1[1] }
      )
      .setXYSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Line)
      .setLineColor(Color.RED);

    // Overlay the iris data points
    plotData(
      irisData,
      irisClassification,
      chart,
      new String[] { seriesName + " versicolor", seriesName + " virginica" }
    );
    double[] rtrn = new double[] {
      xValues[0],
      yValues1[0],
      xValues[1],
      yValues1[1],
    };
    return rtrn;
  }

  public static double[] calculateAllOutputs(
    double[] petalLen,
    double[] petalWidth,
    double[] weights,
    double bias
  ) {
    if (petalLen.length != petalWidth.length) {
      throw new IllegalArgumentException(
        "Petal length and width arrays must have the same length"
      );
    }

    double[] funcOut = new double[petalLen.length];

    for (int i = 0; i < petalLen.length; i++) {
      double[] inputFeatures = { petalLen[i], petalWidth[i] };
      for (int j = 0; j < inputFeatures.length; j++) {
        //funcOut[i] += (-weights[0] * inputFeatures[j] - bias) / weights[1];
        funcOut[i] = computeOutput(inputFeatures, weights, bias);
      }
    }

    return funcOut;
  }

  private static void saveChartAsImage(
    BaseChart chart,
    String directoryPath,
    String filename
  ) {
    try {
      // Concatenate directory path and filename
      String filePath = directoryPath + File.separator + filename;

      // Output the chart to a PNG file
      chart.makeChart(filePath);
      System.out.println("Chart created: " + filePath);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  //3a
  public static double calculateMSE(
    double[][] irisData,
    int[] patternClasses,
    double[] weights,
    double bias
  ) {
    double totalError = 0.0;
    int numSamples = irisData.length;

    for (int i = 0; i < numSamples; i++) {
      double[] input = { irisData[i][2], irisData[i][3] };
      int targetClass = patternClasses[i];

      // Compute the output of the neural network for the input
      double output = (computeOutput(input, weights, bias));

      // Convert target class to expected output
      double expectedOutput = targetClass == 0 ? 0.0 : 1.0;

      // Calculate the error as the squared difference between predicted and expected output
      double error = output - expectedOutput;

      // Square the error and accumulate
      totalError += error * error;
    }

    // Calculate mean squared error
    double meanSquaredError = totalError / numSamples;
    return meanSquaredError;
  }

  public static void updateAndPlot(
    double[][] data,
    int[] classes,
    double[] weights,
    double bias,
    double learningRate,
    int numIterations,
    org.knowm.xchart.XYChart chart
  ) {
    weights[0] -= .3;
    weights[1] -= .2;
    for (int iter = 0; iter < numIterations; iter++) {
      double gradientW1 = 0.0;
      double gradientW2 = 0.0;
      double gradientBias = 0.0;

      for (int i = 0; i < data.length; i++) {
        double[] input = { data[i][0], data[i][1] };
        int targetClass = classes[i];

        double output = computeOutput(input, weights, bias);
        double error = (output - (targetClass == 0 ? 0.0 : 1.0));

        gradientW1 += error * output * (1 - output) * input[0];
        gradientW2 += error * output * (1 - output) * input[1];
        gradientBias += error * output * (1 - output);
      }
      //System.out.println("Weights: " + weights[0] + ", " + weights[1]);
      //System.out.println("Bias: " + bias);
      weights[0] -= (learningRate * gradientW1);
      weights[1] -= (learningRate * gradientW2);
      bias -= learningRate * gradientBias;
      String[] stringClasses = new String[classes.length];
      for (int i = 0; i < classes.length; i++) {
        if (classes[i] == 1) {
          stringClasses[i] = "versicolor";
        } else {
          stringClasses[i] = "virginica";
        }
      }

      plotDecisionBoundary(
        data,
        stringClasses,
        chart,
        weights,
        bias,
        "iteration " + iter
      );
    }
  }

  public static double calculateMean(double[] data) {
    double sum = 0;
    for (double point : data) {
      sum += point; // Assuming the mean is calculated for the first value (petal length)
    }
    return sum / data.length;
  }

  public static String[] predictPointsWDB(
    double[] lineBounds,
    double[][] irisData,
    String[] irisClassification,
    org.knowm.xchart.XYChart chart,
    double[] weights,
    double bias,
    double[] petalLen,
    double[] petalWid
  ) {
    String[] predictions = new String[irisData.length];
    //double[] lineBounds  = plotDecisionBoundary(irisData, irisClassification,chart,weights,bias);
    double[] slopeForm = getLineEquation(lineBounds);
    for (int i = 0; i < irisData.length; i++) {
      predictions[i] =
        classifyPoint(new double[] { petalLen[i], petalWid[i] }, slopeForm);
    }
    return predictions;
  }

  public static double[] getLineEquation(double[] coordinates) {
    // Extracting coordinates from the array
    double x1 = coordinates[0];
    double y1 = coordinates[1];
    double x2 = coordinates[2];
    double y2 = coordinates[3];

    // Calculate the slope
    double slope = (y2 - y1) / (x2 - x1);

    // Calculate the y-intercept (b) using one of the points and the slope (y = mx + b)
    double yIntercept = y1 - slope * x1;

    // Store the slope and y-intercept values in an array
    double[] equation = { slope, yIntercept };

    return equation;
  }

  public static String classifyPoint(double[] point, double[] lineEquation) {
    double x0 = point[0];
    double y0 = point[1];
    double m = lineEquation[0];
    double b = lineEquation[1];

    double calculatedY;
    if (Double.isFinite(m)) {
      // For non-vertical lines
      calculatedY = m * x0 + b;
    } else {
      // For vertical lines
      calculatedY = b;
    }

    if (y0 > calculatedY) {
      return "virginica"; // Above the line
    } else {
      return "versicolor"; // Below or on the line
    }
  }

  public static void updateAndPlotWithGradientDescent(
    double[][] data,
    int[] classes,
    double[] weights,
    double bias,
    double learningRate,
    double convergenceThreshold,
    org.knowm.xchart.XYChart dataAndBoundaryChart,
    org.knowm.xchart.XYChart learningCurveChart
  ) {
    String[] stringClasses = new String[classes.length];
    for (int i = 0; i < classes.length; i++) {
      if (classes[i] == 1) {
        stringClasses[i] = "versicolor";
      } else {
        stringClasses[i] = "virginica";
      }
    }
    Random random = new Random();
    random.setSeed(457);
    double[] randWeights = { random.nextDouble(), random.nextDouble() };
    double randbias = random.nextDouble();
    int iter = 0;
    List<Double> mseList = new ArrayList<>();
    double prevMSE = Double.MAX_VALUE; // Initialize with a large value

    while (true) {
      iter++;
      double gradientW1 = 0.0;
      double gradientW2 = 0.0;
      double gradientBias = 0.0;

      for (int i = 0; i < data.length; i++) {
        double[] input = { data[i][0], data[i][1] };
        int targetClass = classes[i];

        double output = computeOutput(input, weights, bias);
        double error = output - (targetClass == 0 ? 0.0 : 1.0);

        gradientW1 += error * output * (1 - output) * input[0];
        gradientW2 += error * output * (1 - output) * input[1];
        gradientBias += error * output * (1 - output);
      }

      weights[0] -= learningRate * gradientW1;
      weights[1] -= learningRate * gradientW2;
      bias -= learningRate * gradientBias;

      double mse = calculateMSE(data, classes, weights, bias);
      mseList.add(mse);

      // Check convergence
      if (Math.abs(prevMSE - mse) < convergenceThreshold || iter > 1650) {
        System.out.println("Convergence reached at iteration: " + iter);
        break; // Exit the loop if convergence criteria are met
      }
      prevMSE = mse;
      if (iter == 1 || iter % 500 == 0 || iter == 1650) {
        plotDecisionBoundary2(
          data,
          stringClasses,
          dataAndBoundaryChart,
          weights,
          bias,
          "iteration " + iter
        );
      }
    }
    // Plot learning curve
    double[] iterations = new double[mseList.size()];
    for (int i = 0; i < mseList.size(); i++) {
      iterations[i] = i + 1;
    }
    learningCurveChart
      .addSeries(
        "Learning Curve",
        iterations,
        mseList.stream().mapToDouble(Double::doubleValue).toArray()
      )
      .setXYSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Line)
      .setLineColor(Color.BLUE);
  }

  public static double[] plotDecisionBoundary2(
    double[][] irisData,
    String[] irisClassification,
    org.knowm.xchart.XYChart chart,
    double[] weights,
    double bias,
    String seriesName
  ) {
    // Define x-values for the decision boundary range
    double[] xValues = { 1, 5 }; //idk what these should be maybe the range of the data? but of which len or wid idk

    // Calculate y-values for decision boundary 1 using weights1 and bias1
    double[] yValues1 = new double[xValues.length];
    for (int i = 0; i < xValues.length; i++) {
      // Calculate y = (-w1 * x - b) / w2 for each x-value
      yValues1[i] = (-weights[0] * xValues[i] - bias) / weights[1];
      //yValues1[i] = computeOutput(xValues, weights, bias);
    }

    chart
      .addSeries(
        seriesName,
        new double[] { xValues[0], yValues1[0] },
        new double[] { xValues[1], yValues1[1] - 2 }
      )
      .setXYSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Line)
      .setLineColor(Color.RED);

    // Overlay the iris data points
    plotData(
      irisData,
      irisClassification,
      chart,
      new String[] { seriesName + " versicolor", seriesName + " virginica" }
    );
    double[] rtrn = new double[] {
      xValues[0],
      yValues1[0],
      xValues[1],
      yValues1[1],
    };
    return rtrn;
  }
}
