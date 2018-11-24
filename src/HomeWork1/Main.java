package HomeWork1;

import javafx.util.Pair;

public class Main {

    private static final char INDEL = '-';

    private static double maxScore;
    private static boolean firstRun = true;

    public static void main(String[] args) {
        // write your code here
        String x = "ATAAGGCATTGACCGTATTGCCAA";
        String y = "CCCATAGGTGCGGTAGCC";

        // Global

        String[] global = hirschberg(x, y);
        printSequence(global);

        // Local

        Pair<Integer, Integer>[] indices = findStartAndEndIndices(x, y);
        int xStart = indices[0].getKey();
        int xEnd = indices[1].getKey();
        int yStart = indices[0].getValue();
        int yEnd = indices[1].getValue();

        String[] local = hirschberg(x.substring(xStart, xEnd), y.substring(yStart, yEnd));
        printSequence(local);
    }

    private static void printSequence(String[] global) {
        for (int i = 0; i < global.length; i++) {
            System.out.println(global[i]);
        }
        System.out.println("Score: " + maxScore);
    }

    public static Pair<Integer, Integer>[] findStartAndEndIndices(String x, String y) {

        int columns = y.length() + 1;
        int rows = x.length() + 1;

        Pair<Integer, Integer> endIndex = new Pair<>(0, 0);
        double maxValue = -Double.MAX_VALUE;
        Pair<Integer, Integer> maxStartingPoint = new Pair<>(0, 0);

        double[][] lastTwoRows = new double[2][columns];
        Pair<Integer, Integer>[][] startingPoints = new Pair[2][columns];

        lastTwoRows[0][0] = 0;

        for (int i = 1; i < columns; i++) {
            lastTwoRows[0][i] = 0;
            startingPoints[0][i] = new Pair<>(0, i);
        }

        for (int i = 1; i < rows; i++) {
            lastTwoRows[1][0] = 0;
            for (int j = 1; j < columns; j++) {
                char xChar = x.charAt(i - 1);
                char yChar = y.charAt(j - 1);
                double mScore = lastTwoRows[0][j - 1] + scoringFunction(xChar, yChar);
                double xAndIndel = lastTwoRows[0][j] + scoringFunction(xChar, INDEL);
                double yAndIndel = lastTwoRows[1][j - 1] + scoringFunction(yChar, INDEL);

                if (mScore >= 0 && mScore >= xAndIndel && mScore >= yAndIndel) {
                    lastTwoRows[1][j] = mScore;
                    startingPoints[1][j] = startingPoints[0][j - 1];

                } else if (xAndIndel >= 0 && xAndIndel > mScore && xAndIndel > yAndIndel) {
                    lastTwoRows[1][j] = xAndIndel;
                    startingPoints[1][j] = startingPoints[0][j];
                } else if (yAndIndel >= 0) {
                    lastTwoRows[1][j] = yAndIndel;
                    startingPoints[1][j] = startingPoints[1][j - 1];
                } else {
                    lastTwoRows[1][j] = 0;
                    startingPoints[1][j] = new Pair<>(i, j);
                }

                if (lastTwoRows[1][j] > maxValue) {
                    endIndex = new Pair<>(i, j);
                    maxValue = lastTwoRows[1][j];
                    maxStartingPoint = new Pair<>(startingPoints[1][j].getKey(), startingPoints[1][j].getValue());
                }

            }
            lastTwoRows[0] = lastTwoRows[1].clone();
            startingPoints[0] = startingPoints[1].clone();
        }

        Pair<Integer, Integer>[] answer = new Pair[2];
        answer[0] = maxStartingPoint;
        answer[1] = endIndex;

        maxScore = maxValue;
        return answer;
    }

    /**
     * @param x - Sequence x
     * @param y - Sequence y
     * @return Needleman - wounch matrix last rows of values
     */

    public static double[] NWScore(String x, String y, boolean firstIteration) {

        int columns = y.length() + 1;
        int rows = x.length() + 1;

        double[][] lastTwoRows = new double[2][columns];

        lastTwoRows[0][0] = 0;

        for (int i = 1; i < columns; i++) {
            lastTwoRows[0][i] = lastTwoRows[0][i - 1] + scoringFunction(y.charAt(i - 1), INDEL);
        }

        for (int i = 1; i < rows; i++) {
            lastTwoRows[1][0] = lastTwoRows[0][0] + scoringFunction(x.charAt(i - 1), INDEL);
            for (int j = 1; j < columns; j++) {
                char xChar = x.charAt(i - 1);
                char yChar = y.charAt(j - 1);
                double mScore = lastTwoRows[0][j - 1] + scoringFunction(xChar, yChar);
                double xAndIndel = lastTwoRows[0][j] + scoringFunction(xChar, INDEL);
                double yAndIndel = lastTwoRows[1][j - 1] + scoringFunction(yChar, INDEL);

                lastTwoRows[1][j] = max(mScore, xAndIndel, yAndIndel);
            }
            lastTwoRows[0] = lastTwoRows[1].clone();

        }

        if (firstIteration) {
            maxScore = lastTwoRows[1][columns - 1];
            firstRun = false;
        }
        return lastTwoRows[1].clone();
    }

    // The alignment itself

    public static String[] hirschberg(String x, String y) {
        String[] alignment = new String[] { "", "" };
        // base cases

        if (x.length() == 0) {
            for (int i = 0; i < y.length(); i++) {
                alignment[0] += INDEL;
                alignment[1] += y.charAt(i);

            }
        } else if (y.length() == 0) {
            for (int i = 0; i < x.length(); i++) {
                alignment[1] += INDEL;
                alignment[0] += x.charAt(i);
            }
        } else if (x.length() == 1 || y.length() == 1) {
            return needlemanWunch(x, y);
        } else {
            int xMid = (int) Math.floor(x.length() / 2);

            String xStartToMid = x.substring(0, xMid);
            String xMidToEnd = x.substring(xMid);

            double[] scoreL = NWScore(xStartToMid, y, firstRun);
            double[] scoreR = NWScore(reverse(xMidToEnd), reverse(y), firstRun);

            int yMid = argMax(scoreL, revArray(scoreR));
            String yStartToMid = y.substring(0, yMid);
            String yMidToEnd = y.substring(yMid);

            String[] firstAlignment = hirschberg(xStartToMid, yStartToMid);
            String[] secondAlignment = hirschberg(xMidToEnd, yMidToEnd);

            alignment[0] = firstAlignment[0] + secondAlignment[0];
            alignment[1] = firstAlignment[1] + secondAlignment[1];
        }

        return alignment;
    }

    private static String[] needlemanWunch(String x, String y) {
        int rows = x.length() + 1;
        int columns = y.length() + 1;

        double[][] matrix = new double[rows][columns];
        Pair<Integer, Integer>[][] pMatrix = new Pair[rows][columns];

        matrix[0][0] = 0;
        pMatrix[0][0] = null;
        for (int j = 1; j < columns; j++) {
            matrix[0][j] = matrix[0][j - 1] + scoringFunction(y.charAt(j - 1), INDEL);
            pMatrix[0][j] = new Pair<>(0, j - 1);
        }

        for (int i = 1; i < rows; i++) {
            matrix[i][0] = matrix[i - 1][0] + scoringFunction(x.charAt(i - 1), INDEL);
            pMatrix[i][0] = new Pair<>(i - 1, 0);
        }

        for (int i = 1; i < rows; i++) {
            for (int j = 1; j < columns; j++) {
                char xChar = x.charAt(i - 1);
                char yChar = y.charAt(j - 1);
                double mScore = matrix[i - 1][j - 1] + scoringFunction(xChar, yChar);
                double xAndIndel = matrix[i - 1][j] + scoringFunction(xChar, INDEL);
                double yAndIndel = matrix[i][j - 1] + scoringFunction(yChar, INDEL);

                if (mScore >= xAndIndel && mScore >= yAndIndel) {
                    matrix[i][j] = mScore;
                    pMatrix[i][j] = new Pair<>(i - 1, j - 1);
                } else if (xAndIndel > mScore && xAndIndel > yAndIndel) {
                    matrix[i][j] = xAndIndel;
                    pMatrix[i][j] = new Pair<>(i - 1, j);
                } else {
                    matrix[i][j] = yAndIndel;
                    pMatrix[i][j] = new Pair<>(i, j - 1);
                }

                matrix[i][j] = max(mScore, xAndIndel, yAndIndel);
            }
        }

        String firstSequence = "";
        String secondSequence = "";

        Pair<Integer, Integer> pointer = pMatrix[rows - 1][columns - 1];

        int currentFirstIndex = rows - 1;
        int currentSecondIndex = columns - 1;

        while (pointer != null) {
            int nextFirstIndex = pointer.getKey();
            int nextSecondIndex = pointer.getValue();

            if (currentFirstIndex - 1 == nextFirstIndex && currentSecondIndex - 1 == nextSecondIndex) { // this is a
                                                                                                        // match /
                                                                                                        // mismatch
                firstSequence = x.charAt(currentFirstIndex - 1) + firstSequence;
                secondSequence = y.charAt(currentSecondIndex - 1) + secondSequence;
            } else if (currentFirstIndex == nextFirstIndex && currentSecondIndex - 1 == nextSecondIndex) { // this is an
                                                                                                           // indel
                firstSequence = INDEL + firstSequence;
                secondSequence = y.charAt(currentSecondIndex - 1) + secondSequence;
            } else if (currentFirstIndex - 1 == nextFirstIndex && currentSecondIndex == nextSecondIndex) { // this is an
                                                                                                           // indel
                secondSequence = INDEL + secondSequence;
                firstSequence = x.charAt(currentFirstIndex - 1) + firstSequence;
            }

            pointer = pMatrix[nextFirstIndex][nextSecondIndex];
            currentFirstIndex = nextFirstIndex;
            currentSecondIndex = nextSecondIndex;
        }

        return new String[] { firstSequence, secondSequence };
    }

    private static int argMax(double[] scoreL, double[] scoreR) {
        int index = 0;
        double max = -Double.MAX_VALUE;

        double[] sumOfScores = new double[scoreL.length];
        for (int i = 0; i < scoreL.length; i++) {
            sumOfScores[i] = scoreL[i] + scoreR[i];
            if (sumOfScores[i] > max) {
                max = sumOfScores[i];
                index = i;
            }
        }

        return index;

    }

    public static double scoringFunction(char a, char b) {
        double score;

        if (a == b) {
            score = 2;
        } else if (a == INDEL || b == INDEL) {
            score = -3;
        } else {
            return -2;
        }

        return score;
    }

    public static double max(double a, double b, double c) {
        return Math.max(Math.max(a, b), c);
    }

    public static String reverse(String s) {
        return new StringBuilder(s).reverse().toString();
    }

    public static double[] revArray(double[] arr) {
        for (int i = 0; i < arr.length / 2; i++) {
            double temp = arr[i];
            arr[i] = arr[arr.length - i - 1];
            arr[arr.length - i - 1] = temp;
        }

        return arr;
    }

    public class Tuple<X, Y> {
        public final X x;
        public final Y y;

        public Tuple(X x, Y y) {
            this.x = x;
            this.y = y;
        }
    }
}
