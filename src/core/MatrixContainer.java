package core;

import java.util.HashMap;
import kurs.BigDecimalMatrix;
/**
 * Created by Lenovo on 06.05.2018.
 */
public class MatrixContainer {
    private static BigDecimalMatrix                     g;
    private static HashMap<Integer, BigDecimalMatrix> generators;
    private static HashMap<Integer, BigDecimalMatrix>   pSlashMatrices;
    private static HashMap<Integer, BigDecimalMatrix>   phiMatrices;
    private static HashMap<Integer, BigDecimalMatrix>   piVectors;
    private static HashMap<Integer, BigDecimalMatrix>   kMatrix;
    private static HashMap<Integer, BigDecimalMatrix>   pMatrix;

    static {
        generators      = new HashMap<>();
        pSlashMatrices  = new HashMap<>();
        phiMatrices     = new HashMap<>();
        piVectors       = new HashMap<>();

    }

    private MatrixContainer() {
    }

    public static BigDecimalMatrix getG() {
        return g;
    }

    public static void setG(BigDecimalMatrix g) {
        MatrixContainer.g = g;
    }

    public static HashMap<Integer, BigDecimalMatrix> getGenerators() {
        return generators;
    }

    public static HashMap<Integer, BigDecimalMatrix> getPSlashMatrices() {
        return pSlashMatrices;
    }

    public static HashMap<Integer, BigDecimalMatrix> getPhiMatrices() {
        return phiMatrices;
    }

    public static HashMap<Integer, BigDecimalMatrix> getPiVectors() {
        return piVectors;
    }

    public static HashMap<Integer, BigDecimalMatrix> getKMatrix() {
        return kMatrix;
    }

    public static HashMap<Integer, BigDecimalMatrix> getPMatrix() {
        return pMatrix;
    }


    public static void reInit() {
        g               = null;
        generators      = new HashMap<>();
        pSlashMatrices  = new HashMap<>();
        phiMatrices     = new HashMap<>();
        piVectors       = new HashMap<>();
        kMatrix         = new HashMap<>();
        pMatrix         = new HashMap<>();
    }
}
