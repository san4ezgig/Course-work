package core;

/**
 * Created by Lenovo on 12.05.2018.
 */

import java.math.BigDecimal;

import kurs.BigDecimalMatrix;


public class GMatrixCreator {

    private GeneratorCreator creator;
    private BigDecimal accuracy;

    public GMatrixCreator(GeneratorCreator creator) {
        this.creator = creator;
        this.accuracy = creator.getAccuracy();
    }

    public BigDecimalMatrix create(BigDecimalMatrix g0) {
        int size = g0.getWidth();
        BigDecimalMatrix gPrev, cur, prev;
        BigDecimalMatrix gCur = new BigDecimalMatrix(size, new BigDecimal(1), 10);
        try {
            do {
                int j = 0;
                gPrev = gCur.clone();
                cur = BigDecimalMatrix.zeroMatrix(size);
                do {
                    prev = cur.clone();
                    cur = cur.add(creator.create(1, j).multiply(powMatrix(gPrev, j)));
                    j++;
                } while (cur.subtract(prev).squaredEuclidianNorm().doubleValue() < accuracy.doubleValue());
                gCur = cur.clone();
            } while (gCur.subtract(gPrev).squaredEuclidianNorm().doubleValue() < accuracy.doubleValue());
            MatrixContainer.setG(gCur.clone());
        } catch (CloneNotSupportedException e) {
            System.out.println(e.getMessage());
        }
        return gCur;
    }

    //Approved
    private BigDecimalMatrix powMatrix(BigDecimalMatrix a, int b) {
        if (b < 0) {
            throw new IllegalArgumentException("Negative argument in factorial function.");
        }
        if (b == 0) {
            return BigDecimalMatrix.identity(a.getWidth());
        }
        return a.multiply(powMatrix(a, b - 1));
    }
}

