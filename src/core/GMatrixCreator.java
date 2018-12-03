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
                int k = 0;
                gPrev = gCur.clone();
                cur = BigDecimalMatrix.zeroMatrix(size);
                do {
                    prev = cur.clone();
                    cur = cur.add(per(creator.create(1, k)).multiply(powMatrix(gPrev, k)));
                    //System.out.println(per(creator.create(1, k)));
                    k++;
                } while (cur.subtract(prev).squaredEuclidianNorm().compareTo(accuracy) > 0);
                gCur = cur.clone();
            } while (gCur.subtract(gPrev).squaredEuclidianNorm().compareTo(accuracy) > 0);
            MatrixContainer.setG(gCur.clone());
        } catch (CloneNotSupportedException e) {
            System.out.println(e.getMessage());
        }
       // System.out.println(gCur.getWidth());
        return gCur;
    }

    public static BigDecimalMatrix per(BigDecimalMatrix[][] a) {
        int size = a.length * 2;
        BigDecimalMatrix cur = BigDecimalMatrix.zeroMatrix(size);
        for (int i = 0; i < a.length; i++) {
            for(int j = 0; j < a[i].length; j++) {
                BigDecimalMatrix val = a[i][j];
                //System.out.println(a[i][j]);
                for (int k = 0; k <= 1; k++) {
                    for (int r = 0; r <= 1; r++) {
                        cur.setElement(k + i*2, r + j*2, val.getElement(k, r));
                    }
                }
            }
        }
        //System.out.println(cur);
        return cur;
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

