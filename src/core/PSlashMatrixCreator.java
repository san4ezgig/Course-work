package core;

import kurs.BigDecimalMatrix;

import java.math.BigDecimal;
import java.util.Objects;


public class PSlashMatrixCreator {
    private GeneratorCreator creator;
    private BigDecimalMatrix g;
    private BigDecimal accuracy;

    //Approved
    public PSlashMatrixCreator(GeneratorCreator creator, BigDecimalMatrix g) {
        this.g = g;
        this.creator = creator;
        this.accuracy = creator.getAccuracy();
    }

    //Approved

    /*public BigDecimalMatrix create(Integer i, Integer j) {
        if(j < i) {
            throw new IllegalArgumentException("Index i is bigger than index j in PSlashMatrixCreator.");
        }

        int l = j + 1;
        int size = g.getWidth();
        BigDecimalMatrix result = GMatrixCreator.per(creator.create(i, j));
        BigDecimalMatrix val = BigDecimalMatrix.zeroMatrix(size);
        BigDecimalMatrix sum = BigDecimalMatrix.zeroMatrix(size);
        do {
            val = GMatrixCreator.per(creator.create(i, l)).multiply(powMatrix(g, l - j));
            sum.add(val);

            l++;
        } while (Math.sqrt(val.squaredEuclidianNorm().doubleValue()) >= accuracy.doubleValue());
        return result.add(sum);
    }*/

    public BigDecimalMatrix create(Integer i, Integer j) {
        if(j < i) {
            throw new IllegalArgumentException("Index i is bigger than index j in PSlashMatrixCreator.");
        }

        int l = j + 1;
        int size = g.getWidth();
        int hash = i > 1 ? Objects.hash(1, j - i + 1) : Objects.hash(0, j);

        BigDecimalMatrix result = BigDecimalMatrix.zeroMatrix(size);

        try {
            if(MatrixContainer.getPSlashMatrices().containsKey(hash)) {
                result = MatrixContainer.getPSlashMatrices().get(hash).clone();
            } else {
                BigDecimalMatrix prev;
                BigDecimalMatrix cur = BigDecimalMatrix.zeroMatrix(size);
                do {
                    prev = cur.clone();
                    cur = cur.add(GMatrixCreator.per(creator.create(i, l)).multiply(powMatrix(g, l - j)));
                    l++;
                } while (cur.subtract(prev).squaredEuclidianNorm().compareTo(accuracy) > 0);
                result = GMatrixCreator.per(creator.create(i, j)).add(cur);
                MatrixContainer.getPSlashMatrices().put(hash, result.clone());
            }
        } catch (CloneNotSupportedException e) {
            System.out.println(e.getMessage());
        }
        return result;
    }

    //Approved
    public BigDecimal getAccuracy() {
        return accuracy;
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

