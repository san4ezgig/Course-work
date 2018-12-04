package kurs;

import java.math.BigDecimal;

/**
 * Created by Lenovo on 04.12.2018.
 */
public class BigDecimalBlockMatrix extends AbstractMatrix<BigDecimalMatrix> implements Cloneable {
    private Matrix<BigDecimalMatrix> matrix;
    private int scale;

    public BigDecimalBlockMatrix(BigDecimalMatrix[][] matr) {
        matrix = new ArrayMatrix<BigDecimalMatrix>(matr);
    }

    public BigDecimalBlockMatrix(int size, BigDecimalMatrix matr, int scale)
    {
        matrix = new ArrayMatrix<BigDecimalMatrix>(size, size);
        for (int i = 0; i < matrix.getHeight(); i++)
            for (int j = 0; j < matrix.getWidth(); j++)
                matrix.setElement(i, j, i == j ? matr : new BigDecimalMatrix(2, 2, 6));
        this.scale = scale;
    }

    public BigDecimalBlockMatrix(Matrix<BigDecimalMatrix> matrix, int scale, boolean makeCopy)
    {
        this.matrix = makeCopy ? new ArrayMatrix<BigDecimalMatrix>(matrix) : matrix;
        this.scale = scale;
    }

    @Override
    public BigDecimalBlockMatrix clone() throws CloneNotSupportedException {
        return new BigDecimalBlockMatrix(matrix, 6, true);
    }
}
