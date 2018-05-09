package kurs;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.Scanner;

/**
 *\
 * @author Rizar
 */
public class BigDecimalMatrix extends AbstractMatrix<BigDecimal> implements Cloneable
{
    private Matrix<BigDecimal> matrix;
    private int scale;

    /**
     * Constructs zero matrix.
     * @param height
     * @param width
     * @param scale
     */
    public BigDecimalMatrix(BigDecimal[][] matr, int scale) {
        matrix = new ArrayMatrix<BigDecimal>(matr);
        this.scale = scale;
    }

    public BigDecimalMatrix(int height, int width, int scale)
    {
        matrix = new ArrayMatrix<BigDecimal>(height, width);
        for (int i = 0; i < height; i++)
            for (int j = 0; j < width; j++)
                matrix.setElement(i, j, BigDecimal.ZERO);
        this.scale = scale;
    }

    /**
     * Constructs scalar matrix.
     * @param size
     * @param scalar
     * @param scale
     */
    public BigDecimalMatrix(int size, BigDecimal scalar, int scale)
    {
        matrix = new ArrayMatrix<BigDecimal>(size, size);
        for (int i = 0; i < matrix.getHeight(); i++)
            for (int j = 0; j < matrix.getWidth(); j++)
                matrix.setElement(i, j, i == j ? scalar : BigDecimal.ZERO);
        this.scale = scale;
    }

    /**
     * If makeCopy, construct matrix shallow copy, else uses matrix data.
     * @param matrix
     * @param scale
     */
    public BigDecimalMatrix(Matrix<BigDecimal> matrix, int scale, boolean makeCopy)
    {
        this.matrix = makeCopy ? new ArrayMatrix<BigDecimal>(matrix) : matrix;
        this.scale = scale;
    }

    /**
     * If makeCopy, construct matrix shallow copy, else uses matrix data.
     * @param matrix
     */
    public BigDecimalMatrix(BigDecimalMatrix matrix, boolean makeCopy)
    {
        this.matrix = makeCopy ? new ArrayMatrix<BigDecimal>(matrix) : matrix;
        this.scale = matrix.scale;
    }

    public static BigDecimalMatrix eCol(int size, BigDecimal value) {
        BigDecimalMatrix col = new BigDecimalMatrix(size, 1, 30);
        for(int i = 0; i < size; i++) {
            col.setElement(i, 0, value);
        }
        return col;
    }

    public static BigDecimalMatrix eRow(int size, BigDecimal value) {
        BigDecimalMatrix col = new BigDecimalMatrix(1, size, 30);
        for(int i = 0; i < size; i++) {
            col.setElement(0, i, value);
        }
        return col;
    }


    public static BigDecimalMatrix identity(int size) {
        BigDecimalMatrix identity = new BigDecimalMatrix(size, size, 30);
        for(int i = 0; i < size; i++) {
            for(int j = 0; j < size; j++) {
                if(i == j) {
                    identity.setElement(i , j, BigDecimal.ONE);
                } else {
                    identity.setElement(i, j, BigDecimal.ZERO);
                }
            }
        }
        return identity;
    }

    public static BigDecimalMatrix zeroMatrix(int size) {
        BigDecimalMatrix identity = new BigDecimalMatrix(size, size, 30);
        for(int i = 0; i < size; i++) {
            for(int j = 0; j < size; j++) {
                identity.setElement(i, j, BigDecimal.ZERO);
            }
        }
        return identity;
    }

    public void setScale(int scale)
    {
        this.scale = scale;
        for (int i = 0; i < getHeight(); i++)
            for (int j = 0; j < getWidth(); j++)
                setElement(i, j, getElement(i, j).setScale(scale));
    }

    public int getScale()
    {
        return scale;
    }

    public BigDecimalMatrix multiply(BigDecimal scalar)
    {
        BigDecimalMatrix res = new BigDecimalMatrix(this, true);
        for (int i = 0; i < res.getHeight(); i++)
            for (int j = 0; j < res.getWidth(); j++)
                res.setElement(i, j, res.getElement(i, j).multiply(scalar).setScale(scale, RoundingMode.HALF_UP));
        return res;
    }

    public BigDecimalMatrix multiply(BigDecimalMatrix matrix)
    {
        if (getWidth() != matrix.getHeight())
            throw new InconsistentMatrixesException(this, matrix);
        BigDecimalMatrix res = new BigDecimalMatrix(getHeight(), matrix.getWidth(), scale);
        for (int i = 0; i < res.getHeight(); i++)
            for (int j = 0; j < res.getWidth(); j++)
                for (int k = 0; k < getWidth(); k++)
                    res.setElement(i, j, res.getElement(i, j).add(getElement(i, k).multiply(matrix.getElement(k, j))).setScale(scale, RoundingMode.HALF_UP));
        return res;
    }

    public BigDecimalMatrix add(BigDecimalMatrix matrix)
    {
        BigDecimalMatrix res = new BigDecimalMatrix(getHeight(), getWidth(), scale);
        for (int i = 0; i < getHeight(); i++)
            for (int j = 0; j < getWidth(); j++)
                res.setElement(i, j, getElement(i, j).add(matrix.getElement(i, j)).setScale(scale, RoundingMode.HALF_UP));
        return res;
    }
    
    public BigDecimalMatrix subtract(BigDecimalMatrix matrix)
    {
        BigDecimalMatrix res = new BigDecimalMatrix(getHeight(), getWidth(), scale);
        for (int i = 0; i < getHeight(); i++)
            for (int j = 0; j < getWidth(); j++)
                res.setElement(i, j, getElement(i, j).subtract(matrix.getElement(i, j)).setScale(scale, RoundingMode.HALF_UP));
        return res;
    }

    public void multiplyRow(int i, BigDecimal value)
    {
        for (int j = 0; j < matrix.getWidth(); j++)
            matrix.setElement(i, j, matrix.getElement(i, j).multiply(value).setScale(scale, RoundingMode.HALF_UP));
    }

    public void multiplyRowAndAddToRow(int i1, int i2, BigDecimal value)
    {
        for (int j = 0; j < matrix.getWidth(); j++)
            matrix.setElement(i2, j, matrix.getElement(i1, j).multiply(value).add(matrix.getElement(i2, j)).setScale(scale, RoundingMode.HALF_UP));
    }

    public void zeroFirstColumn()
    {
        for (int i = 1; i < matrix.getHeight(); i++)
        {
            BigDecimal k = matrix.getElement(i, 1).divide(matrix.getElement(1, 1), scale, RoundingMode.HALF_UP).negate();
            multiplyRowAndAddToRow(1, i, k);
        }
    }

    public BigDecimal cubicNorm()
    {
        BigDecimal max = BigDecimal.ZERO;
        for (int i = 0; i < getHeight(); i++)
        {
            BigDecimal cur = BigDecimal.ZERO;
            for (int j = 0; j < getWidth(); j++)
                cur = cur.add(getElement(i, j).abs());
            if (cur.compareTo(max) == 1)
                max = cur;
        }
        return max;
    }

    public BigDecimal squaredEuclidianNorm()
    {
        BigDecimal res = BigDecimal.ZERO;
        for (int i = 0; i < getHeight(); i++)
            for (int j = 0; j < getWidth(); j++)
                res = res.add(getElement(i, j).multiply(getElement(i, j)));
        return res;
    }

    public int getHeight()
    {
        return matrix.getHeight();
    }

    public int getWidth()
    {
        return matrix.getWidth();
    }

    public BigDecimal getElement(int i, int j)
    {
        return matrix.getElement(i, j);
    }

    public void setElement(int i, int j, BigDecimal value)
    {
        matrix.setElement(i, j, value.setScale(scale, RoundingMode.HALF_UP));
    }

    @Override
    public String toString()
    {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < getHeight(); i++)
        {
            for (int j = 0; j < getWidth(); j++)
            {
                sb.append(getElement(i, j).toPlainString());
                if (j < getWidth())
                    sb.append(" ");
            }
            if (i < getHeight())
                sb.append("\n");
        }
        return sb.toString();
    }

    /**
     * Creates and returns a copy of this object.  The precise meaning
     * of "copy" may depend on the class of the object. The general
     * intent is that, for any object {@code x}, the expression:
     * <blockquote>
     * <pre>
     * x.clone() != x</pre></blockquote>
     * will be true, and that the expression:
     * <blockquote>
     * <pre>
     * x.clone().getClass() == x.getClass()</pre></blockquote>
     * will be {@code true}, but these are not absolute requirements.
     * While it is typically the case that:
     * <blockquote>
     * <pre>
     * x.clone().equals(x)</pre></blockquote>
     * will be {@code true}, this is not an absolute requirement.
     *
     * By convention, the returned object should be obtained by calling
     * {@code super.clone}.  If a class and all of its superclasses (except
     * {@code Object}) obey this convention, it will be the case that
     * {@code x.clone().getClass() == x.getClass()}.
     *
     * By convention, the object returned by this method should be independent
     * of this object (which is being cloned).  To achieve this independence,
     * it may be necessary to modify one or more fields of the object returned
     * by {@code super.clone} before returning it.  Typically, this means
     * copying any mutable objects that comprise the internal "deep structure"
     * of the object being cloned and replacing the references to these
     * objects with references to the copies.  If a class contains only
     * primitive fields or references to immutable objects, then it is usually
     * the case that no fields in the object returned by {@code super.clone}
     * need to be modified.
     *
     * The method {@code clone} for class {@code Object} performs a
     * specific cloning operation. First, if the class of this object does
     * not implement the interface {@code Cloneable}, then a
     * {@code CloneNotSupportedException} is thrown. Note that all arrays
     * are considered to implement the interface {@code Cloneable} and that
     * the return type of the {@code clone} method of an array type {@code T[]}
     * is {@code T[]} where T is any reference or primitive type.
     * Otherwise, this method creates a new instance of the class of this
     * object and initializes all its fields with exactly the contents of
     * the corresponding fields of this object, as if by assignment; the
     * contents of the fields are not themselves cloned. Thus, this method
     * performs a "shallow copy" of this object, not a "deep copy" operation.
     *
     * The class {@code Object} does not itself implement the interface
     * {@code Cloneable}, so calling the {@code clone} method on an object
     * whose class is {@code Object} will result in throwing an
     * exception at run time.
     *
     * @return a clone of this instance.
     * @throws CloneNotSupportedException if the object's class does not
     *                                    support the {@code Cloneable} interface. Subclasses
     *                                    that override the {@code clone} method can also
     *                                    throw this exception to indicate that an instance cannot
     *                                    be cloned.
     * @see Cloneable
     */
    @Override
    public BigDecimalMatrix clone() throws CloneNotSupportedException {
        return new BigDecimalMatrix(matrix, scale, true);
    }

    public static BigDecimalMatrix scan(Scanner sc, int height, int width, int scale)
    {
        BigDecimalMatrix res = new BigDecimalMatrix(height, width, scale);
        for (int i = 0; i < height; i++)
            for (int j = 0; j < width; j++)
                res.setElement(i, j, sc.nextBigDecimal());
        return res;
    }

    public BigDecimalMatrix inverse() {
        BigDecimalMatrix a = new BigDecimalMatrix(this, true);
        if(getWidth() != getHeight()) {
            throw new UnsupportedOperationException("Unsupported inversion for non square matrices.");
        }
        int n = getWidth();
        BigDecimalMatrix x = BigDecimalMatrix.zeroMatrix(n);
        BigDecimalMatrix b = BigDecimalMatrix.identity(n);
        int index[] = new int[n];

        // Transform the matrix into an upper triangle
        a = gaussian(a, index);

        // Update the matrix b[i][j] with the ratios stored
        for (int i=0; i<n-1; ++i) {
            for (int j = i + 1; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    b.setElement(index[j], k, b.getElement(index[j], k).subtract(a.getElement(index[j], i).multiply(b.getElement(index[i], k))));
                }
            }
        }

        // Perform backward substitutions
        for (int i=0; i<n; ++i) {
            x.setElement(n - 1, i, b.getElement(index[n-1], i).divide(a.getElement(index[n-1], n - 1), RoundingMode.HALF_UP));
            for (int j=n-2; j>=0; --j) {
                x.setElement(j, i, b.getElement(index[j], i));
                for (int k=j+1; k<n; ++k) {
                    x.setElement(j, i, x.getElement(j, i).subtract(a.getElement(index[j], k).multiply(x.getElement(k, i))));
                }
                x.setElement(j, i, x.getElement(j, i).divide(a.getElement(index[j], j), RoundingMode.HALF_UP));
            }
        }
        return x;
    }

    private BigDecimalMatrix gaussian(BigDecimalMatrix a, int index[]) {
        int n = index.length;
        BigDecimal c[] = new BigDecimal[n];

        // Initialize the index
        for (int i = 0; i < n; ++i) {
            index[i] = i;
        }

        // Find the rescaling factors, one from each row
        for (int i = 0; i < n; ++i) {
            BigDecimal c1 = BigDecimal.ZERO;
            for (int j = 0; j < n; ++j) {
                BigDecimal c0 = a.getElement(i, j).abs();
                if (c0.compareTo(c1) > 0) {
                    c1 = c0;
                }
            }
            c[i] = c1;
        }

        // Search the pivoting element from each column
        int k = 0;
        for (int j = 0; j < n - 1; ++j) {
            BigDecimal pi1 = BigDecimal.ZERO;
            for (int i = j; i < n; ++i) {
                BigDecimal pi0 = a.getElement(index[i], j).abs();
                pi0 = pi0.divide(c[index[i]], RoundingMode.HALF_UP);
                if (pi0.compareTo(pi1) > 0) {
                    pi1 = pi0;
                    k = i;
                }
            }

            // Interchange rows according to the pivoting order
            int itmp = index[j];
            index[j] = index[k];
            index[k] = itmp;
            for (int i = j + 1; i < n; ++i) {
                BigDecimal pj = a.getElement(index[i], j).divide(a.getElement(index[j], j), RoundingMode.HALF_UP);

                // Record pivoting ratios below the diagonal
                a.setElement(index[i], j, pj);

                // Modify other elements accordingly
                for (int l = j + 1; l < n; ++l) {
                    a.setElement(index[i], l, a.getElement(index[i], l).subtract(pj.multiply(a.getElement(index[j], l))));
                }
            }
        }

        return a;
    }
}
