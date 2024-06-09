use crate::InterpolationValue;

#[derive(Debug)]
pub enum MatrixValidationError {
    MatrixShapeError,
}


pub(crate) struct TridiagonalMatrix<V: InterpolationValue> {
    upper_diagonal: Vec<V>,
    diagonal: Vec<V>,
    lower_diagonal: Vec<V>,
    size: usize,
}

impl<V: InterpolationValue> TridiagonalMatrix<V> {
    pub fn try_new(upper_diagonal: Vec<V>,
                   diagonal: Vec<V>,
                   lower_diagonal: Vec<V>) -> Result<Self, MatrixValidationError> {
        if !(upper_diagonal.len() == lower_diagonal.len() && upper_diagonal.len() + 1 == diagonal.len()) {
            return Err(MatrixValidationError::MatrixShapeError);
        }
        Ok(Self {
            size: diagonal.len(),
            upper_diagonal,
            diagonal,
            lower_diagonal,
        })
    }

    // Solve Ax = b.
    pub fn try_solve(self, b: &[V]) -> Result<Vec<V>, MatrixValidationError> {
        thomas(
            self.size,
            self.lower_diagonal.as_slice(),
            self.diagonal.as_slice(),
            self.upper_diagonal.as_slice(),
            b,
        )
    }
}

fn thomas<V: InterpolationValue>(matrix_size: usize, lower_diagonal: &[V], diagonal: &[V],
                                 upper_diagonal: &[V], b: &[V]) -> Result<Vec<V>, MatrixValidationError> {
    if b.len() != matrix_size {
        return Err(MatrixValidationError::MatrixShapeError);
    }
    let mut x = b.to_vec();
    let mut scratch = Vec::with_capacity(matrix_size);
    scratch.push(upper_diagonal[0] / diagonal[0]);
    x[0] /= diagonal[0];

    /* loop from 1 to X - 1 inclusive */
    for ix in 1..matrix_size {
        if ix < matrix_size - 1 {
            scratch.push(upper_diagonal[ix] / (diagonal[ix] - lower_diagonal[ix - 1] * scratch[ix - 1]));
        }
        x[ix] = (x[ix] - lower_diagonal[ix - 1] * x[ix - 1]) / (diagonal[ix] - lower_diagonal[ix - 1] * scratch[ix - 1]);
    }

    /* loop from X - 2 to 0 inclusive */
    for ix in (0..matrix_size - 2).rev() {
        let temp = scratch[ix] * x[ix + 1];
        x[ix] -= temp;
    }
    Ok(x)
}