use nalgebra::{Point3};
use crate::{HermiteSplineError, InterpolationValue};
use crate::tridiagonal_matrix::TridiagonalMatrix;
use num_traits::Zero;

pub struct NaturalCubicSpline<V: InterpolationValue> {
    points: Vec<Point3<V>>,
}

impl<V: InterpolationValue> NaturalCubicSpline<V> {
    pub fn try_new(raw_points: &[(V, V)]) -> Result<Self, HermiteSplineError<V>> {
        if raw_points.len() < 3 {
            return Err(HermiteSplineError::InsufficientPointsError(
                raw_points.len(),
            ));
        }
        let mut du = Vec::with_capacity(raw_points.len() - 1);
        let mut d = Vec::with_capacity(raw_points.len());
        let mut dl = Vec::with_capacity(raw_points.len() - 1);
        for i in 0..raw_points.len() {
            if i == 0 {
                du.push(V::zero());
                d.push(V::one());
            } else if i + 1 == raw_points.len() {
                d.push(V::one());
                dl.push(V::zero());
            } else {
                let h = raw_points[i].0 - raw_points[i - 1].0;
                let h_next = raw_points[i + 1].0 - raw_points[i].0;
                du.push(h_next / V::from_i8(6).unwrap());
                d.push((h + h_next) / V::from_i8(3).unwrap());
                dl.push(h / V::from_i8(6).unwrap());
            }
        }

        let mut b = Vec::with_capacity(raw_points.len());
        for i in 0..raw_points.len() {
            if i == 0 || i + 1 == raw_points.len() {
                b.push(V::zero());
            } else {
                b.push(
                    (raw_points[i + 1].1 - raw_points[i].1) / (raw_points[i + 1].0 - raw_points[i].0)
                        - (raw_points[i].1 - raw_points[i - 1].1) / (raw_points[i].0 - raw_points[i - 1].0)
                )
            }
        }

        let matrix = TridiagonalMatrix::try_new(du, d, dl).unwrap();
        let derivatives = matrix.try_solve(&b).unwrap();


        let mut temp = raw_points[0].0;
        let mut points = Vec::new();
        for (&(x, y), dydx) in raw_points.iter().zip(derivatives) {
            let point = Point3::new(x, y, dydx);
            if point.x < temp {
                return Err(HermiteSplineError::PointOrderError);
            }
            temp = point.x;
            points.push(point);
        }


        Ok(Self { points })
    }

    pub fn try_value(self, x: V) -> Result<V, HermiteSplineError<V>> {
        match self
            .points
            .binary_search_by(|point| point.x.partial_cmp(&x).unwrap())
        {
            Ok(pos) => Ok(self.points[pos].y),
            Err(pos) => {
                if pos.is_zero() {
                    return Err(HermiteSplineError::OutOfLowerBound(x));
                }
                if pos > self.points.len() {
                    return Err(HermiteSplineError::OutOfUpperBound(x));
                }
                let pos = pos - 1;
                let point = &self.points[pos];
                let next_point = &self.points[pos + 1];
                let h = next_point.x - point.x;
                let six = V::from_i8(6).unwrap();
                Ok(
                    (next_point.x - x) * (next_point.x - x) * (next_point.x - x) / six / h * point.z
                        + (x - point.x) * (x - point.x) * (x - point.x) / six / h * next_point.z
                        + (next_point.x - x) * (point.y / h - h / six * point.z)
                        + (x - point.x) * (next_point.y / h - h / six * next_point.z)
                )
            }
        }
    }
}

#[cfg(test)]
mod tests {
    #[cfg(feature = "decimal")]
    use rust_decimal::Decimal;

    use crate::natural_cubic_spline::NaturalCubicSpline;


    #[test]
    fn test_f64() {
        let points = [(0.0, 1.0), (0.5, 0.5), (1.0, 0.0)];
        let interpolator = NaturalCubicSpline::try_new(&points).unwrap();
        let val = interpolator.try_value(0.75).unwrap();
        assert_eq!(val, 0.25_f64);
    }

    #[cfg(feature = "decimal")]
    #[test]
    fn test_decimal() {
        let points = [
            (Decimal::new(0, 0), Decimal::new(1, 0)),
            (Decimal::new(5, 1), Decimal::new(5, 1)),
            (Decimal::new(1, 0), Decimal::new(0, 0)),
        ];
        let interpolator = NaturalCubicSpline::try_new(&points).unwrap();
        let val = interpolator.try_value(Decimal::new(75, 2)).unwrap();
        assert_eq!(val, Decimal::from_str_exact("0.25").unwrap());
    }
}