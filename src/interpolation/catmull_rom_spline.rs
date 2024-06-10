use crate::HermiteSplineError;
use crate::InterpolationValue;
use nalgebra::{Matrix4, Vector4};
use num_traits::Zero;
use std::ops::Mul;

struct Point2<V> {
    pub x: V,
    pub y: V,
}

pub struct CatmullRomSpline<V: InterpolationValue> {
    points: Vec<Point2<V>>,
}

impl<V: InterpolationValue> CatmullRomSpline<V> {
    /// Constructs a new CatmullRomSpline from a slice of raw points.
    ///
    /// # Arguments
    ///
    /// * `raw_points` - A slice of tuples containing raw points `(x, y)` where `x` is the x-coordinate and `y` is the y-coordinate.
    ///
    /// # Returns
    ///
    /// * `Result<Self, HermiteSplineError<V>>` - A `Result` that either contains the constructed `HermiteSpline` or an `HermiteSplineError`.
    ///
    /// # Errors
    ///
    /// * `HermiteSplineError::InsufficientPointsError(n)` - If the number of `raw_points` is less than 3, where `n` is the number of raw_points.
    /// * `HermiteSplineError::PointOrderError` - If the x-coordinates of the `raw_points` are not in ascending order.
    ///
    /// # Example
    ///
    /// ```
    /// use spline_interpolation::interpolation::catmull_rom_spline::CatmullRomSpline;
    ///
    /// let raw_points = [(0.0, 0.0), (1.0, 1.0), (2.0, 0.0)];
    /// let spline = CatmullRomSpline::try_new(&raw_points);
    /// assert!(spline.is_ok());
    /// ```
    pub fn try_new(raw_points: &[(V, V)]) -> Result<Self, HermiteSplineError<V>> {
        if raw_points.len() < 3 {
            return Err(HermiteSplineError::InsufficientPointsError(
                raw_points.len(),
            ));
        }
        let mut temp = raw_points[0].0;
        let mut points = Vec::new();
        for &(x, y) in raw_points {
            let point = Point2 { x, y };
            if point.x < temp {
                return Err(HermiteSplineError::PointOrderError);
            }
            temp = point.x;
            points.push(point);
        }
        Ok(Self { points })
    }
    /// Tries to find the value `x` in the Hermite spline.
    ///
    /// # Arguments
    ///
    /// * `x`: The value to find in the Hermite spline.
    ///
    /// # Returns
    ///
    /// * `Ok(V)`: If the value `x` is found in the Hermite spline, returns the corresponding value `V`.
    /// * `Err(HermiteSplineError<V>)`: If the value `x` is not found, returns an error indicating whether `x` is out of the lower or upper bound of the spline.
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
                let delta = (x - point.x) / h;
                let delta2 = delta * delta;
                let delta3 = delta2 * delta;
                let d = Vector4::new(delta3, delta2, delta, V::one());
                Ok((d.transpose()
                    * if pos == 0 {
                        let next_next_point = &self.points[pos + 2];
                        let next_h = next_next_point.x - next_point.x;
                        let beta = h / (h + next_h);
                        Matrix4::new(
                            V::zero(),
                            V::one() - beta,
                            -V::one(),
                            beta,
                            V::zero(),
                            -V::one() + beta,
                            V::one(),
                            -beta,
                            V::zero(),
                            -V::one(),
                            V::one(),
                            V::zero(),
                            V::zero(),
                            V::one(),
                            V::zero(),
                            V::zero(),
                        )
                        .mul(Vector4::new(
                            V::zero(),
                            point.y,
                            next_point.y,
                            next_next_point.y,
                        ))
                    } else if pos + 2 == self.points.len() {
                        let prev_point = &self.points[pos - 1];
                        let prev_h = next_point.x - prev_point.x;
                        let alpha = h / (h + prev_h);
                        Matrix4::new(
                            -alpha,
                            V::one(),
                            -V::one() * alpha,
                            V::zero(),
                            V::from_i8(2).unwrap() * alpha,
                            V::from_i8(-2).unwrap(),
                            V::from_i8(2).unwrap() - V::from_i8(2).unwrap() * alpha,
                            V::zero(),
                            -alpha,
                            V::zero(),
                            alpha,
                            V::zero(),
                            V::zero(),
                            V::one(),
                            V::zero(),
                            V::zero(),
                        )
                        .mul(Vector4::new(
                            prev_point.y,
                            point.y,
                            next_point.y,
                            V::zero(),
                        ))
                    } else {
                        let prev_point = &self.points[pos - 1];
                        let prev_h = next_point.x - prev_point.x;
                        let alpha = h / (h + prev_h);
                        let next_next_point = &self.points[pos + 2];
                        let next_h = next_next_point.x - next_point.x;
                        let beta = h / (h + next_h);
                        Matrix4::new(
                            -alpha,
                            V::from_i8(2).unwrap() - beta,
                            V::from_i8(-2).unwrap() + alpha,
                            beta,
                            V::from_i8(2).unwrap() * alpha,
                            beta - V::from_i8(3).unwrap(),
                            V::from_i8(3).unwrap() - V::from_i8(2).unwrap() * alpha,
                            -beta,
                            -alpha,
                            V::zero(),
                            alpha,
                            V::zero(),
                            V::zero(),
                            V::one(),
                            V::zero(),
                            V::zero(),
                        )
                        .mul(Vector4::new(
                            prev_point.y,
                            point.y,
                            next_point.y,
                            next_next_point.y,
                        ))
                    })
                .x)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    #[cfg(feature = "decimal")]
    use rust_decimal::Decimal;

    use crate::interpolation::catmull_rom_spline::CatmullRomSpline;

    #[test]
    fn test_f64() {
        let points = [(0.0, 1.0), (0.5, 0.5), (1.0, 0.0)];
        let interpolator = CatmullRomSpline::try_new(&points).unwrap();
        let val = interpolator.try_value(0.75).unwrap();
        assert!((val - 0.27083333333333337_f64).abs() < f64::EPSILON);
    }

    #[cfg(feature = "decimal")]
    #[test]
    fn test_decimal() {
        let points = [
            (Decimal::new(0, 0), Decimal::new(1, 0)),
            (Decimal::new(5, 1), Decimal::new(5, 1)),
            (Decimal::new(1, 0), Decimal::new(0, 0)),
        ];
        let interpolator = CatmullRomSpline::try_new(&points).unwrap();
        let val = interpolator.try_value(Decimal::new(75, 2)).unwrap();
        assert_eq!(
            val,
            Decimal::from_str_exact("0.2708333333333333333333333333").unwrap()
        );
    }
}
