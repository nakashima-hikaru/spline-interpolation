use crate::{HermiteSplineError, InterpolationValue};
use nalgebra::{Matrix4, Vector4};
use num_traits::Zero;

struct Point3<V> {
    pub x: V,
    pub y: V,
    pub dydx: V,
}

pub struct HermiteSpline<V: InterpolationValue> {
    points: Vec<Point3<V>>,
    m: Matrix4<V>,
}

impl<V: InterpolationValue> HermiteSpline<V> {
    /// Creates a new instance of `HermiteSpline` from a slice of raw points.
    ///
    /// # Arguments
    ///
    /// * `raw_points` - A slice of tuples representing the raw points of the spline.
    ///                 Each tuple should contain three elements: the x-coordinate, the y-coordinate,
    ///                 and the derivative of y with respect to x (dy/dx).
    ///
    /// # Returns
    ///
    /// Returns a `Result` containing the constructed `HermiteSpline` on success,
    /// or a `HermiteSplineError` if the raw points are not in ascending order based on x-coordinate.
    ///
    /// # Example
    ///
    /// ```
    /// use spline_interpolation::hermite_spline::HermiteSpline;
    ///
    /// let raw_points = [(0.0, 0.0, 1.0), (1.0, 1.0, 2.0), (2.0, 0.0, -1.0)];
    ///
    /// let spline = HermiteSpline::try_new(&raw_points);
    /// assert!(spline.is_ok());
    /// ```
    pub fn try_new(raw_points: &[(V, V, V)]) -> Result<Self, HermiteSplineError<V>> {
        let mut temp = raw_points[0].0;

        let mut points = Vec::new();
        for &(x, y, dydx) in raw_points {
            let point = Point3 { x, y, dydx };
            if point.x < temp {
                return Err(HermiteSplineError::PointOrderError);
            }
            temp = point.x;
            points.push(point);
        }
        let m = Matrix4::new(
            V::from_i8(2).unwrap(),
            V::from_i8(-2).unwrap(),
            V::one(),
            V::one(),
            V::from_i8(-3).unwrap(),
            V::from_i8(3).unwrap(),
            V::from_i8(-2).unwrap(),
            -V::one(),
            V::zero(),
            V::zero(),
            V::one(),
            V::zero(),
            V::one(),
            V::zero(),
            V::zero(),
            V::zero(),
        );
        Ok(Self { points, m })
    }

    /// Tries to evaluate the interpolated value of Hermite spline at a given point x.
    ///
    /// # Arguments
    ///
    /// * `x` - The point at which to evaluate the Hermite spline.
    ///
    /// # Returns
    ///
    /// Returns `Ok(value)` if the evaluation was successful, where `value` is the evaluated value of the Hermite spline at `x`.
    /// Returns `Err(error)` if the evaluation was unsuccessful, where `error` is an enum variant indicating the reason for failure.
    ///
    /// # Errors
    ///
    /// Returns `OutOfLowerBound(x)` if `x` is less than the minimum x-coordinate value of any point in the Hermite spline.
    /// Returns `OutOfUpperBound(x)` if `x` is greater than the maximum x-coordinate value of any point in the Hermite spline.
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
                let d = Vector4::new(delta3, delta2, delta, V::from_i8(1).unwrap());
                let f = Vector4::new(point.y, next_point.y, point.dydx * h, next_point.dydx * h);
                Ok((d.transpose() * self.m * f).x)
            }
        }
    }
}
