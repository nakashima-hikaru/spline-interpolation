use nalgebra::{Matrix4, Vector4};
use num_traits::Zero;


struct Point<V> {
    pub(crate) x: V,
    pub(crate) y: V,
    pub(crate) dydx: V,
}

struct HermiteSpline {
    points: Vec<Point<f64>>,
}

impl HermiteSpline {
    const M: Matrix4<f64> = Matrix4::new(2., -2., 1., 1., -3., 3., -2., -1., 0., 0., 1., 0., 1., 0., 0., 0.);
    fn new(raw_points: &[(f64, f64, f64)]) -> Self {
        let mut temp = raw_points[0].0;

        let mut points = Vec::new();
        for raw_point in raw_points {
            let point = Point {
                x: raw_point.0,
                y: raw_point.1,
                dydx: raw_point.2,
            };
            if point.x < temp {
                panic!("points must be sorted");
            }
            temp = point.x;
            points.push(point);
        }
        Self { points }
    }

    fn value(self, x: f64) -> f64 {
        let pos = self.points.binary_search_by(|point| point.x.total_cmp(&x));
        let pos = if pos.is_ok() {
            pos.unwrap()
        } else {
            let ret = pos.unwrap_err();
            if ret.is_zero() {
                panic!("need extrapolation");
            }
            if ret > self.points.len() {
                panic!("need extrapolation");
            }
            ret - 1
        };
        let point = &self.points[pos];
        let next_point = &self.points[pos + 1];
        let h = next_point.x - point.x;
        let delta = (x - point.x) / h;
        let delta2 = delta * delta;
        let delta3 = delta2 * delta;
        let d = Vector4::new(delta3, delta2, delta, 1.);
        let f = Vector4::new(point.y, next_point.y, point.dydx * h, next_point.dydx * h);
        (d.transpose() * Self::M * f).x
    }
}
