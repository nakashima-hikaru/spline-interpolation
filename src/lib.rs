use num_traits::{FromPrimitive, Num};
use std::fmt::Debug;
use std::ops::{AddAssign, MulAssign, Neg};
#[cfg(feature = "decimal")]
use rust_decimal::Decimal;
use thiserror::Error;

pub mod catmull_rom_spline;
pub mod hermite_spline;

pub trait InterpolationValue:
'static
+ Num
+ Copy
+ PartialOrd
+ Neg<Output=Self>
+ Debug
+ Sized
+ AddAssign
+ MulAssign
+ FromPrimitive
{}

impl InterpolationValue for f32 {}

impl InterpolationValue for f64 {}

#[cfg(feature = "decimal")]
impl InterpolationValue for Decimal {}

#[derive(Error, Debug)]
pub enum HermiteSplineError<V: InterpolationValue> {
    #[error("points must be sorted")]
    PointOrderError,
    #[error("out of lower bound: {0}")]
    OutOfLowerBound(V),
    #[error("out of upper bound: {0}")]
    OutOfUpperBound(V),
    #[error("length of inputs: {0} is not enough points for construction")]
    InsufficientPointsError(usize),
}
