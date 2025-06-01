/// Creates an iterator over refs or owned values, using parallel iterators if `parallel` feature is enabled.
///
/// Usage:
/// - `iter!(expr, ref)` for reference iterators (`iter()` or `par_iter()`).
/// - `iter!(expr, owned)` for owned iterators (`into_iter()` or `par_into_iter()`).
#[macro_export]
macro_rules! iter {
    ($e:expr, ref) => {{
        #[cfg(feature = "parallel")]
        {
            use rayon::prelude::*;
            let result = $e.par_iter();
            result
        }
        #[cfg(not(feature = "parallel"))]
        {
            let result = $e.iter();
            result
        }
    }};
    ($e:expr, owned) => {{
        #[cfg(feature = "parallel")]
        {
            use rayon::prelude::*;
            let result = $e.into_par_iter();
            result
        }
        #[cfg(not(feature = "parallel"))]
        {
            let result = $e.into_iter();
            result
        }
    }};
}
