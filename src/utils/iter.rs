#[cfg(feature = "parallel")]
use rayon::prelude::*;

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

// Define an extension trait for seamless integration
#[cfg(feature = "parallel")]
pub trait UnifiedFold: Sized {
    fn unified_fold<F>(self, init: Self::Item, f: F) -> Self::Item
    where
        F: Fn(Self::Item, Self::Item) -> Self::Item + Send + Sync + Clone,
        Self: ParallelIterator, // Base on std Iterator
        Self::Item: Send;
}

#[cfg(feature = "parallel")]
impl<I> UnifiedFold for I
where
    I: ParallelIterator, // Iterator itself needs Send for par_bridge
    I::Item: Clone + Send + Sync,
{
    fn unified_fold<F>(self, init: I::Item, f: F) -> I::Item
    where
        F: Fn(I::Item, I::Item) -> I::Item + Send + Sync + Clone,
        I::Item: Send,
    {
        use rayon::iter::ParallelIterator;

        self.fold_with(init.clone(), f.clone())
            .reduce(|| init.clone(), f)
    }
}

#[cfg(not(feature = "parallel"))]
pub trait UnifiedFold: Sized {
    fn unified_fold<T, F>(self, init: T, f: F) -> T
    where
        F: Fn(T, Self::Item) -> T + Send + Sync + Clone,
        Self: Iterator, // Base on std Iterator
        Self::Item: Send;
}

#[cfg(not(feature = "parallel"))]
impl<I> UnifiedFold for I
where
    I: Iterator,
    I::Item: Send,
{
    fn unified_fold<T, F>(self, init: T, f: F) -> T
    where
        F: Fn(T, I::Item) -> T + Send + Sync + Clone,
        Self: Iterator, // Base on std Iterator
        I::Item: Send,
    {
        self.fold(init, f)
    }
}
