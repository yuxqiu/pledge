use ark_std::rand::{SeedableRng, rngs::StdRng};
use ark_test_curves::bls12_381::Fr;
use criterion::{Criterion, criterion_group, criterion_main};

use ark_ec::{CurveGroup, VariableBaseMSM};
use ark_ff::UniformRand;
use ark_test_curves::bls12_381::G1Projective;
use pledge::utils::msm::pippenger;

type F = Fr;

fn pippenger_bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("pippenger");

    // our impl is slower than `arkworks` because of `msm_bigint_wnaf`
    // (which is used when negation of group element is fast)
    //
    // See Signed Bucket Indexes for Multi-Scalar Multiplication (MSM)
    // - https://hackmd.io/@drouyang/signed-bucket-index
    // - Note, signed bucket index is actually different from WNAF, which is commonly
    //   used to speed up fixed MSM (in which the scalar mul can be speed up by using
    //   a pre-built lookup table). WNAF reduces the Hamming Weight of the numbers.
    //   - https://hackmd.io/@drouyang/S1GgWQvoj

    // `arkworks` impl of `msm_bigint_wnaf` is not very optimal
    // as it still allocates a bucket of size 2^c rather than 2^(c-1)
    //
    // Fixing it to use correct size requires negating the scalar and the corresponding point
    // (Yrrid's solution).

    let num: usize = 1024;
    let mut rng = StdRng::seed_from_u64(42);
    let bases: Vec<_> = (0..num)
        .map(|_| G1Projective::rand(&mut rng).into_affine())
        .collect();
    let scalars: Vec<F> = (0..num).map(|_| F::rand(&mut rng)).collect();
    group.bench_function(format!("pippenger MSM arkworks ({})", num), |b| {
        b.iter(|| {
            let _ = G1Projective::msm(&bases, &scalars).unwrap().into_affine();
        });
    });
    group.bench_function(format!("pippenger MSM ({})", num), |b| {
        b.iter(|| {
            let _ = pippenger::<G1Projective>(&bases, &scalars)
                .unwrap()
                .into_affine();
        });
    });

    let num: usize = 8192;
    let mut rng = StdRng::seed_from_u64(42);
    let bases: Vec<_> = (0..num)
        .map(|_| G1Projective::rand(&mut rng).into_affine())
        .collect();
    let scalars: Vec<F> = (0..num).map(|_| F::rand(&mut rng)).collect();
    group.bench_function(format!("pippenger MSM arkworks ({})", num), |b| {
        b.iter(|| {
            let _ = G1Projective::msm(&bases, &scalars).unwrap().into_affine();
        });
    });
    group.bench_function(format!("pippenger MSM ({})", num), |b| {
        b.iter(|| {
            let _ = pippenger::<G1Projective>(&bases, &scalars)
                .unwrap()
                .into_affine();
        });
    });

    group.finish();
}

criterion_group!(benches, pippenger_bench);
criterion_main!(benches);
