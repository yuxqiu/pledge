use ark_poly::{DenseUVPolynomial, Polynomial, univariate::DensePolynomial};
use ark_poly_commit::kzg10;
use ark_std::rand::{SeedableRng, rngs::StdRng};
use ark_test_curves::bls12_381::{Bls12_381, Fr};
use criterion::{Criterion, criterion_group, criterion_main};
use poly_commit::kzg::KZG;

fn kzg_bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("KZG");

    let mut rng = StdRng::seed_from_u64(42);
    let poly: DensePolynomial<Fr> = DensePolynomial::rand(1000, &mut rng);

    // Not a fair comparison for now
    group.bench_function("KZG setup", |b| {
        b.iter(|| {
            KZG::<DensePolynomial<_>, Bls12_381>::setup(
                poly.degree().try_into().unwrap(),
                &mut rng,
            );
        });
    });
    group.bench_function("KZG setup (arkworks)", |b| {
        b.iter(|| {
            kzg10::KZG10::<Bls12_381, DensePolynomial<_>>::setup(poly.degree(), true, &mut rng)
        });
    });
    group.finish();
}

criterion_group!(benches, kzg_bench);
criterion_main!(benches);
