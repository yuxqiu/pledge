# **pledge** 📜

*An experimental playground for core components of ZKP*

## Overview

**pledge** is my personal playground experimenting with ZKP. As for now, it includes **polynomial commitment schemes** and **sumcheck**.

This project is intended as a sandbox for trying out ideas, comparing different commitment schemes, and understanding their tradeoffs in terms of efficiency, security, and implementation complexity.

## Getting Started

Clone the repository:

```bash
git clone https://github.com/yuxqiu/pledge.git
cd pledge
```

Run tests:

```bash
cargo test
```

## Roadmap

* PCS
    * [x] Support for KZG (Kate-Zaverucha-Goldberg)
        * [ ] Batch Opening
        * [ ] Zero Knowledge
    * [ ] Inner product argument (IPA) commitments
    * [ ] FRI (Fast Reed-Solomon IOP of Proximity) commitments
    * [ ] ECC based PCS
        * [x] Reed-Solomon Code
        * [ ] Linear Time Encodable Code
        * [ ] PCS
* Sumcheck
    * [x] CTY11
    * [x] VSBW13
    * [ ] [Blendy](https://github.com/compsec-epfl/efficient-sumcheck)
        * [x] Sumcheck for Multilinear Polynomial
        * [ ] Sumcheck for product of Multilinear Polynomial
* Others
    * [ ] WASM bindings for browser-based experiments

## License

[MIT License](./LICENSE)