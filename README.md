# **pledge** 📜

*An experimental playground for polynomial commitment schemes*

## Overview

**pledge** is my personal playground experimenting with **polynomial commitment schemes**.

This project is intended as a sandbox for trying out ideas, comparing different commitment schemes (e.g., KZG, IPA, FRI), and understanding their tradeoffs in terms of efficiency, security, and implementation complexity.

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

* [x] Support for KZG (Kate-Zaverucha-Goldberg)
    * [ ] Batch Opening
* [ ] Inner product argument (IPA) commitments
* [ ] FRI (Fast Reed-Solomon IOP of Proximity) commitments
* [ ] WASM bindings for browser-based experiments

## License

[MIT License](./LICENSE)