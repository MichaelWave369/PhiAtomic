# ΦAtomic · Recursive Nuclear Topology

An open research prototype exploring whether recursive mathematical features add predictive information beyond conventional nuclear composition models. No claim that atoms are Mandelbrot sets, that the universe is recursive, or that an orbit describes real nuclear dynamics.

## Run

Open `index.html` in any modern browser. The lab is fully offline and has no build step or external assets.

For the benchmark, Python 3.10+:

```sh
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python phiatomic.py --output result.json
```

On Windows, `py` can substitute for `python`. The output includes per-fold MAE, pooled MAE, and the literal recursion escape tests. `data/demo.csv` contains 26 selected illustrative nuclides with rounded binding energies; it is not a validated snapshot of AME. `docs/demo-result.json` is the actual initial local run.

## Architecture

- `index.html`: standalone orbit and Mandelbrot parameter viewer.
- `phiatomic.py`: literal, complex and centered mappings, liquid-drop and shell baselines, grouped held-out ridge benchmarks.
- `data/demo.csv`: small rounded demonstration fixture.
- `tests/`: regression tests; GitHub Actions CI included.
- `docs/RESEARCH.md`: hypothesis, limitations and validation plan.

## Real data

Use evaluated binding energy per nucleon from [AME2020](https://doi.org/10.1088/1674-1137/abddaf) or [IAEA LiveChart](https://iaea-nds.github.io/lc_api_notebook/). Normalize a verified snapshot to columns `name,z,n,binding_mev_per_nucleon`. Convert keV/nucleon to MeV/nucleon by dividing by 1000. Record source URL, evaluation date, checksum, uncertainty and measured/extrapolated status separately. Do not mix isomers, duplicate ground states, or silently treat unknown values as zero. Run `python phiatomic.py --data path/to/normalized.csv` after checking the data. No full evaluated dataset has been imported in v0.1.

## License

MIT for original project code and prose; third-party source data remains subject to its original terms.
