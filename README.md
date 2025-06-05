# RTE Datasets

Datasets for training and evaluating DeepRTE model, see [repo](https://github.com/mazhengcn/deeprte).

[**Raw Datasets**](#raw-datasets) | [**Build Training Datasets**](#build-training-datasets)

## Raw Datasets

The original or "raw" data are generated using a deterministic solver written in Matlab and Python. All the data are stored in the following nameing scheme: `g*-sigma_a*-sigma_t*`, e.g.,

```text
train
├── g0.1-sigma_a3-sigma_t6
│   ├── config.m
│   └── g0.1-sigma_a3-sigma_t6.mat
├── g0.5-sigma_a3-sigma_t6
│   ├── config.m
│   └── g0.5-sigma_a3-sigma_t6.mat
└── g0.8-sigma_a3-sigma_t6
    ├── config.m
    └── g0.8-sigma_a3-sigma_t6.mat
```

Each `*.mat` files has the following keys and array shapes:

| Key            | Array Shape     | Description                                          |
| -------------- | --------------- | ---------------------------------------------------- |
| `list_Psi`     | `[2M, I, J, N]` | Solutions of RTE (labels)                            |
| `list_psiL`    | `[M, J, N]`     | Left boundary values                                 |
| `list_psiR`    | `[M, J, N]`     | Right boundary values                                |
| `list_psiB`    | `[M, I, N]`     | Bottom boundary values                               |
| `list_psiT`    | `[M, I, N]`     | Top boundary values                                  |
| `list_sigma_a` | `[I, J, N]`     | Absorption cross section                             |
| `list_sigma_T` | `[I, J, N]`     | Total cross section                                  |
| `ct`, `st`  | `[1, M]`        | Angular quadrature points                            |
| `omega`        | `[1, M]`        | Angular quadrature weights                           |

where `N` is number of samples, `M` is number of angular quadrature points, `I` and `J` is the number of mesh points in $x$ and $y$ direction respectively.

Inference (test) and pretraining datasets are hosted on Huggingface: https://huggingface.co/datasets/mazhengcn/rte-dataset. Download datasets to `DATA_DIR` with (ensure `huggingface-cli` is installed; if you followed the setup above, it is already included):

```bash
huggingface-cli download mazhengcn/rte-dataset \
    --exclude=interim/* \
    --repo-type=dataset \
    --local-dir=${DATA_DIR}
```

The resulting folder structure should be (for inference, only datasets under `raw/test` are needed):

```text
${DATA_DIR}
├── processed
│   └── tfds      # Processed TFDS dataset for pretraining.
├── raw
│   ├── test      # Raw MATLAB dataset for test/inference.
│   └── train     # Raw MATLAB dataset for pretraining using grain.
└── README.md
```

## Build Training Datasets

To build a tfds or grain dataset, run the following command from the root of the repository:

```bash
bash build_dataset.sh ${RAW_DATA_DIR} ${TFDS_DIR} ${GRAIN_DIR}
```

By default, this script reads MATLAB data from `RAW_DATA_DIR=assets/raw_data` and generates tfds and grain datasets under `TFDS_DIR=assets/tfds` and `GRAIN_DIR=assets/grain`.

Each dataset element for training has the following structure:

```python
FEATURES = {
    "phase_coords": (NUM_PHASE_COORDS, 2 * NUM_DIM),
    "boundary_coords": (NUM_BOUNDARY_COORDS, 2 * NUM_DIM),
    "boundary_weights": (NUM_BOUNDARY_COORDS),
    "position_coords": (NUM_POSITION_COORDS, NUM_DIM),
    "velocity_coords": (NUM_VELOCITY_COORDS, NUM_DIM),
    "velocity_weights": (NUM_VELOCITY_COORDS),
    "boundary": (NUM_BOUNDARY_COORDS),
    "sigma": (NUM_POSITION_COORDS, 2),
    "scattering_kernel": (NUM_PHASE_COORDS, NUM_VELOCITY_COORDS),
    "self_scattering_kernel": (NUM_VELOCITY_COORDS, NUM_VELOCITY_COORDS),
}
```

where

```python
NUM_POSITION_COORDS = number of position coordinates
NUM_VELOCITY_COORDS = number of velocity coordinates
NUM_PHASE_COORDS = number of phase coordinates
NUM_BOUNDARY_COORDS = number of boundary coordinates
```
