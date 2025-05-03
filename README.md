# Radiative transfer (transport) Dataset

Dataset for training and evaluating deeprte model.

[**Raw Datasets**](#raw-datasets) | [**Build Training Datasets**](#build-training-datasets)

## Raw Datasets

The original or "raw" data are generated using a deterministic solver written in Matlab (this may change in future). All the data are stored under `assets/raw_data` folder in the following nameing scheme: `g*-sigma_a*-sigma_t*`, e.g.,

```bash
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
| `list_Psi`     | `[2M, I, J, N]` | numerical solutions as training labels               |
| `list_psiL`    | `[M, J, N]`     | left boundary values                                 |
| `list_psiR`    | `[M, J, N]`     | right boundary values                                |
| `list_psiB`    | `[M, I, N]`     | bottom boundary values                               |
| `list_psiT`    | `[M, I, N]`     | top boundary values                                  |
| `list_sigma_a` | `[I, J, N]`     | absorption coefficient function                      |
| `list_sigma_T` | `[I, J, N]`     | total coefficient function                           |
| `ct` and `st`  | `[1, M]`        | discrete coordinates (quadratures) in velocity space |
| `omega`        | `[1, M]`        | weights of velocity coordinates                      |

## Build Training Datasets

To build a tfds or grain dataset, run the following command from the root of the repository:

```bash
bash build_dataset.sh ${RAW_DATA_DIR} ${TFDS_DIR} ${GRAIN_DIR}
```

By default, this script reads MATLAB data from `RAW_DATA_DIR=assets/raw_data` and generates tfds and grain datasets under `TFDS_DIR=assets/tfds` and `GRAIN_DIR=assets/grain`.

Each dataset element for training has the following structure:

```bash
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

```bash
NUM_POSITION_COORDS = number of position coordinates
NUM_VELOCITY_COORDS = number of velocity coordinates
NUM_PHASE_COORDS = number of phase coordinates
NUM_BOUNDARY_COORDS = number of boundary coordinates
```
