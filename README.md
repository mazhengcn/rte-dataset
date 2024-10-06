# Radiative transfer (transport) Dataset

Dataset for training and evaluating deeprte model.


## MATLAB Datasets

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

## Build tfrecord (for tfds) or array_record (for grain) dataset for training

To build tfds or grain dataset, run the following command from root of the repo

```bash
bash build_dataset.sh ${RAW_DATA_DIR} ${TFDS_DIR} ${GRAIN_DIR}
```

which by default will read matlab data from `RAW_DATA_DIR=assets/raw_data` and generate tfds and grain dataset under `TFDS_DIR=assets/tfds` and `GRAIN_DIR=assets/grain`.

The dataset for training has following structure:

```
{
 'boundary': (8, 1920),
 'boundary_coords': (8, 1920, 4),
 'boundary_weights': (8, 1920),
 'phase_coords': (8, 128, 4),
 'position_coords': (8, 1600, 2),
 'psi_label': (8, 128),
 'scattering_kernel': (8, 128, 24),
 'self_scattering_kernel': (8, 24, 24),
 'sigma': (8, 1600, 2),
 'velocity_coords': (8, 24, 2),
 'velocity_weights': (8, 24)
}
```
