# Copyright 2022 Zheng Ma
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""rte_dataset dataset."""

import os
import pathlib

import numpy as np
import tensorflow as tf
import tensorflow_datasets as tfds

from rte_dataset.builders import features, pipeline, processings

os.environ["NO_GCE_CHECK"] = "true"
tfds.core.utils.gcs_utils._is_gcs_disabled = True

features.register_feature("psi_label", tf.float32, [features.NUM_PHASE_COORDS])
features.register_feature(
    "boundary_scattering_kernel",
    tf.float32,
    [features.NUM_BOUNDARY_COORDS, features.NUM_VELOCITY_COORDS],
)
FEATURES = features.FEATURES
PHASE_FEATURE_AXIS = {
    k: FEATURES[k][1].index(features.NUM_PHASE_COORDS) - len(FEATURES[k][1])
    for k in FEATURES
    if features.NUM_PHASE_COORDS in FEATURES[k][1]
}
BOUNDARY_FEATURE_AXIS = {
    k: FEATURES[k][1].index(features.NUM_BOUNDARY_COORDS) - len(FEATURES[k][1])
    for k in FEATURES
    if features.NUM_BOUNDARY_COORDS in FEATURES[k][1]
}


def _get_config_names(file):
    config_path = pathlib.Path(__file__).parent / file
    with open(config_path, "r") as f:
        return f.read().splitlines()


class Builder(tfds.core.GeneratorBasedBuilder):
    """DatasetBuilder for rte_dataset dataset."""

    VERSION = tfds.core.Version("0.0.2")
    RELEASE_NOTES = {
        "0.0.2": "Fix weights",
    }
    BUILDER_CONFIGS = [
        tfds.core.BuilderConfig(name=name) for name in _get_config_names("CONFIGS.txt")
    ]
    MANUAL_DOWNLOAD_INSTRUCTIONS = """
        Please download the raw dataset to project_root/data/raw_data
    """

    def _info(self) -> tfds.core.DatasetInfo:
        """Returns the dataset metadata."""
        feature_info_dict = {}
        for k, v in FEATURES.items():
            _, shape = v
            new_shape = [None if isinstance(k, str) else k for k in shape]
            feature_info_dict[k] = tfds.features.Tensor(
                shape=new_shape, dtype=np.float32, encoding="zlib"
            )
        self.metadata_dict = tfds.core.MetadataDict(
            {
                "phase_feature_axis": PHASE_FEATURE_AXIS,
                "boundary_feature_axis": BOUNDARY_FEATURE_AXIS,
            }
        )
        return self.dataset_info_from_configs(
            features=tfds.features.FeaturesDict(feature_info_dict),
            metadata=self.metadata_dict,
            homepage="https://github.com/mazhengcn/deeprte",
        )

    def _split_generators(self, dl_manager: tfds.download.DownloadManager):
        """Returns SplitGenerators."""

        datasets_dir = dl_manager.manual_dir / self.builder_config.name
        filenames = [
            file.name for file in datasets_dir.iterdir() if file.name.endswith(".mat")
        ]

        data_pipeline = pipeline.DataPipeline(datasets_dir, filenames)
        raw_data = data_pipeline.process(normalization=True)

        self.metadata_dict.update(
            {
                "normalization": tf.nest.map_structure(
                    lambda x: str(x), raw_data["normalization"]
                )
            }
        )

        return {
            "train": self._generate_examples(raw_data),
        }

    def _generate_examples(self, raw_data):
        """Yields examples."""

        num_examples = raw_data["shape"]["num_examples"]

        for i in range(num_examples):
            np_example = {
                **tf.nest.map_structure(lambda x: x[i], raw_data["functions"]),
                **raw_data["grid"],
            }
            tensor_dict = processings.np_to_tensor_dict(
                np_example, raw_data["shape"], FEATURES.keys()
            )
            yield i, tf.nest.map_structure(lambda x: x.numpy(), tensor_dict)
