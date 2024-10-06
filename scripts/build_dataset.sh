#!/bin/bash
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
set -e

RAW_DATA_DIR=${1:-"assets/raw_data"}
TFDS_DIR=${2:-"assets/tfds"}
GRAIN_DIR=${3:-"assets/grain"}

find "${RAW_DATA_DIR}/train" -mindepth 1 -maxdepth 1 -type d -exec basename {} \; > rte_dataset/builders/tfds/rte/CONFIGS.txt

TFDS_ARGS="--data_dir=${TFDS_DIR} --manual_dir=${RAW_DATA_DIR}/train"

# Build tfrecord dataset
tfds build rte_dataset/builders/tfds/rte ${TFDS_ARGS}
# Convert to also obtain array_record dataset for grain
tfds convert_format --root_data_dir=assets/tfds --out_file_format="array_record" --out_dir="${GRAIN_DIR}"
