#!/usr/bin/env python
# coding: utf-8
from pathlib import Path

from autorad.utils.preprocessing import get_paths_with_separate_folder_per_case
from autorad.data.dataset import ImageDataset
from autorad.feature_extraction.extractor import FeatureExtractor
import logging

def feature_extraction(input_folder, output_file):
    """
    Extracts features from images in the input folder and saves the results to the output file.

    Args:
        input_folder (str): Path to the folder containing the input images.
        output_file (str): Path to the output file where the extracted features will be saved.

    Returns:
        None
    """
    paths_df = get_paths_with_separate_folder_per_case(f"{input_folder}", relative=True)
    image_dataset = ImageDataset(
       paths_df,
       ID_colname="ID",
       root_dir=f"{input_folder}",
    )
    extractor = FeatureExtractor(image_dataset, extraction_params="MR_default.yaml")
    feature_df = extractor.run()
    feature_df.ID = feature_df.ID.astype(int)
    feature_df = feature_df.sort_values("ID")
    feature_df.ID = feature_df.ID.astype(str)
    df = feature_df
    df.to_csv(f"{output_file}", index = False)

def test_segmentation(input_folder):
    paths_df = get_paths_with_separate_folder_per_case(f"{input_folder}", relative=True)
    logging.getLogger().setLevel(logging.CRITICAL)
    image_dataset = ImageDataset(
       paths_df,
       ID_colname="ID",
       root_dir=f"{input_folder}",
    )
    image_dataset.plot_examples(n=10, window=None)

feature_extraction(Path(snakemake.input[0]).parent.parent, snakemake.output[0])