#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse
import numpy as np
import pickle
import pandas as pd
from sklearn.ensemble import ExtraTreesClassifier

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Predict modification sites using a pre-trained model.")
    parser.add_argument(
        "--features_file",
        type=str,
        required=True,
        help="Path to the input features CSV file."
    )
    parser.add_argument(
        "--model_file",
        type=str,
        required=True,
        help="Path to the pre-trained model file (.pkl)."
    )
    parser.add_argument(
        "--output_file",
        type=str,
        required=True,
        help="output_file to save the prediction results."
    )
    parser.add_argument(
        "--mod",
        type=str,
        required=True,
        help="Nead predict which modification."
    )    
    args = parser.parse_args()
    
    features_file = args.features_file
    model_file = args.model_file
    output_file = args.output_file
    mod = args.mod

    # Validate input files
    if not os.path.isfile(features_file):
        raise FileNotFoundError(f"Features file not found: {features_file}")
    if not os.path.isfile(model_file):
        raise FileNotFoundError(f"Model file not found: {model_file}")
    
    # Load the pre-trained model
    with open(model_file, "rb") as f:
        model = pickle.load(f)
    
    # Read the features file
    df = pd.read_csv(
        features_file,
        header=None,
        names=[
            "transcript_ID", "position", "base_type", "coverage", "ins", "ins_len",
        "del", "del_len", "del_site", "mis", "mis_A", "mis_C", "mis_G",
        "mis_T", "base_qual_mean", "base_qual_STD", "base_qual_count_0"
        ]
    )
    
    # Filter rows where base_type is "T"
    df = df[df["base_type"] == mod].copy()
    
    # Select columns for site information
    site_info = df[["transcript_ID", "position", "base_type", "coverage"]].reset_index(drop=True)
    
    # Drop unnecessary columns
    df.drop(["transcript_ID", "position", "base_type", "coverage", "mis_T"], axis=1, inplace=True)
    
    # Ensure there are no missing values and data types are correct
    df = df.fillna(0)
    numeric_columns = ["ins", "ins_len", "del", "del_len", "del_site", "mis",
        "mis_A", "mis_C", "mis_G", "base_qual_mean",
        "base_qual_STD", "base_qual_count_0"]
    df[numeric_columns] = df[numeric_columns].apply(pd.to_numeric, errors='coerce').fillna(0)
    
    # Predict probabilities
    try:
        ET1_predict_proba = model.predict_proba(df)[:, 1]
    except Exception as e:
        raise RuntimeError(f"Error during prediction: {e}")
    
    # Create prediction dataframe
    pred_df = pd.DataFrame(ET1_predict_proba, columns=["modi"])
    
    # Combine site information with predictions
    site_info = site_info.reset_index(drop=True).join(pred_df)
    
    site_info.to_csv(output_file, index=False, header=False)

if __name__ == "__main__":
    main()
