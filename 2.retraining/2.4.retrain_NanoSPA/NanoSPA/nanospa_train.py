import argparse
import os
import pandas as pd
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve,auc,precision_recall_curve
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.preprocessing import StandardScaler
import pandas as pd

def parse_arguments():
    parser = argparse.ArgumentParser(description="Train a model and save the results")
    parser.add_argument(
        "--features_file",
        type=str,
        required=True,
        help="Path to the input features CSV file"
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Directory to save the output results"
    )
    return parser.parse_args()

def main():
    # Parse command-line arguments
    args = parse_arguments()
    
    features_file = args.features_file
    #bed_file = args.bed_file
    output_dir = args.output_dir

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    df = pd.read_csv(
        features_file,
        header=None,
        names=[
        "transcript_ID", "position", "base_type", "coverage", "ins", "ins_len",
        "del", "del_len", "del_site", "mis", "mis_A", "mis_C", "mis_G",
        "mis_T", "base_qual_mean", "base_qual_STD", "base_qual_count_0","label"
        ]
    )
    
    df_label_1 = df[df['label'] == 1]
    df_label_0_sample = df[df['label'] == 0].sample(n=5000, random_state=42)
    df_new = pd.concat([df_label_1, df_label_0_sample])
    
    df_new.reset_index(drop=True, inplace=True)

    # Features and labels
    X = df_new[[
        "ins", "ins_len", "del", "del_len", "del_site", "mis",
        "mis_A", "mis_C", "mis_G", "base_qual_mean",
        "base_qual_STD", "base_qual_count_0"
    ]]
    y = df_new["label"]

    # Split the dataset
    X_train, X_temp, y_train, y_temp = train_test_split(
        X, y, test_size=0.4, random_state=42, stratify=y
    )
    X_val, X_test, y_val, y_test = train_test_split(
        X_temp, y_temp, test_size=0.5, random_state=42, stratify=y_temp
    )

    # Train the model
    ext_model = ExtraTreesClassifier(
        n_estimators=200, criterion="gini", max_depth=None,
        min_samples_split=2, random_state=42
    )
    ext_model.fit(X_train, y_train)

    # Validation set prediction
    y_val_pred_proba = ext_model.predict_proba(X_val)[:, 1]
    val_auc = roc_auc_score(y_val, y_val_pred_proba)
    print(f"Validation AUC: {val_auc:.4f}")

    # Test set prediction
    y_test_pred_proba = ext_model.predict_proba(X_test)[:, 1]
    test_auc = roc_auc_score(y_test, y_test_pred_proba)
    print(f"Test AUC: {test_auc:.4f}")

    # Plot ROC Curve
    fpr, tpr, _ = roc_curve(y_test, y_test_pred_proba)
    plt.figure()
    plt.plot(fpr, tpr, label=f"Test AUC = {test_auc:.4f}")
    plt.plot([0, 1], [0, 1], color='gray', linestyle='--')
    plt.xlabel("False Positive Rate (1 - Specificity)")
    plt.ylabel("True Positive Rate (Sensitivity)")
    plt.title("ROC Curve")
    plt.legend(loc="lower right")
    roc_path = os.path.join(output_dir, 'test_ROCAUC_train.pdf')
    plt.savefig(roc_path, bbox_inches='tight', dpi=300)
    plt.close()
    
    precision, recall, thresholds = precision_recall_curve(y_test, y_test_pred_proba)
    pr_auc = auc(recall, precision)
    print(f"Test PRAUC: {pr_auc:.4f}")
    plt.figure()
    plt.plot(recall, precision, label=f"Test AUC = {pr_auc:.4f}")
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('PRAUC Curve')
    plt.legend(loc="lower left")
    roc_path = os.path.join(output_dir, 'test_PRAUC_train.pdf')
    plt.savefig(roc_path, bbox_inches='tight', dpi=300)
    plt.close()

    # Save the model
    model_path = os.path.join(output_dir, "model.pkl")
    with open(model_path, "wb") as model_file:
        pickle.dump(ext_model, model_file)

if __name__ == "__main__":
    main()