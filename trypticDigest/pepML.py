from deeplc import DeepLC
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

orig_df = pd.DataFrame()

def ml():
    global orig_df
    df = pd.read_csv("peptides.csv")
    orig_df = df.copy()
    # start to transform df into something DeepLC can use as input
    df = df[["Peptide"]].rename(columns={"Peptide": "seq"})
    valid = (
        # DeepLC can't interpret non-standard amino acids like U
        ~df["seq"].str.contains("U") & ~df["seq"].str.contains("B") &
        # DeepLC doesnt work well with peptides larger than 60 AA's
        df["seq"].str.len().between(4, 60)

    )
    df = df[valid]
    df = df.dropna(subset=["seq"])
    df["seq"] = df["seq"].astype(str)
    df["seq"] = df["seq"].astype(str).str.strip()

    df["modifications"] = pd.Series("", index=df.index, dtype="string")

    # make predictions
    dlc = DeepLC()
    preds_uncal = dlc.make_preds(seq_df=df, calibrate=False)

    # append predictions to the dataframe used to make them, then join it to the original unaltered dataframe
    df["pred"] = preds_uncal
    orig_df = orig_df.join(df["pred"])

    # plotting time!!
    plt.figure(figsize=(12, 4))
    sns.scatterplot(
        x=np.arange(len(preds_uncal)),
        y=preds_uncal,
        s=5,
        alpha=0.5
    )
    plt.title("Predicted Retention Time for Individual Peptides")
    plt.xlabel("Individual Peptide")
    plt.ylabel("Predicted Retention Time")
    plt.xticks([])
    plt.tight_layout()
    plt.savefig("predicted_rt.png")
    plt.show()

    plt.figure(figsize=(12, 4))
    sns.scatterplot(
        x=orig_df["Length"],
        y=orig_df["pred"],
        s=75,
        alpha=0.5
    )
    plt.title("Predicted Retention Time Plotted Against Peptide Lengths")
    plt.xlabel("Length of Peptide")
    plt.ylabel("Predicted Retention Time")
    plt.savefig("compare.png")
    plt.show()

if __name__ == "__main__":
    ml()