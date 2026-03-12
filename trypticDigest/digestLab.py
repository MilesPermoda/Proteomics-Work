# imports
import csv
import pandas as pd
from Bio import SeqIO
from pyteomics import parser
import seaborn as sns
import matplotlib.pyplot as plt

#global variables
pep_df = pd.DataFrame(columns=["Peptide"])
unique_df = pd.DataFrame(columns=["Peptide"])
rep_df = pd.DataFrame(columns=["Peptide"])

# parses the proteome(fasta input) and digests all the proteins into individual peptides
# and adds these peptides to a dataframe
def parse_prot():
    for record in SeqIO.parse("uniprotMethanococcus1-15-26.fasta", "fasta"):
        sequence = str(record.seq)
        # use pyteomics to generate tryptic peptides
        tryptic_peptides = parser.cleave(sequence, "trypsin")
        for pep in tryptic_peptides:
            temp_df = pd.DataFrame([{
                "Peptide": pep
            }])
            global pep_df
            pep_df = pd.concat([pep_df, temp_df], ignore_index=True)

# takes the full pep dataframe and counts the number of unique and non-unique peptides
# and uses them to calculate the percentage.
def percent_unique():
    total_unique_peptides = pep_df['Peptide'].nunique()
    non_unique_peptides = (
            pep_df['Peptide'].value_counts() > 1
    ).sum()

    percentage = 100 * non_unique_peptides / total_unique_peptides
    print("percentage:", percentage)

# creates a new dataframe with any repeat peptides
# and appends a unique vs non unique classifier to both it
# and the unique dataframe, then calculates the length of their peptides before
# concatenating both into a combined dataframe to be used for graphing and such.
# finally, this combined dataframe is exported as a csv file.
def dataframes():
    global rep_df
    global unique_df
    unique_df = pep_df[~pep_df['Peptide'].duplicated(keep=False)]
    rep_df = pep_df[pep_df['Peptide'].duplicated(keep=False)].drop_duplicates()
    unique_df['Length'] = unique_df['Peptide'].str.len()
    unique_df['Group'] = 'Unique'
    rep_df['Length'] = rep_df['Peptide'].str.len()
    rep_df['Group'] = 'Non_Unique'
    together = pd.concat([unique_df, rep_df], ignore_index=True)

    together.to_csv(
        "peptides.csv",
        index=False,
        sep=",",
        encoding="utf-8",
        na_rep="NA"
    )

    grapher(together)

# Plots a violin plot and box plot of the unique vs non-unique variables
# and their length distributions, as well as providing a second plot
# with a zoomed in view for ease of understanding/visualizaion
def grapher(df):
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    # create violinplots
    sns.violinplot(x="Group", y="Length", data=df, inner="quartile", ax=ax[0])
    ax[0].set_title("Violin Plot of Lengths of Unique and Non-Unique Peptides")
    ax[0].set_xlabel("Peptide Groupings")
    ax[0].set_ylabel("Length of Peptides")

    sns.violinplot(x="Group", y="Length", data=df, inner="quartile", ax=ax[1])
    ax[1].set_title("Zoomed in View (Cutoff at y=20)")
    ax[1].set_xlabel("Peptide Groupings")
    ax[1].set_ylabel("Length of Peptides")
    ax[1].set_ylim(0, 20)
    plt.show()

    #create boxplots 
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    sns.boxplot(data=df, x='Group', y='Length', ax=ax[0])
    ax[0].set_title("Boxplots Showcasing Length Variance of Unique and Non-Unique Peptides")
    ax[0].set_xlabel("Peptide Groupings")
    ax[0].set_ylabel("Length of Peptides")

    sns.boxplot(data=df, x='Group', y='Length', ax=ax[1])
    ax[1].set_title("Zoomed in View (Cutoff at y=20)")
    ax[1].set_xlabel("Peptide Groupings")
    ax[1].set_ylabel("Length of Peptides")
    plt.ylim(0, 20)
    plt.show()

if __name__ == "__main__":
    parse_prot()
    percent_unique()
    dataframes()
    file = "pepML.py"
    # runs the machine learning script!
    exec(open(file).read())


