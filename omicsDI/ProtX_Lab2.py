# imports
import requests
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

count = 0
base_url = None
species_name = None
params = None
datasets = None
def setup(): # sets up the program for the 3 species
    global count
    global base_url
    base_url = "https://www.omicsdi.org/ws/dataset/search"
    if (count == 0):
        global species_name
        species_name = "Methanococcus maripaludis"

        #filter by 'omics_type:Proteomics' to exclude non-proteomics data
        global params
        params = {
            "query": f'"{species_name}" AND omics_type:"Proteomics"',
            "size": 100,  # Number of results to return --> increase this if you expect more
            "start": 0,
        }
    elif (count == 1):
        species_name = "Shewanella oneidensis"


        params = {
            "query": f'"{species_name}" AND omics_type:"Proteomics"',
            "size": 100,  # Number of results to return --> increase this if you expect more
            "start": 0,
        }
    elif (count == 2):
        species_name = "Bacillus amyloliquefaciens"


        params = {
            "query": f'"{species_name}" AND omics_type:"Proteomics"',
            "size": 100,  # Number of results to return --> increase this if you expect more
            "start": 0,
        }
    else:
        print("Out of bounds error")
        exit()
    count = count + 1

def getter(): # gets the correct datasets
    global datasets
    datasets = None

    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        data = response.json()


        datasets = data.get("datasets", [])

        if not datasets:
            print(f"No datasets found for {species_name}.")

        print(f"Found {data['count']} datasets across all repositories.\n")

    except requests.exceptions.RequestException as e:
        print(f"Error connecting to OmicsDI: {e}")
    parser(datasets)
    # pprint.pprint(datasets)
def parser(data): # parses through the datasets and creates dataframes with relevant information
    count_df = pd.DataFrame(columns=['Accession','Source','Species Group', 'Year'])
    for ds in data:
        # OmicsDI returns the source DB and the ID
        accession = ds["id"]
        database = ds["source"]
        title = ds["title"]
        date = ds.get("publicationDate")
        if date:
            year = int(date[:4])
        else:
            year = None

        temp_df = pd.DataFrame([{
            "Accession": accession,
            "Source": database,
            "Species Group": count,
            "Year": year
        }])
        count_df = pd.concat([count_df, temp_df], ignore_index=True)

    # calls the counting functions on the dataframe
    counter(count_df)
    yearCounter(count_df)
def counter(df): # counts the occurrences of the databases and graphs them
    print("Source counts for species  " + species_name)
    valCounts = df['Source'].value_counts()
    valCounts = valCounts.reset_index()
    valCounts.columns = ['Source', 'Count']
    #plot time!!
    plt.figure(figsize=(8, 6))
    sns.barplot(data=valCounts, x='Source', y='Count')
    plt.title("Number of Datasets per Database for " + species_name)
    plt.xlabel('Database')
    plt.ylabel('Number of Datasets')
    plt.show()
def yearCounter(df): # graphs the occurrences of the years and graphs them
    print("Year counts for species  " + species_name)
    yearCounts = df['Year'].value_counts()
    yearCounts = yearCounts.reset_index()
    yearCounts.columns = ['Year', 'Count']
    #plot
    plt.figure(figsize=(8, 6))
    sns.barplot(data=yearCounts, x='Year', y='Count')
    plt.title("Datasets per Year for " + species_name)
    plt.xlabel('Year')
    plt.ylabel('Number of Datasets')
    plt.show()
if __name__ == "__main__":
    # iteration 1
    setup()
    getter()
    # iteration 2
    setup()
    getter()
    # iteration 3
    setup()
    getter()
