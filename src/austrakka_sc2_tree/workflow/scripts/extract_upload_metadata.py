import pandas as pd

df = pd.read_csv(snakemake.input.nextclade_tsv, sep="\t")

columns = {}

for col in snakemake.config['extract']:
    try:
        new_col, operation = col.split(":")
    except ValueError:
        new_col = col
        operation = f"row['{col}']"
    columns[new_col] = operation
# New DataFrame to store the results
upload_df = pd.DataFrame()


# Function to apply transformations based on config
def apply_transformations(row):
    transformed_row = {}
    for new_col, operation in columns.items():
        transformed_row[new_col] = eval(operation)
    return pd.Series(transformed_row)

# Applying the transformations to the DataFrame
upload_df = df.apply(apply_transformations, axis=1)

upload_df.to_csv(snakemake.output.metadata_csv, index=False)
