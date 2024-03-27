import numpy as np
import similaritymeasures
import pandas as pd


def calculate_mae(input_path, output_path):
    df = pd.read_csv(input_path)
    files = list(df.file.unique())
    mae_values = []
    for count1, file1 in enumerate(files):
        for count2, file2 in enumerate(files):
            if count1 < count2:
                subset1 = df[df.file == file1]
                subset2 = df[df.file == file2]
                if len(subset1.columns) != len(subset2.columns):
                    raise ValueError("Subsets should have the same amount of columns, unless something went very wrong")
                columns = list(subset1.sort_values("ID"))
                for column in range(5, len(columns)):
                    array_1 = np.array([[x, y] for x, y in zip(subset1.loc[:, 'ID'], subset1.loc[:, columns[column]])])
                    array_2 = np.array([[x, y] for x, y in zip(subset2.loc[:, 'ID'], subset2.loc[:, columns[column]])])
                    mae = similaritymeasures.mae(array_1, array_2)
                    measures = [columns[column], file1, file2, mae]
                    measures_names = ["Feature_name","Dataframe_1","Dataframe_2", "mae"]
                    mae_values.append(dict( (this,that) for this, that in zip(measures_names, measures)))
    df = pd.DataFrame(mae_values, columns=measures_names)
    df.to_csv(output_path)

calculate_mae(snakemake.input[0], snakemake.output[0])