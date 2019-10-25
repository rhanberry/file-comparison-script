import pandas as pd
import os
import math
import sklearn
import argparse

#This script requires python 3.x to run properly, you will receive syntax errors in 2.x
#the first file in the arguments should be the 'test' file, and the second should be the 'reference' file

def print_statement(s):
    print("\n" + "*" * 100)
    print(s.upper())
    print("*" * 100 + "\n")


# read in data
def read_data(filepath, delimiter):
    print_statement("reading file %s" % (filepath))
    return pd.read_csv(filepath,
                       delimiter=delimiter,
                       na_values=['null'],
                       dtype=str)


# read counts
def get_data_counts(data):
    print_statement("getting record count")
    count_row = data.shape[0]  # gives number of row count
    count_col = data.shape[1]
    counts = {"row_count": count_row, "column_count": count_col}
    print("row and column counts:\n %s" % (counts))
    return counts


def compare_counts(counts1, counts2):
    print_statement("comparing record counts")
    return {"subtract_column_count": counts1.get("column_count") - counts2.get("column_count"),
            "subtract_row_count": counts1.get("row_count") - counts2.get("row_count")}


# read columns
def get_data_columns(data):
    print_statement("reading columns from dataset")
    columns = data.columns.tolist()
    print("the list of headers from data %s" % columns)
    cols = {"columns": ",".join(columns)}
    return cols


def compare_columns(cols1, cols2):
    print_statement("comparing columns from two datasets")
    return {"additional_columns": ",".join(
        list(set(cols1.get("columns").split(",")) - set(cols2.get("columns").split(",")))),
        "missing_columns": ",".join(
            list(set(cols2.get("columns").split(",")) - set(cols1.get("columns").split(","))))}


# calculate number nulls for columns
def get_data_number_of_nulls(data):
    print_statement("calculating the number of nulls for each column in the dataframe")
    column_number_nulls = {"number_nulls__" + key: val for key, val in dict(data.isnull().sum()).items()}
    print("the number of nulls for each column\n %s" % column_number_nulls)
    return column_number_nulls


def compare_null_counts(null_count1, null_count2):
    print_statement("comparing the number of nulls for each column in two dataframes")
    comparisons = {}
    for key, val in null_count1.items():
        for key2, val2 in null_count2.items():
            if (key == key2):
                comparisons[key] = val - val2
    comparisons = {**{"compare_" + key: math.nan for key in list(set(null_count1.keys()) - set(comparisons.keys()))},
                   **{"compare_" + key: math.nan for key in list(set(null_count2.keys()) - set(comparisons.keys()))},
                   **{"compare_" + key: val for key, val in comparisons.items()}}
    return comparisons

#uncomment below to compute average column lengths in your comparison files
# calculate average column length
# def get_data_average_column_lengths(data):
#     print_statement("calculating the average column length for each column in the dataframe")
#     column_number_not_nulls = dict(data.notna().sum())
#     print(column_number_not_nulls)
#     avg_column_lengths_prep = data.fillna('').astype(str).apply(lambda x: x.str.len()).sum(skipna=True)
#     avg_column_lengths = {col: val / col_sum_val
#                           for col, val in avg_column_lengths_prep.items()
#                           for col_sum, col_sum_val in column_number_not_nulls.items() if col == col_sum}
#     avg_column_lengths = {"average_column_length__" + key: val for key, val in avg_column_lengths.items()}
#     print("the average column length for each column\n %s" % avg_column_lengths)
#     return avg_column_lengths

#uncomment below to compare average column lengths in your comparison files
# def compare_average_column_lengths(avg_column_lengths1, avg_column_lengths2):
#     print_statement("comparing the average column lengths for each column in two dataframes")
#     comparisons = {}
#     for key, val in avg_column_lengths1.items():
#         for key2, val2 in avg_column_lengths2.items():
#             if (key == key2):
#                 if (val and val2):
#                     comparisons[key] = val - val2
#                 else:
#                     comparisons[key] = math.nan
#     comparisons = {
#         **{"compare_" + key: math.nan for key in list(set(avg_column_lengths1.keys()) - set(comparisons.keys()))},
#         **{"compare_" + key: math.nan for key in list(set(avg_column_lengths2.keys()) - set(comparisons.keys()))},
#         **{"compare_" + key: val for key, val in comparisons.items()}}
#     return comparisons


def compare_joined_dataframes(df1, df2, df_keys, returnSampledData=False):
    print_statement("performing an outer join on keys %s on the two dataframes" % (df_keys))
    df1.columns = [c.replace(' ', '_') for c in df1.columns]
    df2.columns = [c.replace(' ', '_') for c in df2.columns]
    joined_df = pd.merge(df1, df2, left_on=df_keys, right_on=df_keys, how='outer', indicator=True)
    additional_rows = joined_df.query('_merge == "left_only"').drop('_merge', 1)
    missing_rows = joined_df.query('_merge == "right_only"').drop('_merge', 1)
    matching_rows = joined_df.query('_merge == "both"').drop('_merge', 1)
    mismatching_columns = {}
    for column in list(set(data.columns.tolist() + data_compare.columns.tolist()) - set(df_keys)):
        if column + "_x" not in matching_rows.columns.tolist() or column + "_y" not in matching_rows.columns.tolist():
            mismatching_columns[column] = math.nan
        else:
            mismatch_rows_preped = matching_rows.dropna(subset=[column + "_x"]).dropna(subset=[column + "_y"])
            mismatching_columns[column] = mismatch_rows_preped.query("%s_x != %s_y" % (column, column))

    joined_counts = {
        "missing_rows": missing_rows.shape[0],
        "additional_rows": additional_rows.shape[0]}
    mismatched_counts = {"mismatch_count__" + col: math.nan if isinstance(val, float) else val.shape[0] for col, val
                         in mismatching_columns.items()}
    joined_counts_ret = {**joined_counts, **mismatched_counts}
    if (not returnSampledData):
        return joined_counts_ret
    else:
        sample_discrep = {
            "missing_rows": sklearn.utils.shuffle(missing_rows).head(100000),
            "additional_rows": sklearn.utils.shuffle(additional_rows).head(100000)}
        mismatched_counts_sample = {
            "mismatch_sample__" + col: math.nan if isinstance(val, float) else sklearn.utils.shuffle(val).head(100000) for
            col, val in mismatching_columns.items()}
        sample_discrep_ret = {**sample_discrep, **mismatched_counts_sample}
        return joined_counts_ret, sample_discrep_ret


# write out data
def write_data_summaries(output_filepath, all_summary_data):
    print_statement("writing out regression test data summary to file %s" % (output_filepath))
    all_summary_data = {key: [val] for key, val in all_summary_data.items()}
    print(all_summary_data)
    summary_df = pd.DataFrame.from_dict(all_summary_data)
    summary_df.index.name = "index"
    print(summary_df)
    if os.path.isfile(output_filepath):
        header = False
    else:
        header = True
    summary_df.to_csv(output_filepath, mode="a", header=header)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--delimiter", required=False,
                        help="this argument will determine the delimiter type in the input file/s, default will be a comma")
    parser.add_argument("-f", "--files", nargs='+', required=True,
                        help="either one or two files for which we are running a regression test on, if one is listed we will just collect metrics, if two are listed the two files will be compared")
    parser.add_argument("-pk", "--primary_keys", nargs='+', required=False,
                        help="the primary keys which will be used in the outer join to compare two files")
    parser.add_argument("-t", "--tag", required=True, help="tag used to logically separate the regression test output between logically different files")
    parser.add_argument("-s", "--printing_discrepancy_sample", required=False, action='store_true', help="boolean of whether we would like to collect a sample of the data of the mismatched data between two files")
    args = parser.parse_args()

    # set variables initially
    comparing = None
    delimiter = ","
    filepath = filepath_compare = None
    output_directory = "./output/{tag}"
    output_filepath = "./output/{tag}/data_summary.csv"
    output_filepath_compare = "./output/{tag}/data_summary_compare.csv"
    discrepancy_sample_directory = "./output/{tag}/sample"
    primary_keys = None
    tag = None
    printing_discrepancy_sample = False

    # using input variables now to set variables
    if args.delimiter:
        delimiter = args.delimiter
    filepath = args.files[0]
    if len(args.files) > 1:
        comparing = True
        filepath_compare = args.files[1]
    if args.primary_keys:
        primary_keys = args.primary_keys
    tag = args.tag
    output_directory = output_directory.replace("{tag}", tag)
    output_filepath = output_filepath.replace("{tag}", tag)
    output_filepath_compare = output_filepath_compare.replace("{tag}", tag)
    discrepancy_sample_directory = discrepancy_sample_directory.replace("{tag}", tag)
    if args.printing_discrepancy_sample:
        printing_discrepancy_sample = True

    # printing the regression test plan
    if not comparing:
        print(
                    "We are going to perform regression testing by collecting metrics on file %s and we will write our output to %s" % (
            filepath, output_filepath))
    else:
        print(
                    "We are going to perform regression testing by collecting metrics and comparing files %s and %s and we will write our output to %s" % (
            filepath, filepath_compare, output_filepath_compare))
    print("will we be using delimiter '%s'" % (delimiter))
    if primary_keys:
        print("we will be using the following primary keys: %s" % (primary_keys))
    else:
        print("we will not be using any primary keys")

    # performing tests:

    # reading the actual data
    data = read_data(filepath, delimiter)
    if comparing:
        data_compare = read_data(filepath_compare, delimiter)
    else:
        data_compare = None

    # collecting the file counts
    counts = get_data_counts(data)
    if comparing:
        counts_compare = get_data_counts(data_compare)
        comparison_counts = compare_counts(counts, counts_compare)
    else:
        counts_compare = comparison_counts = None

    # collecting column names
    columns = get_data_columns(data)
    if comparing:
        columns_compare = get_data_columns(data_compare)
        comparison_columns = compare_columns(columns, columns_compare)
    else:
        columns_compare = comparison_columns = None

    # collecting the number of nulls in each column
    number_of_nulls = get_data_number_of_nulls(data)
    if comparing:
        number_of_nulls_compare = get_data_number_of_nulls(data_compare)
        comparisons_null_counts = compare_null_counts(number_of_nulls, number_of_nulls_compare)
    else:
        number_of_nulls_compare = comparisons_null_counts = None

    #uncomment below to pring column length comparisons in output file
    #  # collecting the average column length of each column
    # average_column_lengths = get_data_average_column_lengths(data)
    # if comparing:
    #     average_column_lengths_compare = get_data_average_column_lengths(data_compare)
    #     comparisons_average_column_lengths = compare_average_column_lengths(average_column_lengths,
    #                                                                         average_column_lengths_compare)
    # else:
    #     average_column_lengths_compare = comparisons_average_column_lengths = None

    # comparing the joined datasets from two files (if comparing)
    if comparing:
        if printing_discrepancy_sample:
            joined_counts, joined_sample = compare_joined_dataframes(data, data_compare, primary_keys, True)
        else:
            joined_counts = compare_joined_dataframes(data, data_compare, primary_keys)
    else:
        joined_counts = joined_sample = None

    # writing the regression test output to a file
    # making output directory if not exists
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)
    # printing for a single file
    if not comparing:
        all_summary_data = {**counts, **columns, **number_of_nulls} # add **average_column_lengths if you wish to compare those values in output file
        write_data_summaries(output_filepath, all_summary_data)
    # printing for two files (ie comparison)
    else:
        all_summary_data_compare = {**comparison_columns, **comparison_columns, **comparisons_null_counts,
                                    **joined_counts} # add **comparisons_average_column_lengths if you wish to compare those values in output file
        write_data_summaries(output_filepath_compare, all_summary_data_compare)
        # if we are comparing and we want to print the discrepancy we will do so here
        if printing_discrepancy_sample:
            print("printing out additional, missing, and mismatching data to sample directory")
            if not os.path.isdir(discrepancy_sample_directory):
                os.makedirs(discrepancy_sample_directory)
            for file, df in joined_sample.items():
                if type(df) != pd.core.frame.DataFrame or df.empty:
                    continue
                else:
                    df.index.name = "index"
                    filepath = os.path.join(discrepancy_sample_directory, file)
                    if os.path.isfile(filepath):
                        header = False
                    else:
                        header = True
                    df.to_csv(filepath, mode="a", header=header, sep=delimiter)
