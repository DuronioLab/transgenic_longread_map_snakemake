import sys
import pandas as pd

def baseNameFormat(df, cols, delim="-"):
    # Input: dataframe, cols = colnames in order of basename
    # delim = string delimiter for basename column values

    # Checks:
    #
    # -- cols are members of df:
    nonExistColnames = [c for c in cols if c not in list(df)]
    if len(nonExistColnames) != 0:
        errorString = "Some columns in 'cols' don't exist:\n" + " ".join(nonExistColnames)
        raise ValueError(errorString)

    ###
    # make format string:
    baseNameFormat = delim.join(["{}" for c in cols])
    return (baseNameFormat)


def baseNameFromCols(df, cols, formatString, outputColName="baseName"):
    # Appends basename column to dataframe
    # Input:
    df[outputColName] = df[cols].apply(lambda x: formatString.format(*x), axis=1)
    return (df)


def addBaseName(df, cols, delim="-", outputColName="baseName"):
    # Takes sampleSheet & list of columns as input,
    # returns sampleSheet with 'basename' column appended
    formatString = baseNameFormat(df, cols, delim)
    df = baseNameFromCols(df, cols, formatString, outputColName=outputColName)
    return (df)


# def addExt(df, ext, baseNameCol = 'basename'):
#    # input sampleSheet, return extension of file

def makeSampleSheets(sampleInfoPath, baseNameCol, delim, fileDelimiter="\t"):
    df = pd.read_table(sampleInfoPath, delimiter=fileDelimiter)
    df = addBaseName(df, baseNameCol, delim)
    print(df.head())
    keep_cols = baseNameCol.copy()
    keep_cols.append(baseNameCol)
    return df


def main(sampleInfoPath, baseNameCol, delim='-'):
    # format & write sampleSheet
    sampleSheet = makeSampleSheets(sampleInfoPath, baseNameCol, delim)
    filename = "../sampleSheet.tsv"
    sampleSheet.to_csv(filename, sep="\t", index=False)


if __name__ == "__main__":
    # Drop call from argv
    sys.argv.pop(0)
    # Path to configFile is first argument
    path = sys.argv.pop(0)
    # Remaining function calls are the id variables for the basename
    idcols = sys.argv

    main(path, idcols, '-')