import pickle

__all__ = ['dump_stuff', 'read_pickle']
# ------------------------------------------------------------------------------

def dump_stuff(data_to_save, filename, outdir):
    """

    This function dumps the data_to_save as a pickle in outdir.

    Required Inputs
    ---------------
    * data_to_save: data to save.
    * filename: str: name of the file
                    (if .pickle extension is missing, it will be added).
    * outdir: str: path to the folder where the data is to be saved.

    """
    if not filename.__contains__('.pickle'):
        filename += '.pickle'
    with open('%s/%s'%(outdir, filename), 'wb') as handle:
        pickle.dump(obj=data_to_save, file=handle,
                    protocol=pickle.HIGHEST_PROTOCOL)

# ------------------------------------------------------------------------------
def read_pickle(filename):
    """

    This function reads the pickle: filename.

    Required Inputs
    ---------------
    * filename: str: name of the pickle file (with extension).

    Returns
    --------
    * Pickled data

    """
    with open(filename, 'rb') as handle:
        return pickle.load(file=handle)
