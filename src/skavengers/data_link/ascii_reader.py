import pandas as pd

class ASCIIReader:
    def __init__(self, path, separator=' '):
        self.path = path
        self.separator = separator

    @property
    def data(self):
        return pd.read_csv(self.path, sep=self.separator)
