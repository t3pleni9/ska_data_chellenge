import pandas as pd
import sqlite3

class SQLReader:
    def __init__(self, path, separator=None):
        self.connection = sqlite3.connect(":memory:")
        self.cursor = self.connection.cursor()

        with open(path) as sql_file:
            sql_as_string = sql_file.read()

        sql_as_string = sql_as_string.split(
            'SET SQL_MODE = "NO_AUTO_VALUE_ON_ZERO";'
        )[1].replace(
            'KEY (`id`)', ''
        ).replace(
            'NOT NULL', ''
        ).replace(
            'nan', 'NULL'
        ).replace(
            'PRIMARY ,', 'PRIMARY KEY (`id`)'
        ).replace("COMMENT='SoFiA source catalogue; created with SoFiA version 2.3.0';", '')

        create_table = sql_as_string.split('INSERT')[0]
        insert_values = 'INSERT' + sql_as_string.split('INSERT')[1]

        self.table_name = 'SoFiA-Catalogue'
        self.cursor.executescript(create_table)
        self.cursor.executescript(insert_values)

    @property
    def data(self):
        return pd.read_sql(f'SELECT * from `{self.table_name}`', self.connection)

