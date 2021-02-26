from astropy.coordinates import SkyCoord
from astropy.table import Table, join, join_skycoord
import astropy.units as u
from .data_link import SQLReader, ASCIIReader
import pandas as pd

class Catalogue:
    READER_TYPE = {
        'sql': SQLReader,
        'ascii': ASCIIReader
    }

    def __init__(self, path, input_type, separator=' '):
        if not input_type in self.READER_TYPE:
            raise Exception('Unrecognized input type. use sql or ascii')
        self.reader = self.READER_TYPE[input_type](path, separator)

    @staticmethod
    def get_sky_coord(data):
        return SkyCoord(ra=data['ra'], dec=data['dec'], unit='deg', frame='icrs')

    def match_sources_within(self, catalogue, separation):
        df_1 = self.reader.data
        df_1 = Table.from_pandas(self.reader.data)

        df_2 = catalogue.reader.data
        df_2 = Table.from_pandas(catalogue.reader.data)
        
        df_1['sky_coord'] = self.get_sky_coord(df_1)
        df_2['sky_coord'] = self.get_sky_coord(df_2)

        df_1['idx_1'] = range(len(df_1))
        df_2['idx_2'] = range(len(df_2))
        
        return join(
            df_1[['sky_coord', 'idx_1']],
            df_2[['sky_coord', 'idx_2']],
            join_funcs={'sky_coord': join_skycoord(separation * u.arcsec)}
        )



#    c = Catalogue('~/Downloads/sky_dev_truthcat_v1.1.txt',  'ascii')
#    c1 = Catalogue('./output/sky_fits_sources_cat.sql', 'sql')                                                   
#    c1.match_sources_within(c, 1) #sources within 1 arcsec

