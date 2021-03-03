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

    def __init__(self, path=None, input_type=None, separator=' ', data=None):
        if data is None:
            if not input_type in self.READER_TYPE:
                raise Exception('Unrecognized input type. use sql or ascii')
            reader = self.READER_TYPE[input_type](path, separator)
            self.data = reader.data
        else:
            self.data = data
        
        

    @staticmethod
    def get_sky_coord(data):
        return SkyCoord(ra=data['ra'], dec=data['dec'], unit='deg', frame='icrs')

    def merge_catalogues(self, catalogues, separation):
        merged_catalogue = self

        for catalogue in catalogues:
            matched_sources = merged_catalogue.match_sources_within(catalogue, separation)
            data = merged_catalogue.data
            
            idx_2 = matched_sources['idx_2']

            unmatched_idx_2 = [x for x in range(len(catalogue.data)) if x not in idx_2]
            unmatched_sources = catalogue.data.iloc[unmatched_idx_2]

            merged_catalogue_data = data.append(unmatched_sources)

            merged_catalogue = Catalogue(data=merged_catalogue_data)

        
        return merged_catalogue

    
    def match_sources_within(self, catalogue, separation):
        df_1 = Table.from_pandas(self.data)
        df_2 = Table.from_pandas(catalogue.data)
        
        df_1['sky_coord'] = self.get_sky_coord(df_1)
        df_2['sky_coord'] = self.get_sky_coord(df_2)

        df_1['idx_1'] = range(len(df_1))
        df_2['idx_2'] = range(len(df_2))
        
        
        return join(
            df_1[['sky_coord', 'idx_1']],
            df_2[['sky_coord', 'idx_2']],
            join_funcs={'sky_coord': join_skycoord(separation * u.arcsec)}
        )
