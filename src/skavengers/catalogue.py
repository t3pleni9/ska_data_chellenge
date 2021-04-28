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
        conflicting_sources = None

        for catalogue in catalogues:
            matched_sources = all_sources.match_sources_within(catalogue, separation)
            idx_1 = matched_sources['idx_1']
            idx_2 = matched_sources['idx_2']

            data = merged_catalogue.data

            conflicting_sources = matched_sources.to_pandas() \
                if conflicting_sources is None \
                else conflicting_sources.append(matched_sources.to_pandas())
            
            unmatched_idx_1 = [x for x in range(len(data)) if x not in idx_1]
            unmatched_idx_2 = [x for x in range(len(catalogue.data)) if x not in idx_2]

            merged_catalogue_data = data.iloc[unmatched_idx_1]
            unmatched_sources = catalogue.data.iloc[unmatched_idx_2]

            merged_catalogue_data = merged_catalogue_data.append(unmatched_sources)

            merged_catalogue = Catalogue(data=merged_catalogue_data)

        
        return merged_catalogue, conflicting_sources

    
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




def get_catalogues(sep):
   importlib.reload(catalogue)
   cat1 = catalogue.Catalogue('sofia_output_001_cat.sql', 'sql')
   cat2 = [catalogue.Catalogue(f"sofia_output_0{n:02}_cat.sql", 'sql') for n in range(2, 17)]
   return cat1.merge_catalogues(cat2, sep)    
