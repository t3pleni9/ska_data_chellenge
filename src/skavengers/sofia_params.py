class SofiaParams:
    def __init__(self, par_filename, par_file_directory, output_filename, output_directory):
        self.par_filename = par_filename
        self.par_file_directory = par_file_directory
        self.output_filename = output_filename
        self.output_directory = output_directory
        self.output_type = 'sql'

    @property
    def output_extension(self):
        extensions = {
            'sql': 'sql',
            'ascii': 'txt'
        }

        return extensions[self.output_type]
