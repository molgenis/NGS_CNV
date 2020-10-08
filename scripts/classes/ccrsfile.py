class CcrsFile:
    def __init__(self, fileloc):
        self.ccrs_los = fileloc
        self.header_data = {}
        self.ccrs_data = {}

    def read_ccrs_file(self):
        header_tags = ("@", "CONTIG")
        try:
            with open(self.fileloc, 'r') as infile:
                for fileline in infile:
                    if fileline.startswith(header_tags):
                        
                    else:
                        
        except IOError:
            

    def write_header(self, outfileloc, headerdata):
        self.__write_ccrs__(outfileloc, True, headerdata, False, None)

    def write_entries(self, outfileloc, ccrs_entries):
        self.__write_ccrs__(outfileloc, False, None, True, ccrs_entries)

    def write_ccrs_file(self, outfileloc, headerdata=self.header_data, ccrsdata=self.ccrs_data):
        self.__write_ccrs__(outfileloc, True, headerdata, True, ccrsdata)

    def __write_ccrs__(outfileloc, writeheader=True, ccrsheader=self.header_data, writebody=True, ccrsbody=self.ccrs_data):
        file_written = False
        if writeheader or writebody:
            try:
                with open(outfileloc, 'w') as outfile:
                    if write_header:
                        if ccrsheader is not None:
                            for headerfield in ccrsheader:
                                
                    if write_body:
                        if ccrsbody is not None:
                            
                file_written = True
            except IOError:
                print("")
            finally:
                return file_written
        return file_written
