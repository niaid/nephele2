from enum import Enum
import os
import re
import traceback
import json
import pandas
import numpy
import sh
from nephele2.infra.utils.nephele_logger import create_logger

logger = create_logger()


class Headers(Enum):
    FIRST_COL = '#SampleID'
    LAST_COL = 'Description'
    TREATMENT = 'TreatmentGroup'
    FORWARD_READ = 'ForwardFastqFile'
    REVERSE_READ = 'ReverseFastqFile'

    BASIC_HEADERS = [FIRST_COL, TREATMENT, LAST_COL]
    PE_HEADERS = [FIRST_COL, FORWARD_READ, REVERSE_READ, TREATMENT, LAST_COL]
    SE_DEMULTIPLEX_HEADERS = [FIRST_COL, FORWARD_READ, TREATMENT, LAST_COL]
    SE_WGS_HEADERS = [FIRST_COL, FORWARD_READ]
    PE_WGS_HEADERS = [FIRST_COL, FORWARD_READ, REVERSE_READ]
    SE_QC_HEADERS = [FIRST_COL, FORWARD_READ]
    PE_QC_HEADERS = [FIRST_COL, FORWARD_READ, REVERSE_READ]
    DS_ANALYSIS_HEADERS = [FIRST_COL, FORWARD_READ, REVERSE_READ, TREATMENT,
                           LAST_COL]


class MapType(Enum):
    STANDARD = (1, Headers.BASIC_HEADERS.value)
    PAIRED_END = (2, Headers.PE_HEADERS.value)
    SE_DEMULTIPLEX = (4, Headers.SE_DEMULTIPLEX_HEADERS.value)
    PE_WGS = (5, Headers.PE_WGS_HEADERS.value)
    SE_WGS = (6, Headers.SE_WGS_HEADERS.value)
    PE_QC = (7, Headers.PE_QC_HEADERS.value)
    SE_QC = (8, Headers.SE_QC_HEADERS.value)
    DS_ANALYSIS = (9, Headers.DS_ANALYSIS_HEADERS.value)

    def __init__(self, index, required_headers):
        self.index = index
        self.required_headers = required_headers


class map_validator(object):

    def __init__(self, filepath, map_type, allowed_exts, column_map=None, headers_list=None, listing_file=None):
        """Instantiates the map validator by reading in the map file using
        pandas and setting class variables related to the file and error lists.

        The filepath and map_type arguments are required.  If the column_map
        argument is included, filepath is assumed to refer to the file the data
        from the map will be written to.  Otherwise, filepath is assumed to be
        the file containing the data to be read. The headers_list argument is
        only valid when used with the column_map argument and dictates the
        order of the columns in the resulting dataframe.  See
        https://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.html
        allowed_exts is required for validating the read file names. It is
        required, and should be a list.  listing_file is required for
        validating file existence if the files are coming from a remote host.
        Since we'll download these files at runtime, this validator expects a
        file containing the names of the files found in the remote location.
        """
        self.errs = {}
        self.warnings = {}
        try:
            # Read in the map file, either from a file or from a dict, set the
            # self path variables
            if column_map:
                self.df = pandas.DataFrame(column_map, columns=headers_list)
                self.df.replace('', numpy.nan, inplace=True)
                self._create_corrected_filename(filepath)
            else:
                self._read_map_file(filepath)
            # set the type of analysis
            if not isinstance(map_type, MapType):
                raise TypeError(
                    "The map_type must be an instance of MapType Enum.")
            else:
                self.map_type = map_type
            if self.map_type == MapType.DS_ANALYSIS:
                self.ignore_fq_files = True
            else:
                self.ignore_fq_files = False
            # set the allowed file extensions
            self.allowed_file_extensions = allowed_exts
            # set the location of the .listing file if uploading from remote host later
            self.remote_file_list = listing_file
        except Exception:
            logger.error("Failed to init map validation " +
                         str(traceback.format_exc()))
            raise

    @staticmethod
    def is_ascii(test_string):
        """Tests if a string is ASCII or not.
        
        Args:
            test_string (str): the string to check
        
        Returns:
            boolean: True if ASCII else False
        """
        try:
            test_string.encode('ascii')
        except UnicodeEncodeError:
            return False
        return True

    def _create_corrected_filename(self, filepath):
        """Takes a file name and returns it with _corrected added before the extension.
        
        Args:
            filepath (str): the filepath to alter
        
        Returns:
            str: the modified filepath or the filepath if it already contained _corrected.
        
        Raises:
            Exception: if there was an error parsing the original filename
        """
        try:
            self.dir = os.path.dirname(filepath)
            filebasename, ext = os.path.splitext(os.path.basename(filepath))
            if re.match('.*_corrected$', filebasename):
                self.fixed_filename = filepath
            else:
                self.fixed_filename = self.dir+"/"+filebasename+"_corrected.txt"
        except Exception:
            raise

    def _read_map_file(self, filepath):
        """Reads in a text or Excel file using pandas.  For non-Excel files, pandas will
        read the file in universal newline mode and attempt to auto-determine the separator.
        """
        try:
            filebasename, ext = os.path.splitext(os.path.basename(filepath))
            self._create_corrected_filename(filepath)
            # read the file (pandas may handle excel for us)
            if 'xl' in ext.lower():
                # TODO: this is working, but I don't know how to unit test it
                self.df = pandas.read_excel(filepath)
            else:
                # read from csv
                # TODO: this seems to work with an actual file, but I don't know how to mock this
                self.df = pandas.read_csv(
                    open(filepath, newline=None), sep=None, engine='python', index_col=False)
        except Exception:
            logger.error("Failed to read map file " +
                         str(traceback.format_exc()))
            raise

    # FIXME: needs some exception handling in case write fails, we'd need to add an error to the validation result object
    def _write_map_file(self):
        """Writes out the contents of the dataframe to a file named
        <original_filename>_corrected.txt.  It will not print the index
        and will use tabs as the separator.
        """
        self.df.to_csv(self.fixed_filename, sep='\t', index=False)

    def _add_err(self, err_type, msg, col_name, row_list):
        """Adds a new error message to the dict of errors.  Each error added
        must have it's own, unique key.
        """
        self.errs.update(
            {err_type: {"msg": msg, "col": col_name, "rows": row_list}})

    def _add_warning(self, warn_type, msg, col_name, row_list):
        """Adds a new error message to the dict of warnings.  Each warning added
        must have it's own, unique key.
        """
        self.warnings.update(
            {warn_type: {"msg": msg, "col": col_name, "rows": row_list}})

    def _get_validation_result(self):
        """Triggers the writing of the corrected file and creates a JSON
        containing information about the state of the mapping file, including whether
        it is valid and any errors and warnings generated duringin validation. It also
        includes a list of metadata columns.

        returns:
        {is_valid: True/False,
         corrected_file: filename,
         Errors: {err_name: { msg: "",
                              col: col_name,
                              rows: []
                            }
                 },
         Warnings: {},
         metadata: []
        }
        """
        # NOTE: adding the list of metadata columns to this object may cause an exception
        # to be logged if the file is junk, but everything will still work
        self._write_map_file()
        validation_result = {}
        if self.errs:
            validation_result['is_valid'] = False
        else:
            validation_result['is_valid'] = True
        validation_result['corrected_file'] = self.fixed_filename
        validation_result['Errors'] = self.errs
        validation_result['Warnings'] = self.warnings
        # TODO: do we want to fetch this here, or should it be set somewhere?
        validation_result['metadata'] = self._get_metadata_columns()
        return json.dumps(validation_result)

    def render_table(self):
        df_styled = self.df.copy()
        df_styled = df_styled.fillna('')
        styler = df_styled.style
        hlist = self.df.columns.tolist()
        # NOTE: we can't style the header here, so we'll do it in JQuery
        for err in self.errs:
            col = self.errs[err]['col']
            if 'header' in self.errs[err]['rows']:
                rows = [x for x in self.errs[err]['rows'] if x != 'header']
            else:
                rows = self.errs[err]['rows']
            if rows:
                styler = styler.set_properties(
                    subset=pandas.IndexSlice[rows, hlist[col]], **{'background-color': 'red'})

        for warning in self.warnings:
            col = self.warnings[warning]['col']
            if 'header' in self.warnings[warning]['rows']:
                # FIXME: style the header
                rows = [x for x in self.warnings[warning]
                        ['rows'] if x != 'header']
            else:
                rows = self.warnings[warning]['rows']
            if rows:
                styler = styler.set_properties(
                    subset=pandas.IndexSlice[rows, hlist[col]], **{'background-color': 'yellow'})
        return styler.render()

    def _validate_headers(self):
        """Checks the header row to ensure that it meets the requirements
        listed in the git wiki.  Automatically removes any spaces, and attempts
        to rescue the #SampleID header from common errors and rescue the other
        headers from capitalization errors.  Adds any errors to the error dict.
        """
        # fix any headers with spaces, just remove them
        self.df.rename(columns=lambda x: x.replace(" ", ""), inplace=True)

        # absolutely horrible hack that allows DS analysis to have SE headers
        # or PE headers
        hlist = self.df.columns.tolist()
        if self.map_type == MapType.DS_ANALYSIS:
            if hlist[2].lower() == 'reversefastqfile':
                self.map_type = MapType.PAIRED_END
            else:
                self.map_type = MapType.SE_DEMULTIPLEX

        requiredh = self.map_type.required_headers

        hlist = self.rename_duplicate_headers()

        # check for empty headers
        pattern = re.compile(r'^Unnamed:\s?\d+$')
        for label in hlist:
            if not label or label.isspace() or pattern.match(label):
                self._add_err("empty_column_label", "Column labels cannot be blank.", hlist.index(
                    label), ["header"])
        # verify the first column header
        if hlist[0] != requiredh[0]:
            if hlist[0].lower() == 'sampleid' or hlist[0].lower() == '#sampleid':
                self.df.rename(
                    columns={hlist[0]: Headers.FIRST_COL.value}, inplace=True)
                self._add_warning("bad_sample_header",
                                  "Corrected bad column name '" +
                                  hlist[0]+"' to "+Headers.FIRST_COL.value+".",
                                  0, ["header"])
            else:
                self._add_err("bad_sample_header",
                              "Bad column name '"+hlist[0]+"'. The first column should be '"+Headers.FIRST_COL.value+"'. \
                              Correct headers for this analysis are <strong>"+str(requiredh)+"</strong>. \
                              All metadata columns should be listed between TreatmentGroup and Description. \
                              See the notes at the bottom of this page for a link to a correct template.",
                              0, ["header"])
        # make sure all other required column headers are present
        # (for non-WGS, don't validate the last one, since there can be many columns
        # between the last two required headers.
        if self.map_type in [MapType.SE_WGS, MapType.PE_WGS, MapType.SE_QC,
                             MapType.PE_QC]:
            last = len(requiredh)
        else:
            last = len(requiredh)-1
        if last > len(hlist):
            self._add_err("not_enough_columns", "Missing required columns. Correct \
                          headers for this analysis are <strong>"+str(requiredh)+"</strong>. \
                          All metadata columns should be listed between TreatmentGroup and Description \
                          for non-WGS analyses, or after the required columns for WGS analyses. \
                          See the notes at the bottom of this page for a link to a correct template.",
                          len(hlist)-1, ["header"])
        else:
            for i in range(1, last):
                if not re.fullmatch(requiredh[i], hlist[i]):
                    if re.fullmatch(requiredh[i].lower(), hlist[i].lower()):
                        self.df.rename(
                            columns={hlist[i]: requiredh[i]}, inplace=True)
                        self._add_warning("bad_column_header_"+hlist[i],
                                          "Corrected bad column name '" +
                                          hlist[i]+"' to "+requiredh[i]+".",
                                          i, ["header"])
                    else:
                        self._add_err("bad_column_header_"+hlist[i],
                                      "Bad column name '"+hlist[i]+"'. Column "+str(i+1)+" should be '"+requiredh[i]+"'. \
                                      Correct headers for this analysis are <strong>"+str(requiredh)+"</strong>. \
                                      All metadata columns should be listed between TreatmentGroup and Description \
                                      for non-WGS analyses, or after the required columns for WGS analyses. \
                                      See the notes at the bottom of this page for a link to a correct template.",
                                      i, ["header"])

        self.validate_metadata_column_headers(hlist, len(requiredh)-1, len(hlist)-1)

        # For non-WGS only, make sure the last column is Description
        if self.map_type not in [MapType.PE_WGS, MapType.SE_WGS, MapType.PE_QC, MapType.SE_QC] and hlist[-1] != Headers.LAST_COL.value:
            if re.fullmatch(hlist[-1].lower(), Headers.LAST_COL.value.lower()):
                self.df.rename(
                    columns={hlist[-1]: Headers.LAST_COL.value}, inplace=True)
                self._add_warning("bad_final_header",
                                  "Corrected bad column name '" +
                                  hlist[-1]+"' to "+Headers.LAST_COL.value+".",
                                  len(hlist)-1, ["header"])
            else:
                self._add_err("bad_final_header",
                              "Bad column name '" +
                              hlist[-1]+"'. The last column should be '" +
                              Headers.LAST_COL.value+"'.",
                              len(hlist)-1, ["header"])

        # ReversePrimer is only used for multiplex data, which we no longer use
        if "ReversePrimer" in hlist:
            self._add_err("invalid_col_header",
                          "ReversePrimer should only be included for Paired End Multiplex analysis. \
                          Nephele no longer supports multiplex data sets.",
                          hlist.index("ReversePrimer"), ["header"])

        if self.map_type in [MapType.SE_WGS, MapType.SE_QC] and "ReverseFastqFile" in hlist:
            self._add_err("invalid_pe_col_header", "ReverseFastqFile is only used for paired end data. \
            Single end data does not require this column, please remove it and try again.",
                          hlist.index("ReverseFastqFile"), ["header"])

    def rename_duplicate_headers(self):
        """Checks for duplicate headers (they should be unique)
        This is only a problem when loading df not from a file. For consistency,
        we'll just rename them, which is what pandas does when it reads in from a file anyway.
        """
        try:
            cols = pandas.Series(self.df.columns)
            for dup in self.df.columns[self.df.columns.duplicated()].unique():
                cols[self.df.columns.get_loc(dup)] = [
                    dup+'.'+str(d_idx) if d_idx != 0 else dup for d_idx in range(self.df.columns.get_loc(dup).sum())]
            self.df.columns = cols
            return self.df.columns.tolist()
        except Exception:
            self._add_err("failed_to_rename_columns", "Error reading in data, failed to make column headers unique", None, None)

    def validate_metadata_column_headers(self, hlist, start, end):
        """Checks all metadata column headers to ensure no bad chars (see 1-d)"""
        for i in range(start, end):
            pattern = r'^[\w\s\.\-\+%;:,\/]+$'
            if not re.match(pattern, hlist[i]):
                self._add_err("bad_column_header_"+hlist[i],
                              "Bad column name '" +
                              hlist[i] +
                              "'. Column names must be composed of only alphanumeric, underscore (\"_\"), period (\".\"), minus sign (\"-\"), plus sign (\"+\"), percentage (\"%\"), space (\" \"), semicolon (\";\"), colon (\":\"), comma (\",\"), and/or forward slash (\"/\") characters.",
                              i, ["header"])

    def _get_metadata_columns(self):
        """Returns the list of metadata columns, which is defined as everything
        beginning at, and including, the TreatmentGroup column through, but not
        including, the Description column (for non-WGS pipes), or from the last
        required column header to the last column (for WGS pipes)."""
        # is this a reasonable way to handle this?
        try:
            hlist = self.df.columns.tolist()
            if self.map_type in [MapType.SE_WGS, MapType.PE_WGS, MapType.SE_QC, MapType.PE_QC]:
                start = len(self.map_type.required_headers)
                # from the start index to the end of the list
                metadata_cols = hlist[start:]
            else:
                start = hlist.index(Headers.TREATMENT.value)
                # from the start index to, but not including, the last index of the list
                metadata_cols = hlist[start:-1]
            return metadata_cols
        except ValueError:
            self._add_err("missing_metadata_columns",
                          "Couldn't find the required column header TreatmentGroup, which must come before Description.", len(hlist)-1, ["header"])
            return []

    def _strip_quotes(self):
        """
        Removes all double-quotes from the dataframe.
        """
        self.df.columns = self.df.columns.str.replace('\"', '')
        self.df.replace({'\"': ''}, regex=True, inplace=True)

    def _remove_empty_cols(self):
        """
        Gets rid of any columns that are completely empty and have no headers.
        """
        hlist = self.df.columns.tolist()
        pattern = re.compile(r'^Unnamed:\s?\d+$')
        # check if any headers are blank
        for label in hlist:
            if not label or label.isspace() or pattern.match(label):
                if not (len(self.df[label].value_counts() > 0)):
                    # remove any completely empty columns
                    self.df.drop(label, axis=1, inplace=True)

    def _has_empty_values(self, col_name):
        """Returns a list of rows with empty values in the given column."""
        null_indexes = numpy.where(pandas.isnull(self.df[col_name]) == True)
        return self.df.iloc[null_indexes].index.values.tolist()

    def _validate_sample_ids(self):
        """Checks values in the #SampleID column to ensure they meet the requirements
        in the git wiki.  If any values are empty, checking is halted and the errors
        are returned.  All leading and trailing spaces are automatically removed.
        """
        # validate SampleID column
        # check for NaN values
        null_ids = self._has_empty_values(Headers.FIRST_COL.value)
        if len(null_ids) > 0:
            self._add_err("null_values_"+Headers.FIRST_COL.value,
                          "Entries in the "+Headers.FIRST_COL.value+" column cannot be blank",
                          0, null_ids)
        else:
            # make sure this column is interpreted as a string or this blows up
            self.df[Headers.FIRST_COL.value] = self.df[Headers.FIRST_COL.value].apply(
                str)
            # add warning if leading or trailing space (or we could remove it and warn)
            # this removes the spaces, but I don't know how to get these
            # indexes, so this may be a silent fix
            self.df[Headers.FIRST_COL.value] = self.df[Headers.FIRST_COL.value].str.strip()
            # validate that sample IDs are unique
            indexes = numpy.where(self.df.duplicated(
                Headers.FIRST_COL.value, keep=False) == True)
            # these are the row indexes
            duplicated_ids = self.df.iloc[indexes].index.values.tolist()
            if len(duplicated_ids) > 0:
                self._add_err("unique_ids",
                              "Sample IDs must be unique",
                              0, duplicated_ids)
            if self.map_type in [MapType.SE_WGS, MapType.PE_WGS]:
                # validate only alphanumeric, dash '-', and underscore '_', not '.' so we can't use \w
                alphanum_pattern = r'^[a-zA-Z0-9\_\-]+$'
            elif self.map_type in [MapType.SE_QC, MapType.PE_QC]:
                alphanum_pattern = r'^[a-zA-Z0-9\_\-\.]+$'
            else:
                # validate only alphanumeric and '.', not underscore '_', so we can't use \w
                alphanum_pattern = r'^[a-zA-Z0-9\.]+$'
            # keep the indexes that don't match, they're bad
            alphanum_indexes = numpy.where(
                self.df[Headers.FIRST_COL.value].str.match(alphanum_pattern) == False)
            # be sure to convert ALL of the arrays to Python native type lists (using tolist())
            bad_ids = self.df.iloc[alphanum_indexes].index.values.tolist()
            if len(bad_ids) > 0:
                msg = "Sample IDs can only contain alphanumeric and dot '.' characters."
                if self.map_type in [MapType.SE_WGS, MapType.PE_WGS]:
                    msg = "Sample IDs can only contain alphanumeric, dash '-', and underscore '_' characters."
                if self.map_type in [MapType.SE_QC, MapType.PE_QC]:
                    msg = "Sample IDs can only contain alphanumeric, periods '.', dash '-', and underscore '_' characters."
                self._add_err("bad_ids",
                              msg,
                              0, bad_ids)
            # can't use str.isdigit() to check if it's numeric because it doesn't recognize decimals

    def _is_valid_categorical_data(self, col_name):
        """Checks that the column has at least two unique values, and that all values
        conform to the character set allowed in the metadata columns.  Converts blanks to NA.
        """
        # blanks are not allowed, they should be entered as "NA"
        # replace all empty values with NaN
        self.df[col_name].replace(
            to_replace=r'^\s+$', value=numpy.NaN, regex=True, inplace=True)
        # add indexes to warnings to tell users their blanks have been altered
        blank_entries = self._has_empty_values(col_name)
        self.df[col_name] = self.df[col_name].replace(numpy.NaN, "NA")
        if len(blank_entries) > 0:
            self._add_warning("added_NA_"+col_name,
                              "Converted empty entries to NA.",
                              self.df.columns.get_loc(col_name), blank_entries)
        # Don't apply str until after converting NaNs
        # Treat the column as a string now or the match fails on numeric data
        self.df[col_name] = self.df[col_name].apply(str)
        # Remove leading and trailing whitespace
        self.df[col_name] = self.df[col_name].str.strip()
        # must have at least 2 unique values
        if len(self.df[col_name].value_counts()) < 2:
            self._add_err("too_few_categories_"+col_name,
                          "Metadata columns must have at least two unique values per column. " +
                          col_name+" does not meet this requirement.",
                          self.df.columns.get_loc(col_name), self.df.index.values.tolist())
        else:
            #check for non-ASCII data
            self.check_for_non_ascii(col_name)
            # Only alphanumerics and a few special chars are allowed in these columns
            pattern = r'^[\w\s\.\-\+%;:,\/]+$'
            indexes = numpy.where(
                self.df[col_name].str.match(pattern) == False)
            bad_seqs = self.df.iloc[indexes].index.values.tolist()
            if len(bad_seqs) > 0:
                self._add_err("bad_entry_"+col_name,
                              col_name +
                              " entries can only contain alphanumeric, underscore (“_”), period (”.”), minus sign (“-”), plus sign (“+”), percentage (“%”), space (” ”), semicolon (”;”), colon (”:”), comma (”,”), and/or forward slash (“/”) characters. Entries must not begin or end with a space.",
                              self.df.columns.get_loc(col_name), bad_seqs)

    def check_for_non_ascii(self, col_name):
        """Checks for non-ASCII data.
        
        Args:
            col_name (str): name of the column to check
        """
        try:
            non_ascii_indexes = numpy.where(
                [map_validator.is_ascii(x) == False for x in self.df[col_name]])
            non_ascii_strs = self.df.iloc[non_ascii_indexes].index.values.tolist()
            if len(non_ascii_strs) > 0:
                self._add_err("non_ascii_str_"+col_name,
                              col_name+" entries cannot contain non-ASCII characters or analysis steps will fail.",
                              self.df.columns.get_loc(col_name), non_ascii_strs)
        except Exception:
            self._add_err("validation_failure", "Unexpected error encountered while checking "+col_name+" for non-ASCII characters. Could not validate column data.", self.df.columns.get_loc(col_name), ["header"])

    def _validate_read_names(self, col_name):
        """Checks that file names are all unique or all identical, and contain
        only allowed characters. Any non- .fq or .fastq extension will throw an error here.
        """
        # make sure there are no empty values
        null_ids = self._has_empty_values(col_name)
        if len(null_ids) > 0:
            self._add_err("null_values_"+col_name,
                          "Entries in the "+col_name+" column cannot be blank",
                          0, null_ids)
        else:
            self.df[col_name] = self.df[col_name].str.strip()
            # Check for duplicates
            indexes = numpy.where(self.df.duplicated(
                col_name, keep=False) == True)
            bad_names = self.df.iloc[indexes].index.values.tolist()
            if len(self.df[col_name].value_counts()) != 1 and len(bad_names) != 0:
                self._add_err("bad_file_names_"+col_name,
                              "Identifiers in the file name columns must be either unique, or identical in the entire column",
                              self.df.columns.get_loc(col_name), bad_names)
            # only accept allowed file extensions, and names can't have spaces
            # then check that the names have valid extensions
            # this builds the pattern "\.ext1|\.ext2" etc.
            ext_set = "|".join([r"\."+ext for ext in [s.replace('.', r'\.')
                                                      for s in self.allowed_file_extensions]])
            pattern = r'\S*(?:'+ext_set+')$'
            indexes = numpy.where(
                self.df[col_name].str.match(pattern) == False)
            bad_names = self.df.iloc[indexes].index.values.tolist()
            if len(bad_names) > 0:
                self._add_err("bad_filename_characters_"+col_name,
                              "The file names should not contain spaces and should end in " +
                              ", ".join(self.allowed_file_extensions),
                              self.df.columns.get_loc(col_name), bad_names)

    def check_files_exist(self, col_name):
        """Checks to make sure that every file listed in the column exists in the directory.
        Adds an error to the error dict if files aren't found.
        """
        # Flag anything that is a null value and don't include those in the existence check
        null_indexes = numpy.where(pandas.isnull(self.df[col_name]))
        null_vals = self.df.iloc[null_indexes].index.values.tolist()
        if len(null_vals) > 0:
            self._add_err("null_values_"+col_name,
                          "There are null values present in your "+col_name+" column, which is not allowed. \
                          File names must match exactly the names of the files you uploaded for analysis.",
                          self.df.columns.get_loc(col_name), null_vals)
        # returns array of True/False: [os.path.isfile(self.dir+"/"+x) for x in df[col_name]]
        # pandas Series returns a vector(?)/dataframe(?),
        # since we didn't pass the index option, a numeric, 0 base index was created
        # numpy.where gets the indexes of the values that match the conditional
        filenames = self.df[col_name][~pandas.isnull(
            self.df[col_name])]  # don't include nulls
        if self.remote_file_list:
            indexes = numpy.where(pandas.Series(
                [self._check_remote_listing(x) for x in filenames]) == False)
        else:
            indexes = numpy.where(pandas.Series(
                [os.path.isfile(self.dir+"/"+x) for x in filenames]) == False)
        missing = self.df.iloc[indexes].index.values.tolist()
        if len(missing) > 0:
            self._add_err("missing_files_"+col_name,
                          "Some files listed in the mapping file were not found in your data set.",
                          self.df.columns.get_loc(col_name), missing)

    def _check_remote_listing(self, filename):
        """
        Called only from check_files_exist when validating remote files from a .listing file.

        params:    filename    name of the file we're looking for
        returns:   True if the file is listed, False otherwise
        """
        try:
            # use this for finding file names in wget .listing files, or any other remote file listing
            sh.grep("-E", '(^|\s)'+filename+'($|\s)',
                    "-s", self.remote_file_list)
            return True
        except Exception:
            return False

    def validate_mapping(self):
        # read the file (pandas may handle excel for us)
        # remove any double quotes
        self._strip_quotes()
        # remove empty rows
        self.df = self.df.dropna(how='all')
        # remove empty cols
        self._remove_empty_cols()
        # FIXME: verify that the number of headers matches the number of columns in the data (no extra tabs, etc.)
        # NOTE: ^ if there are more columns of data than header rows,
        # the data in the last column(s) are simply not read in
        # check the headers
        self._validate_headers()
        if self.errs:
            return self._get_validation_result()

        # validate the columns
        # validate SampleID column
        # return immediately if empty IDs)
        self._validate_sample_ids()
        if "null_values" in self.errs:
            return self._get_validation_result()

        # validate all columns from TreatmentGroup to Description
        # as categorical data
        metadata_cols = self._get_metadata_columns()
        if self.map_type not in [MapType.SE_WGS, MapType.PE_WGS, MapType.PE_QC, MapType.SE_QC]:
            for col in metadata_cols:
                self._is_valid_categorical_data(col)

        # FIXME: I don't like how this is working atm
        # validate any read columns - get the column names that we need
        # to check from the required headers list
        regex = re.compile(".*FastqFile")
        for col_name in [m.group(0) for col in self.df.columns.tolist() for m in [regex.search(col)] if m]:
            self._validate_read_names(col_name)
            if self.ignore_fq_files is False:
                self.check_files_exist(col_name)

        # if paired end, verify that read names exist in only one column
        if "ForwardFastqFile" in self.df.columns and "ReverseFastqFile" in self.df.columns:
            # FIXME: I don't like having to run this twice, and this shouldn't be here in the code
            indexes = numpy.where(self.df['ForwardFastqFile'].isin(
                self.df['ReverseFastqFile']) == True)
            bad_forward_names = self.df.iloc[indexes].index.values.tolist()
            if len(bad_forward_names) > 0:
                # these are only the indexes in the forward read column
                self._add_err("duplicate_forward_file_name",
                              "Filenames cannot appear in both the Forward and Reverse Fastq columns.",
                              self.df.columns.get_loc("ForwardFastqFile"), bad_forward_names)
                # so we need to fetch the indexes in the reverse read column
                indexes = numpy.where(self.df['ReverseFastqFile'].isin(
                    self.df['ForwardFastqFile']) == True)
                bad_rev_names = self.df.iloc[indexes].index.values.tolist()
                self._add_err("duplicate_reverse_file_name",
                              "Filenames cannot appear in both the Forward and Reverse Fastq columns.",
                              self.df.columns.get_loc("ReverseFastqFile"), bad_rev_names)
        return self._get_validation_result()
