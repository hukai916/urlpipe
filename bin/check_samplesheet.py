#!/usr/bin/env python


"""Provide a command line tool to validate and transform tabular samplesheets."""


import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path

logger = logging.getLogger()


class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a previously
            validated and transformed row. The order of rows is maintained.

    """

    VALID_FORMATS = (
        ".fq.gz",
        ".fastq.gz",
    )

    def __init__(
        self,
        sample          = "sample",
        fastq_1         = "fastq_1",
        fastq_2         = "fastq_2",
        single_col      = "single_end",
        allele_number   = "allele_number",
        start_allele_1  = "start_allele_1",
        end_allele_1    = "end_allele_1",
        start_allele_2  = "start_allele_2",
        end_allele_2    = "end_allele_2",
        **kwargs,
    ):
        """
        Initialize the row checker with the expected column names.

        Args:
            sample (str): the name of the column that contains the sample name.
            fastq_1 (str): the name of the column that contains the first (or only) FASTQ file path.
            fastq_2 (str): the name of the column that contains the second (if any).
            start_allele_1 (int): the lower cutoff of the first allele repeat length cutoff.
            end_allele_1 (int): the higher cutoff of the first allele repeat length cutoff.
            start_allele_2 (int): the lower cutoff of the second (if any) allele repeat length cutoff.
            end_allele_2 (int): the higher cutoff of the second allele (if any) repeat length cutoff.
            allele_number (int [1|2]): if 1, 2 repeat length cutoffs are needed; if 2, 4 cutoffs are needed.
            single_col (str): The name of the new column that will be inserted and
                records whether the sample contains single- or paired-end sequencing
                reads (default "single_end"). Not used.

        """
        super().__init__(**kwargs)
        self._sample = sample
        self._fastq_1 = fastq_1
        self._fastq_2 = fastq_2
        self._single_col = single_col
        self._start_allele_1 = start_allele_1
        self._end_allele_1 = end_allele_1
        self._start_allele_2 = start_allele_2
        self._end_allele_2 = end_allele_2
        self._allele_number = allele_number
        self._seen = set()
        self.modified = []

    def validate_and_transform(self, row):
        """
        Perform all validations on the given row and insert the read pairing status.

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).

        """
        self._validate_sample(row)
        self._validate_first(row)
        self._validate_second(row)
        self._validate_pair(row)
        if self._allele_number == 1:
            self._validate_length_cutoff_1(row)
        elif self._allele_number == 2:
            self._validate_length_cutoff_1(row)
            self._validate_length_cutoff_2(row)
        self._seen.add((row[self._sample], row[self._fastq_1]))
        self.modified.append(row)

    def _validate_sample(self, row):
        """Assert that the sample name exists and convert spaces to underscores."""
        if len(row[self._sample]) <= 0:
            raise AssertionError("Sample input is required.")
        # Sanitize samples slightly.
        row[self._sample] = row[self._sample].replace(" ", "_")

    def _validate_first(self, row):
        """Assert that the first FASTQ entry is non-empty and has the right format."""
        if len(row[self._fastq_1]) <= 0:
            raise AssertionError("At least the first FASTQ file is required.")
        self._validate_fastq_format(row[self._fastq_1])

    def _validate_second(self, row):
        """Assert that the second FASTQ entry has the right format if it exists."""
        if len(row[self._fastq_2]) > 0:
            self._validate_fastq_format(row[self._fastq_2])

    def _validate_pair(self, row):
        """Assert that read pairs have the same file extension. Report pair status."""
        if row[self._fastq_1] and row[self._fastq_2]:
            row[self._single_col] = False
            fastq_1_suffix = Path(row[self._fastq_1]).suffixes[-2:]
            fastq_2_suffix = Path(row[self._fastq_2]).suffixes[-2:]
            if fastq_1_suffix != fastq_2_suffix:
                raise AssertionError("FASTQ pairs must have the same file extensions.")
        else:
            row[self._single_col] = True

    def _validate_fastq_format(self, filename):
        """Assert that a given filename has one of the expected FASTQ extensions."""
        if not any(filename.endswith(extension) for extension in self.VALID_FORMATS):
            raise AssertionError(
                f"The FASTQ file has an unrecognized extension: {filename}\n"
                f"It should be one of: {', '.join(self.VALID_FORMATS)}"
            )

    def _validate_length_cutoff_1(self, row):
        """Assert that the length_cutoff_1 are set and are positive integers."""
        if len(row[self._start_allele_1]) <= 0:
            raise AssertionError("start_allele_1 must be set.")
        if len(row[self._end_allele_1]) <= 0:
            raise AssertionError("end_allele_1 must be set.")

        try:
            row[self._start_allele_1] = int(row[self._start_allele_1])
            row[self._end_allele_1] = int(row[self._end_allele_1])
        except ValueError as error:
            logger.critical("length_cutoff_1 must be integer.")
            sys.exit(1)

        if row[self._start_allele_1] < 0 or row[self._end_allele_1] < 0 or row[self._start_allele_1] >= row[self._end_allele_1]:
            raise AssertionError("start_allele_1 must be lower than end_allele_1, and both should be positive integer.")

    def _validate_length_cutoff_2(self, row):
        """Assert that the length_cutoff_2 are set and are positive integers."""
        if len(row[self._start_allele_2]) <= 0:
            raise AssertionError("start_allele_2 must be set.")
        if len(row[self._end_allele_2]) <= 0:
            raise AssertionError("end_allele_2 must be set.")
        try:
            # row[self._start_allele_1] = int(row[self._start_allele_1])
            # row[self._end_allele_1] = int(row[self._end_allele_1])
            row[self._start_allele_2] = int(row[self._start_allele_2])
            row[self._end_allele_2] = int(row[self._end_allele_2])
        except ValueError as error:
            logger.critical("length_cutoff_2 must be integer.")
            sys.exit(1)

        if row[self._start_allele_2] < 0 or row[self._end_allele_2] < 0 or row[self._start_allele_2] >= row[self._end_allele_2]:
            raise AssertionError("start_allele_2 must be lower than end_allele_2, and both should be positive integer.")
        if row[self._start_allele_2] < row[self._end_allele_1]:
            raise AssertionError("The following must hold: start_allele_1 < end_allele_1 < start_allele_2 < end_allele_2.")

    def validate_unique_samples(self):
        """
        Assert that the combination of sample name and FASTQ filename is unique.

        In addition to the validation, also rename all samples to have a suffix of _T{n}, where n is the
        number of times the same sample exist, but with different FASTQ files, e.g., multiple runs per experiment.

        """
        if len(self._seen) != len(self.modified):
            raise AssertionError("The pair of sample name and FASTQ must be unique.")
        seen = Counter()
        for row in self.modified:
            sample = row[self._sample]
            seen[sample] += 1
            # row[self._sample] = f"{sample}_T{seen[sample]}"
            row[self._sample] = f"{sample}"


def read_head(handle, num_lines=10):
    """Read the specified number of lines from the current position in the file."""
    lines = []
    for idx, line in enumerate(handle):
        if idx == num_lines:
            break
        lines.append(line)
    return "".join(lines)


def sniff_format(handle):
    """
    Detect the tabular format.

    Args:
        handle (text file): A handle to a `text file`_ object. The read position is
        expected to be at the beginning (index 0).

    Returns:
        csv.Dialect: The detected tabular format.

    .. _text file:
        https://docs.python.org/3/glossary.html#term-text-file

    """
    peek = read_head(handle)
    handle.seek(0)
    sniffer = csv.Sniffer()
    if not sniffer.has_header(peek):
        logger.critical("The given sample sheet does not appear to contain a header.")
        sys.exit(1)
    dialect = sniffer.sniff(peek)
    return dialect


def check_samplesheet(file_in, file_out, allele_number):
    """
    Check that the tabular samplesheet has the structure expected by nf-core pipelines.

    Validate the general shape of the table, expected columns, and each row. Also add
    an additional column which records whether one or two FASTQ reads were found.

    Args:
        file_in (pathlib.Path): The given tabular samplesheet. The format can be either
            CSV, TSV, or any other format automatically recognized by ``csv.Sniffer``.
        file_out (pathlib.Path): Where the validated and transformed samplesheet should
            be created; always in CSV format.
        allele_number: determines how many length_cutoffs are expected

    Example:
        This function checks that the samplesheet follows the following structure,
        see also the `viral recon samplesheet`_::

            sample,fastq_1,fastq_2,start_allele_1,end_allele_1,start_allele_2,end_allele_2
            SAMPLE_PE,SAMPLE_PE_RUN1_1.fastq.gz,SAMPLE_PE_RUN1_2.fastq.gz
            SAMPLE_PE,SAMPLE_PE_RUN2_1.fastq.gz,SAMPLE_PE_RUN2_2.fastq.gz
            SAMPLE_SE,SAMPLE_SE_RUN1_1.fastq.gz,

    .. _viral recon samplesheet:
        https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv

    """
    if allele_number == 1:
        required_columns = {"sample", "fastq_1", "fastq_2", "start_allele_1", "end_allele_1"}
    elif allele_number == 2:
        required_columns = {"sample", "fastq_1", "fastq_2", "start_allele_1", "end_allele_1", "start_allele_2", "end_allele_2"}
    else:
        logger.critical(f"The params.allele_number must be integer 1 or 2.")
        sys.exit(1)
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_in.open(newline="") as in_handle:
        reader = csv.DictReader(in_handle)
        # Validate the existence of the expected header columns.
        if not required_columns.issubset(reader.fieldnames):
            print(reader.fieldnames)
            req_cols = ", ".join(required_columns)
            logger.critical(f"The sample sheet **must** contain these column headers: {req_cols}.")
            sys.exit(1)
        # Validate each row.
        checker = RowChecker(allele_number = allele_number)
        for i, row in enumerate(reader):
            try:
                checker.validate_and_transform(row)
            except AssertionError as error:
                logger.critical(f"{str(error)} On line {i + 2}.")
                sys.exit(1)
        checker.validate_unique_samples()
    header = list(reader.fieldnames)
    header.insert(1, "single_end")
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_out.open(mode="w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, header, delimiter=",")
        writer.writeheader()
        for row in checker.modified:
            writer.writerow(row)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular input samplesheet in CSV or TSV format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Transformed output samplesheet in CSV format.",
    )
    parser.add_argument(
        "allele_number",
        metavar="ALLELE_NUMBER",
        type=int,
        help="Allele number in integer.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_samplesheet(args.file_in, args.file_out, args.allele_number)


if __name__ == "__main__":
    sys.exit(main())
