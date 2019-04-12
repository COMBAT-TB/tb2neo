"""
Interface to Ensembl ftp
"""
import sys
from ftplib import FTP, all_errors
from gzip import GzipFile
from io import BytesIO

from tb2neo.ncbi import get_fasta

FILE_PATH = '/pub/bacteria/release-39/fasta/bacteria_0_collection/mycobacterium_tuberculosis_h37rv/dna'
FASTA_FILE = "Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna_sm.chromosome.Chromosome.fa.gz"
URL = 'ftp.ensemblgenomes.org'


def download_fasta(url, file_path, fasta_file):
    """
    Download fasta file from Ensembl
    :return: zipped_fasta
    """
    try:
        data_writer = BytesIO()
        ftp_conn = FTP(url, 'anonymous', 'anonymouspass')
        ftp_conn.cwd(file_path)
        # size = ftp_conn.size(fasta_file)
        ftp_conn.retrbinary("RETR {}".format(fasta_file), data_writer.write)
        ftp_conn.quit()
        data_writer.seek(0)
        zipped_file = GzipFile(fileobj=data_writer)
    except all_errors as e:
        raise e
    return zipped_file


def get_nucleotides(strain):
    """
    Gets zipped fasta file and returns seq string
    :return:
    """
    try:
        nucleotides = get_fasta(strain=strain)
        sys.stdout.write(nucleotides[:4])
    except IOError as e:
        raise e
        # fasta_zip = download_fasta(url=URL, file_path=FILE_PATH, fasta_file=FASTA_FILE)
        # seq_rec = Bio.SeqIO.read(fasta_zip, "fasta")
        # nucleotides = str(seq_rec.seq)
    return nucleotides
