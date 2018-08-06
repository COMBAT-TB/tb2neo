"""
Variant Processing
"""
import re
import sys

from Bio.SeqUtils import seq3
from bioservices import KEGG

from gff2neo.dbconn import create_known_mutation_nodes

kegg = KEGG(verbose=False)


def get_drug_info(drug_name):
    """
    Get Drug Info
    :param drug_name:
    :return:
    """
    drug_groups = ["aminoglycosides", "fluoroquinolones"]
    drugbank_dict, drugbank_id, drug_info = dict(), None, dict()
    if 'aminosalicylic_acid' in drug_name:
        drug_name = 'Aminosalicylic acid'
    if drug_name not in drugbank_dict.values() or drug_groups:
        k_drug = kegg.find("drug", str(drug_name))
        drug_id = k_drug.split('\t')[0]
        drug_info = kegg.parse(kegg.get(drug_id))
        if not isinstance(drug_info, int):
            dblinks = drug_info.get("DBLINKS", None)
            if dblinks:
                drugbank_id = dblinks.get("DrugBank", None)
                if drugbank_id:
                    drugbank_dict[drugbank_id] = drug_name
        else:
            for drug_id, name in drugbank_dict.items():
                if drug_name == name:
                    drugbank_id = drug_id

    return drugbank_id


def _process_coll_mutations(in_file, cset_name):
    """
    Load Coll et al. mutation Library
    :param in_file:
    :param cset_name:
    :return:
    """
    sys.stdout.write("\nProcessing {cset}...\n".format(cset=cset_name))
    drug_groups = ["aminoglycosides", "fluoroquinolones"]
    drugbank_id, biotype = None, ''
    vset_name = "doi.org/10.1186/s13073-015-0164-0"
    vset_owner = "Coll, F. et al."
    for line in in_file:
        tab_split = line.strip().split('\t')
        drug_name = tab_split[9].lower()
        drugbank_id = get_drug_info(drug_name)

        variant_pos = tab_split[1]
        # gi|444893469|emb|AL123456.3|	1674484/1674485	inhA	Rv1484	283/284	95	A/T/C/C	ATT/CCT	I/P	ISONIAZID	MUBII-TB-DB
        # check base pair length
        # ref_allele, alt_allele = "", ""
        bps = tab_split[6]
        bp_mid = len(bps.split("/")) / 2
        bps = bps.replace("/", "")
        ref_allele = bps[:bp_mid]
        alt_allele = bps[bp_mid:]

        # amino acid change
        gene_cord = tab_split[4]
        codon_number = tab_split[5]
        # 'L/P' to ['Leu', 'Pro']
        amino_change = [seq3(a, custom_map={"*": "Stop"}, undef_code='-') for a in tab_split[8].strip().split("/")]
        if len(amino_change) > 1 and amino_change[0] == amino_change[1]:
            biotype = 'synonymous'
        elif len(amino_change) > 1 and amino_change[0] is not amino_change[1]:
            biotype = 'non-synonymous'
        elif len(ref_allele) is not len(alt_allele):
            biotype = "indel"

        # ['Leu', 'Pro'] to ['Leu', 123,'Pro']
        amino_change.insert(1, codon_number)
        # ['Leu', 123, 'Pro'] to 'Leu123Pro'
        consequence = ''.join(amino_change) if biotype is not "indel" else ref_allele + gene_cord + alt_allele
        # some string manipulation kung-fu
        sources = tab_split[10].translate(None, "()").replace('"', '')

        if "_promoter" in tab_split[2]:
            promoter = tab_split[2]
            gene_name = tab_split[2].split("_")[0]
            consequence = ref_allele + gene_cord + alt_allele
            biotype = 'promoter'
        else:
            gene_name = tab_split[2]
            promoter = None

        create_known_mutation_nodes(chrom="Chr1", pos=variant_pos, ref_allele=str(ref_allele),
                                    alt_allele=str(alt_allele), gene=gene_name, promoter=promoter,
                                    pk=str(drug_name + consequence +
                                           variant_pos + ref_allele + alt_allele).lower(),
                                    consequence=consequence, vset_name=vset_name, vset_owner=vset_owner,
                                    cset_name=cset_name, sources=sources,
                                    drugbank_id=drugbank_id, drug_name=drug_name, biotype=biotype)


def _process_tbprofiler_mutations(in_file, cset_name):
    """
    Load TBProfiler mutation Library
    :param in_file:
    :param cset_name:
    :return:
    """
    drugbank_dict, drugbank_id, biotype = dict(), None, ''
    vset_name = "https://github.com/jodyphelan/TBProfiler"
    vset_owner = "Coll, F. et al."
    for line in in_file:
        tab_split = line.split('\t')
        drug_name = tab_split[0].lower()
        drugbank_id = get_drug_info(drug_name)

        variant_pos = tab_split[1]
        ref_allele = tab_split[2]
        alt_allele = tab_split[3]
        if "_promoter" in tab_split[4]:
            promoter = tab_split[4]
            gene_name = tab_split[4].split("_")[0]
            biotype = 'promoter'
        else:
            gene_name = tab_split[4]
            promoter = None
        # amino acid change
        consequence = tab_split[5].strip()
        # Pro241Pro to ['Pro', '241', 'Pro']
        amino_change = re.split('(\d+)', consequence)
        if amino_change[0] == amino_change[2] and consequence.isalnum() and any(c.islower() for c in consequence):
            biotype = 'synonymous'
        elif amino_change[0] is not amino_change[2] and consequence.isalnum() and any(c.islower() for c in consequence):
            biotype = 'non-synonymous'
        elif consequence.isupper() and not any(c.islower() for c in consequence) and '-' not in consequence:
            biotype = "indel"

        create_known_mutation_nodes(chrom="Chr1", pos=variant_pos, ref_allele=str(ref_allele),
                                    alt_allele=str(alt_allele), gene=gene_name, promoter=promoter,
                                    pk=str(drug_name + consequence +
                                           variant_pos + ref_allele + alt_allele).lower(),
                                    consequence=consequence, vset_name=vset_name, vset_owner=vset_owner,
                                    cset_name=cset_name,
                                    drugbank_id=drugbank_id, drug_name=drug_name, biotype=biotype)


def process_mutation_file(in_file):
    """
    Process mutation file
    :param in_file:
    :return:
    """

    if in_file and in_file.endswith(".txt"):
        with open(in_file) as in_file:
            next(in_file)
            cset_name = str(in_file.name).split('/')[-1]
            if 'coll' in cset_name:
                # pass
                _process_coll_mutations(in_file=in_file, cset_name=cset_name)
            elif 'drdb' in cset_name:
                # pass
                _process_tbprofiler_mutations(in_file=in_file, cset_name=cset_name)
            elif 'phyresse' in cset_name:
                pass
                # _process_phyresse_mutations(in_file=in_file, cset_name=cset_name)
            elif 'tgstb' in cset_name:
                pass
                # _process_tgstb_mutations(in_file=in_file, cset_name=cset_name)
