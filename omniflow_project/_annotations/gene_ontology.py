import re
import requests


class Ontology:
    """
    class that stores some functionalities to connect phenotypes to nodes
    and to associate information for each node at tissue level

    """
    def __init__(self):
        gene_ontology_url = "https://golr-aux.geneontology.io/solr/select?defType=edismax&qt=standard&indent=on&wt=csv&rows=100000&start=0&fl=bioentity_label&facet=true&facet.mincount=1&facet.sort=count&json.nl=arrarr&facet.limit=25&hl=true&hl.simple.pre=%3Cem%20class=%22hilite%22%3E&hl.snippets=1000&csv.encapsulator=&csv.separator=%09&csv.header=false&csv.mv.separator=%7C&fq=document_category:%22annotation%22&fq=isa_partof_closure:%22GO:0007049%22&fq=taxon_subset_closure_label:%22Homo%20sapiens%22&fq=type:%22protein%22&fq=annotation_class_label:%22G1/S%20transition%20of%20mitotic%20cell%20cycle%22&facet.field=aspect&facet.field=taxon_subset_closure_label&facet.field=type&facet.field=evidence_subset_closure_label&facet.field=regulates_closure_label&facet.field=isa_partof_closure_label&facet.field=annotation_class_label&facet.field=qualifier&facet.field=annotation_extension_class_closure_label&facet.field=assigned_by&facet.field=panther_family_label&q=*:*"
        codes_dict = {"GO:0001837": "epithelial to mesenchymal transition"}
        return

    def modify_url_ontology(self, url, new_go_code, new_description):
        """
        Modifies a given URL by replacing a GO code and a descriptive string
        in very specific locations identified by prefixes.

        Parameters:
        - url (str): The original URL to be modified.
        - new_go_code (str): The new GO code to insert into the URL.
        - new_description (str): The new descriptive string to insert into the URL.

        Returns:
        - str: The modified URL.
        """
        # Define the prefixes that identify where the replacements should occur,
        go_code_prefix = "isa_partof_closure:%22"
        description_prefix = "annotation_class_label:%22"

        # Patterns to find the exact locations for replacements
        go_code_pattern = re.escape(go_code_prefix) + r"GO:\d{7}"
        description_pattern = re.escape(description_prefix) + r".+?%22"

        # Perform the replacements
        url = re.sub(go_code_pattern, go_code_prefix + new_go_code, url, 1)
        # URL-encode the new description
        new_description_encoded = re.sub(r" ", "%20", new_description)
        url = re.sub(description_pattern, description_prefix + new_description_encoded + "%22", url, 1)

        return url

    def fetch_nodes_from_url(self, url):

        """
        fetch the nodes in a list from the given geneontology url
        """

        response = requests.get(url)
        genes = response.text.split("\n")
        genes_unique = list(set(genes))
        return genes_unique
