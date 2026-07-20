# CHDI-Style Gene Dossier API Validation Report - SREBF2

Prepared from the Postman validation workflow for `SREBF2`, with `HTT` intended as a later Huntington disease stress-test gene and `MSH3`/modifier genes intended for later expansion.

This report lists the correct API call formats, the data each call retrieves, the fields to extract, chained identifiers, gotchas, and provenance notes. API keys are shown only as Postman variables.

## Global Identifiers And Variables

Validated SREBF2 identifiers:


| Field             | Human                | Mouse                | Rat                  |
| ----------------- | -------------------- | -------------------- | -------------------- |
| Entrez Gene ID    | `6721`               | `20788`              | `300095`             |
| Gene symbol       | `SREBF2`             | `Srebf2`             | `Srebf2`             |
| Ensembl gene ID   | `ENSG00000198911`    | `ENSMUSG00000022463` | `ENSRNOG00000007400` |
| UniProt accession | `Q12772`             | `Q3U1N2`             | `Q3T1I5`             |
| RefSeq protein    | `NP_004590.2`        | use source response  | use source response  |
| GTEx GENCODE ID   | `ENSG00000198911.11` | n/a                  | n/a                  |
| Mouse MGI ID      | n/a                  | `MGI:107585`         | n/a                  |


Recommended Postman variables:

```text
ncbi_api_key
omim_api_key
biogrid_accesskey
serpapi_api_key
```

Optional convenience variables:

```text
ncbi_base_url = https://eutils.ncbi.nlm.nih.gov/entrez/eutils/
ensembl_base_url = https://rest.ensembl.org
uniprot_base_url = https://rest.uniprot.org
ucsc_api_base_url = https://api.genome.ucsc.edu
gtex_base_url = https://gtexportal.org/api/v2
allen_api_base_url = https://api.brain-map.org/api/v2
pdbe_base_url = https://www.ebi.ac.uk/pdbe/api
alphafold_base_url = https://alphafold.ebi.ac.uk/api
ncbi_datasets_base_url = https://api.ncbi.nlm.nih.gov/datasets/v2
```

NOTE: Keep API keys in environment variables only. Do not export collections with live key values.

---



## Section 1a - General Gene Information And Aliases

Goal: build the Human/Mouse/Rat table with Entrez Gene ID, symbol, name, Ensembl ID, UniProt ID, and aliases.

### 1a.1 NCBI Gene ESearch

Purpose: find the Entrez Gene ID candidates for a gene symbol and organism.

Request format:

```text
GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=SREBF2%5BGene%20Name%5D%20AND%20Homo%20sapiens%5BOrganism%5D&retmode=json&sort=relevance&api_key={{ncbi_api_key}}
```

Mouse example:

```text
GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=SREBF2%5BGene%20Name%5D%20AND%20Mus%20musculus%5BOrganism%5D&retmode=json&sort=relevance&api_key={{ncbi_api_key}}
```

Rat example:

```text
GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=SREBF2%5BGene%20Name%5D%20AND%20Rattus%20norvegicus%5BOrganism%5D&retmode=json&sort=relevance&api_key={{ncbi_api_key}}
```

Data retrieved: candidate Entrez Gene IDs.

Fields to extract:

```text
esearchresult.count
esearchresult.idlist
esearchresult.querytranslation
```

Chained identifiers:

```text
gene_id_candidates = esearchresult.idlist
```

Is it what we need? Yes, for candidate ID discovery.

NOTE: Rat returned both `300095` and retired/interim `404651`. Use the current official record, `300095`.

### 1a.2 NCBI Gene ESummary

Purpose: retrieve official symbol, name, aliases, genomic location, status, and organism metadata.

Request format:

```text
GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=6721&retmode=json&api_key={{ncbi_api_key}}
```

Mouse:

```text
GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=20788&retmode=json&api_key={{ncbi_api_key}}
```

Rat candidates:

```text
GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=300095,404651&retmode=json&api_key={{ncbi_api_key}}
```

Fields to extract:

```text
result.uids
result.{gene_id}.uid
result.{gene_id}.name
result.{gene_id}.description
result.{gene_id}.otheraliases
result.{gene_id}.otherdesignations
result.{gene_id}.nomenclaturesymbol
result.{gene_id}.nomenclaturename
result.{gene_id}.nomenclaturestatus
result.{gene_id}.status
result.{gene_id}.currentid
result.{gene_id}.organism.scientificname
result.{gene_id}.organism.taxid
result.{gene_id}.genomicinfo[*].chraccver
result.{gene_id}.genomicinfo[*].chrstart
result.{gene_id}.genomicinfo[*].chrstop
```

Selection rule:

```text
taxid matches expected species
nomenclaturestatus == Official
currentid is empty
status is empty or not retired
symbol matches expected species casing
```

Is it what we need? Yes, for Entrez ID, symbol, gene name, aliases, and current-record validation.

NOTE: `currentid` is the key field for avoiding retired NCBI records.

### 1a.3 Ensembl Gene Lookup

Purpose: retrieve Ensembl gene ID and genomic coordinates by species and symbol.

Request format:

```text
GET https://rest.ensembl.org/lookup/symbol/homo_sapiens/SREBF2?content-type=application/json
```

Mouse:

```text
GET https://rest.ensembl.org/lookup/symbol/mus_musculus/Srebf2?content-type=application/json
```

Rat:

```text
GET https://rest.ensembl.org/lookup/symbol/rattus_norvegicus/Srebf2?content-type=application/json
```

Recommended headers:

```text
Accept: application/json
Content-Type: application/json
```

Fields to extract:

```text
id
display_name
description
species
biotype
canonical_transcript
assembly_name
seq_region_name
start
end
strand
version
```

Chained identifiers:

```text
ensembl_gene_id = id
canonical_transcript = canonical_transcript
```

Is it what we need? Yes, for Ensembl ID.

NOTE: The typo `aplication/json` causes Ensembl to return YAML/HTML-like output. Use `application/json`.

### 1a.4 UniProtKB Search

Purpose: retrieve reviewed Swiss-Prot UniProt accession for each species.

Human request:

```text
GET https://rest.uniprot.org/uniprotkb/search?query=%28gene_exact:SREBF2%29%20AND%20%28organism_id:9606%29%20AND%20%28reviewed:true%29&format=json&fields=accession,id,gene_names,protein_name,organism_name,organism_id,xref_ensembl
```

Mouse:

```text
GET https://rest.uniprot.org/uniprotkb/search?query=%28gene_exact:Srebf2%29%20AND%20%28organism_id:10090%29%20AND%20%28reviewed:true%29&format=json&fields=accession,id,gene_names,protein_name,organism_name,organism_id,xref_ensembl
```

Rat:

```text
GET https://rest.uniprot.org/uniprotkb/search?query=%28gene_exact:Srebf2%29%20AND%20%28organism_id:10116%29%20AND%20%28reviewed:true%29&format=json&fields=accession,id,gene_names,protein_name,organism_name,organism_id,xref_ensembl
```

Fields to extract:

```text
results[0].primaryAccession
results[0].uniProtkbId
results[0].genes[*].geneName.value
results[0].proteinDescription.recommendedName.fullName.value
results[0].organism.scientificName
results[0].organism.taxonId
results[0].uniProtKBCrossReferences[database="Ensembl"].id
```

Chained identifiers:

```text
uniprot_accession = results[0].primaryAccession
```

Is it what we need? Yes, for UniProt ID.

NOTE: Use `reviewed:true` to avoid unreviewed isoform noise.

---



## Section 1b - UCSC Conservation / Genome Browser Context

Goal: identify the transcript/region and create the conservation/browser image.

### 1b.1 UCSC Search

Purpose: find UCSC hits and coordinates for SREBF2.

```text
GET https://api.genome.ucsc.edu/search?genome=hg38&search=SREBF2
```

Fields to extract:

```text
positionMatches[*].db
positionMatches[*].name
positionMatches[*].description
positionMatches[*].position
```

Chained identifiers:

```text
ucsc_chrom
ucsc_start
ucsc_end
canonical_transcript
```

Validated SREBF2 region:

```text
chr22:41833105-41907305
canonical transcript: ENST00000361204.9
```



### 1b.2 UCSC Track Data

Purpose: retrieve transcript/track records for the selected region.

```text
GET https://api.genome.ucsc.edu/getData/track?genome=hg38&track=knownGene&chrom=chr22&start=41833105&end=41907305
```

Fields to extract:

```text
knownGene[*].name
knownGene[*].chrom
knownGene[*].txStart
knownGene[*].txEnd
knownGene[*].strand
knownGene[*].exonStarts
knownGene[*].exonEnds
knownGene[*].geneName
knownGene[*].geneName2
```

Is it what we need? Yes, for track coordinates and transcript confirmation.

NOTE: UCSC requires both `start` and `end`; omitting `start` causes bad request.

### 1b.3 UCSC Browser Screenshot URL

Purpose: browser image for report.

```text
GET https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr22:41833105-41907305&knownGene=pack&cons100way=full
```

Is it what we need? Yes, for the image/screenshot.

NOTE: The `hgTracks` URL is a browser view, not JSON. Capture the image manually or with a browser automation step later.

---



## Section 1c - Known Structure: CDD And PDBe



### 1c.1 NCBI CDD Batch CD-Search Submit

Purpose: submit SREBF2 RefSeq protein to conserved domain search.

```text
GET https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi?queries=NP_004590.2&db=cdd&smode=auto&useid1=true&maxhit=250&filter=true&evalue=0.01
```

Fields to extract:

```text
cdsid
status
```

Chained identifier:

```text
cdd_cdsid
```

Is it what we need? Yes, but only submits the job.

NOTE: Batch CDD returns text/status, not normal JSON. Poll until `status=0` or completed.

### 1c.2 NCBI CDD Retrieve Domain Hits

Purpose: retrieve conserved domain hits.

```text
GET https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi?cdsid={{cdd_cdsid}}&tdata=hits&dmode=full&qdefl=true&cddefl=true
```

Fields to extract from completed hit output:

```text
query accession
domain accession
domain short name
domain description
from residue
to residue
evalue
bitscore
superfamily/accession if present
```

Validated important SREBF2 domains:

```text
cd18922 - bHLHzip_SREBP2 - residues 326-402
cl00081 - bHLH_SF - residues 326-402
```



### 1c.3 NCBI CDD Alignment JSON

Purpose: retrieve residue-level alignments.

```text
GET https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi?cdsid={{cdd_cdsid}}&tdata=aligns&alnfmt=json
```

Is it what we need? Optional. Use for domain alignment details, not the main table.

NOTE: `alnfmt=json` only applies when `tdata=aligns`.

### 1c.4 PDBe Best Structures

Purpose: find PDB structures mapped to SREBF2 UniProt accession.

```text
GET https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/Q12772
```

Fields to extract:

```text
Q12772[*].pdb_id
Q12772[*].chain_id
Q12772[*].unp_start
Q12772[*].unp_end
Q12772[*].pdb_start
Q12772[*].pdb_end
Q12772[*].coverage
Q12772[*].resolution
Q12772[*].experimental_method
```

Chained identifier:

```text
pdb_id
```

Validated SREBF2 PDB:

```text
1ukl
```



### 1c.5 PDBe UniProt Mapping

Purpose: validate SIFTS chain mappings for selected PDB.

```text
GET https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/1ukl
```

Fields to extract:

```text
1ukl.UniProt.Q12772.name
1ukl.UniProt.Q12772.mappings[*].chain_id
1ukl.UniProt.Q12772.mappings[*].start.residue_number
1ukl.UniProt.Q12772.mappings[*].end.residue_number
1ukl.UniProt.Q12772.mappings[*].unp_start
1ukl.UniProt.Q12772.mappings[*].unp_end
```

NOTE: If you query a PDB that does not map to Q12772, PDBe may return "Requested endpoint does not contain any data."

### 1c.6 PDBe Entry Summary

Purpose: retrieve title, authors, method, release dates.

```text
GET https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/1ukl
```

Fields to extract:

```text
1ukl[0].title
1ukl[0].experimental_method
1ukl[0].experimental_method_class
1ukl[0].entry_authors
1ukl[0].deposition_date
1ukl[0].release_date
1ukl[0].revision_date
```

Is it what we need? Yes, for structure title and provenance.

---



## Section 1d - AlphaFold Protein Structure Prediction

Purpose: retrieve predicted protein model metadata and links for human, mouse, rat.

Human:

```text
GET https://alphafold.ebi.ac.uk/api/prediction/Q12772
```

Mouse:

```text
GET https://alphafold.ebi.ac.uk/api/prediction/Q3U1N2
```

Rat:

```text
GET https://alphafold.ebi.ac.uk/api/prediction/Q3T1I5
```

Fields to extract:

```text
[0].entryId
[0].modelEntityId
[0].uniprotAccession
[0].uniprotId
[0].uniprotDescription
[0].organismScientificName
[0].taxId
[0].globalMetricValue
[0].fractionPlddtVeryLow
[0].fractionPlddtLow
[0].fractionPlddtConfident
[0].fractionPlddtVeryHigh
[0].latestVersion
[0].modelCreatedDate
[0].pdbUrl
[0].cifUrl
[0].bcifUrl
[0].paeImageUrl
[0].plddtDocUrl
[0].paeDocUrl
```

Build viewer link:

```text
https://alphafold.ebi.ac.uk/entry/{uniprot_accession}
```

NOTE: The browser viewer link is not returned directly in the JSON. Build it from the accession.

---



## Section 1e - Homologues In Model Animals

Purpose: retrieve NCBI ortholog reports.

```text
GET https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/6721/orthologs?returned_content=COMPLETE&page_size=1000
```

Fields to extract:

```text
reports[*].gene.gene_id
reports[*].gene.symbol
reports[*].gene.description
reports[*].gene.tax_id
reports[*].gene.common_name
reports[*].gene.organism.scientific_name
reports[*].gene.chromosomes
reports[*].gene.genomic_ranges
reports[*].gene.annotations
```

Is it what we need? Yes, for homologue/ortholog table.

NOTE: Use `page_size=1000` to avoid missing entries when the default page is too small.

---



## Section 2a - Tissue Expression: GTEx



### 2a.1 GTEx Gene Lookup

```text
GET https://gtexportal.org/api/v2/reference/gene?geneId=SREBF2&genomeBuild=GRCh38/hg38
```

Fields to extract:

```text
data[*].gencodeId
data[*].geneSymbol
data[*].entrezGeneId
data[*].description
data[*].chromosome
data[*].start
data[*].end
data[*].tss
data[*].genomeBuild
```

Chained identifier:

```text
gtex_gencode_id = ENSG00000198911.11
```



### 2a.2 GTEx Median Expression

```text
GET https://gtexportal.org/api/v2/expression/medianGeneExpression?gencodeId=ENSG00000198911.11&datasetId=gtex_v8
```

Fields:

```text
data[*].tissueSiteDetailId
data[*].median
data[*].unit
```

Is it what we need? Yes, for tissue expression summary.

### 2a.3 GTEx Sample-Level Expression

```text
GET https://gtexportal.org/api/v2/expression/geneExpression?gencodeId=ENSG00000198911.11&datasetId=gtex_v8
```

Fields:

```text
data[*].tissueSiteDetailId
data[*].data
```

Is it what we need? Yes, for box/violin style plots.

NOTE: GTEx is human only. Use it for human tissue/brain expression, not mouse/rat.

---



## Section 2b - Brain Expression: Allen HBA And BrainRNASeq



### 2b.1 Allen Human Brain Atlas Probe Lookup

Purpose: find microarray probes for SREBF2.

```text
GET https://api.brain-map.org/api/v2/data/query.json?criteria=model::Probe,rma::criteria,%5Bprobe_type$eq'DNA'%5D,products%5Babbreviation$eq'HumanMA'%5D,gene%5Bacronym$eq'SREBF2'%5D,rma::options%5Bonly$eq'probes.id,probes.name,genes.acronym,genes.name,genes.entrez_id'%5D
```

Fields:

```text
msg[*].id
msg[*].name
```

Validated probe IDs:

```text
1051154
1067243
1051153
```



### 2b.2 Allen HBA Microarray Expression

Purpose: retrieve expression for one probe.

```text
GET https://api.brain-map.org/api/v2/data/query.json?criteria=service::human_microarray_expression%5Bprobes%24eq1051154%5D
```

Fields:

```text
msg.probes
msg.samples
msg.expression
```

Is it what we need? Yes, when the service returns JSON.

NOTE: This Allen service was fragile in Postman because criteria encoding is strict. It worked through correctly encoded/curl-style requests. Validate one probe at a time.

### 2b.3 BrainRNASeq Human CSV

```text
GET https://brainrnaseq.org/wp-content/uploads/2022/09/fe-wp-dataset-124.csv
```

Fields:

```text
gene_id
id
astrocytes_fetal_*
astrocytes_mature_*
endothelial_*
microglia_*
neurons_*
oligodendrocytes_*
```



### 2b.4 BrainRNASeq Mouse CSV

```text
GET https://brainrnaseq.org/wp-content/uploads/2022/09/fe-wp-dataset-120.csv
```

Fields:

```text
gene_id
id
astrocytes_*
endothelial_*
microglia_macrophage_*
myelinating_oligodendrocyte_*
neurons_*
newly_formed_oligodendrocyte_*
opc_*
```

NOTE: BrainRNASeq is CSV download, not a JSON API. Parse rows where `gene_id` contains SREBF2/Srebf2.

---



## Section 2c - Brain Cell-Type Expression

Purpose: cell-type/single-nucleus RNA-seq images and qualitative interpretation.

Useful URL:

```text
https://celltypes.brain-map.org/rnaseq/human_m1_10x?selectedVisualization=Scatter+Plot&colorByFeature=Gene+Expression&colorByFeatureValue=SREBF2
```

Is it what we need? Yes for screenshots/visual evidence, not a clean API table.

NOTE: AllenSDK `CellTypesApi` is mainly morphology/ephys/cell metadata. It does not provide the simple gene-expression UMAP data needed for this section through the documented REST workflow.

NOTE: DropViz is portal/file based, not a clean REST API for Postman. Treat as manual screenshot/source validation.

---



## Section 3 - GEO Perturbations That Alter The Gene



### 3.1 GEO Profiles ESearch

Broad mouse gene search:

```text
GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=geoprofiles&term=Srebf2%20AND%20%22Mus%20musculus%22%5BOrganism%5D&retmode=json&retmax=50&sort=relevance&api_key={{ncbi_api_key}}
```

Brain/neuron context:

```text
GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=geoprofiles&term=Srebf2%20AND%20%22Mus%20musculus%22%5BOrganism%5D%20AND%20%28brain%20OR%20neuron%20OR%20hippocampus%20OR%20cortex%20OR%20cerebellum%20OR%20striatum%29&retmode=json&retmax=50&sort=relevance&api_key={{ncbi_api_key}}
```

Perturbation context:

```text
GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=geoprofiles&term=Srebf2%20AND%20%22Mus%20musculus%22%5BOrganism%5D%20AND%20%28brain%20OR%20neuron%20OR%20hippocampus%20OR%20cortex%20OR%20cerebellum%20OR%20striatum%29%20AND%20%28stress%20OR%20fluoxetine%20OR%20antidepressant%20OR%20paraquat%20OR%20mutant%20OR%20knockout%20OR%20treatment%20OR%20exposed%20OR%20disease%29&retmode=json&retmax=50&sort=relevance&api_key={{ncbi_api_key}}
```

Fields:

```text
esearchresult.count
esearchresult.idlist
esearchresult.querytranslation
```

Chained:

```text
geo_profile_id
```



### 3.2 GEO Profiles To GDS Link

Purpose: map profile ID to GEO Dataset ID.

```text
GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=geoprofiles&db=gds&id=97740750&retmode=json&api_key={{ncbi_api_key}}
```

Fields:

```text
linksets[*].ids
linksets[*].linksetdbs[*].links
```

Chained:

```text
gds_uid
```



### 3.3 GEO DataSets ESummary

```text
GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id=4524&retmode=json&api_key={{ncbi_api_key}}
```

Validated examples:

```text
97740750 -> GDS4524
121571950 -> GDS5307
72531839 -> GDS3913
```

Fields:

```text
result.{gds_uid}.accession
result.{gds_uid}.title
result.{gds_uid}.summary
result.{gds_uid}.gpl
result.{gds_uid}.gse
result.{gds_uid}.taxon
result.{gds_uid}.gdstype
result.{gds_uid}.valtype
result.{gds_uid}.ssinfo
result.{gds_uid}.subsetinfo
result.{gds_uid}.n_samples
result.{gds_uid}.samples[*].accession
result.{gds_uid}.samples[*].title
result.{gds_uid}.pubmedids
result.{gds_uid}.ftplink
```

Is it what we need? Yes, for perturbation dataset metadata and sample labels.

NOTE: `esummary.fcgi?db=geoprofiles` failed. Correct chain is `geoprofiles ESearch -> ELink to gds -> GDS ESummary`.

---



## Section 4 - Transcription Factors / Regulators



### 4.1 Harmonizome Gene Associations

```text
GET https://maayanlab.cloud/Harmonizome/api/1.0/gene/SREBF2?showAssociations=true
```

Fields:

```text
symbol
name
synonyms
associations[*].gene.symbol
associations[*].gene.name
associations[*].dataset.name
associations[*].attribute.name
associations[*].attribute.href
associations[*].thresholdValue
associations[*].standardizedValue
```

Filter datasets:

```text
ENCODE Transcription Factor Binding Site Profiles
ENCODE Transcription Factor Targets
ChEA Transcription Factor Binding Site Profiles
ChEA Transcription Factor Targets
JASPAR Predicted Transcription Factor Targets
MotifMap Predicted Transcription Factor Targets
```



### 4.2 Harmonizome Association Detail

```text
GET https://maayanlab.cloud/Harmonizome/api/1.0/gene_set/TEAD4_HepG2_hg19_1/ENCODE+Transcription+Factor+Binding+Site+Profiles?showAssociations=true
```

Fields:

```text
attribute.name
dataset.name
associations[*].gene.symbol
associations[*].standardizedValue
```

Is it what we need? Yes, for TF association table.

NOTE: Binding-site profiles use names like `TEAD4_HepG2_hg19_1`, not just `TEAD4`.

---



## Section 5 - Protein-Protein Interaction Partners



### 5.1 STRING Identifier Mapping

```text
GET https://string-db.org/api/json/get_string_ids?identifiers=SREBF2&species=9606&echo_query=1&caller_identity=gene_dossier_postman
```

Fields:

```text
[0].queryItem
[0].stringId
[0].preferredName
[0].annotation
[0].ncbiTaxonId
```

Chained:

```text
string_id = 9606.ENSP00000354476
```



### 5.2 STRING Interaction Partners

```text
GET https://string-db.org/api/json/interaction_partners?identifiers=9606.ENSP00000354476&species=9606&limit=100&required_score=400&network_type=functional&caller_identity=gene_dossier_postman
```

Fields:

```text
[*].preferredName_A
[*].preferredName_B
[*].stringId_A
[*].stringId_B
[*].score
[*].nscore
[*].fscore
[*].pscore
[*].ascore
[*].escore
[*].dscore
[*].tscore
```



### 5.3 STRING Network Image

```text
GET https://string-db.org/api/highres_image/network?identifiers=9606.ENSP00000354476&species=9606&add_color_nodes=10&add_white_nodes=20&network_flavor=evidence&network_type=functional&hide_disconnected_nodes=1&caller_identity=gene_dossier_postman
```

Is it what we need? Yes, for the network image.

### 5.4 BioGRID Interaction Table

```text
GET https://webservice.thebiogrid.org/interactions/?searchNames=true&geneList=SREBF2&taxId=9606&includeInteractors=true&includeInteractorInteractions=false&selfInteractionsExcluded=true&interSpeciesExcluded=true&format=jsonExtended&max=10000&accesskey={{biogrid_accesskey}}
```

Fields:

```text
BIOGRID_INTERACTION_ID
ENTREZ_GENE_A
ENTREZ_GENE_B
OFFICIAL_SYMBOL_A
OFFICIAL_SYMBOL_B
SYNONYMS_A
SYNONYMS_B
EXPERIMENTAL_SYSTEM
EXPERIMENTAL_SYSTEM_TYPE
AUTHOR
PUBMED_ID
ORGANISM_A
ORGANISM_B
THROUGHPUT
SOURCE_DATABASE
```

Is it what we need? Yes, for curated PPI table and supplementary file.

NOTE: Do not call BioGRID without filters. The unfiltered endpoint returns millions of interactions.

---



## Section 6 - Chemical Perturbations Affecting The Gene: CTD



### 6.1 CTD Batch Query

```text
GET https://ctdbase.org/tools/batchQuery.go?inputType=gene&inputTerms=SREBF2&report=cgixns&actionTypes=ANY&format=tsv
```

Fields:

```text
Input
ChemicalName
ChemicalID
CasRN
GeneSymbol
GeneID
Organism
OrganismID
Interaction
InteractionActions
PubMedIDs
```

Is it what we need? Yes, for chemical-gene perturbation table and top interacting chemicals.

NOTE: CTD returns TSV, not JSON. The gene page is useful for provenance but the batch query is the structured data source.

---



## Section 7 - Chemical Tools And Effects



### 7.1 ChEMBL Target Search

```text
GET https://www.ebi.ac.uk/chembl/api/data/target/search.json?q=SREBF2
```

Fields:

```text
targets[*].target_chembl_id
targets[*].pref_name
targets[*].organism
targets[*].target_type
targets[*].target_components[*].accession
targets[*].target_components[*].target_component_synonyms
```

Chained:

```text
chembl_target_id
```



### 7.2 ChEMBL Assay Search

```text
GET https://www.ebi.ac.uk/chembl/api/data/assay.json?description__icontains=SREBP2&limit=100
```

Additional searches:

```text
GET https://www.ebi.ac.uk/chembl/api/data/assay.json?description__icontains=SREBF2&limit=100
GET https://www.ebi.ac.uk/chembl/api/data/assay.json?description__icontains=sterol%20regulatory%20element-binding%20protein&limit=100
```

Fields:

```text
assays[*].assay_chembl_id
assays[*].description
assays[*].assay_type
assays[*].assay_organism
assays[*].assay_cell_type
assays[*].target_chembl_id
assays[*].document_chembl_id
```

Chained:

```text
assay_chembl_ids
```



### 7.3 ChEMBL Activities

Correct format:

```text
GET https://www.ebi.ac.uk/chembl/api/data/activity.json?assay_chembl_id__in=CHEMBL1827133,CHEMBL4369049&limit=1000
```

Fields:

```text
activities[*].molecule_chembl_id
activities[*].canonical_smiles
activities[*].standard_type
activities[*].standard_relation
activities[*].standard_value
activities[*].standard_units
activities[*].pchembl_value
activities[*].assay_chembl_id
activities[*].document_chembl_id
```

NOTE: Remove the trailing `?` after `limit=1000`. It is invalid.

### 7.4 PubChem BioAssay

Find AIDs by gene ID:

```text
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/target/geneid/6721/aids/JSON
```

Assay description:

```text
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/615679/description/JSON
```

Assay data:

```text
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/615679/CSV
```

Fields:

```text
IdentifierList.AID
PC_AssayContainer[*].assay.descr.name
PC_AssayContainer[*].assay.descr.aid.id
PC_AssayContainer[*].assay.descr.comment
CSV: PUBCHEM_CID
CSV: PUBCHEM_ACTIVITY_OUTCOME
CSV: Standard Type
CSV: Standard Value
CSV: Standard Units
```



### 7.5 Open Targets Tractability

```text
POST https://api.platform.opentargets.org/api/v4/graphql
```

Body:

```json
{
  "query": "query Srebf2Tractability($ensemblId: String!) { target(ensemblId: $ensemblId) { id approvedSymbol approvedName tractability { modality label value } chemicalProbes { id drugId drugFromSourceId targetFromSourceId mechanismOfAction isHighQuality origin probeMinerScore probesDrugsScore control urls { niceName url } } } }",
  "variables": {
    "ensemblId": "ENSG00000198911"
  }
}
```

Fields:

```text
data.target.tractability[*].modality
data.target.tractability[*].label
data.target.tractability[*].value
data.target.chemicalProbes[*]
```

Is it what we need? Yes, for tractability, not for direct inhibitors.

NOTE: DrugBank requires a licensed API key. NCATS Inxight did not validate as useful for SREBF2 in this workflow.

---



## Section 8 - Brain eQTLs: GTEx

All tissues:

```text
GET https://gtexportal.org/api/v2/association/singleTissueEqtl?gencodeId=ENSG00000198911.11&datasetId=gtex_v8&itemsPerPage=100000
```

Brain tissues:

```text
GET https://gtexportal.org/api/v2/association/singleTissueEqtl?gencodeId=ENSG00000198911.11&datasetId=gtex_v8&tissueSiteDetailId=Brain_Amygdala&tissueSiteDetailId=Brain_Anterior_cingulate_cortex_BA24&tissueSiteDetailId=Brain_Caudate_basal_ganglia&tissueSiteDetailId=Brain_Cerebellar_Hemisphere&tissueSiteDetailId=Brain_Cerebellum&tissueSiteDetailId=Brain_Cortex&tissueSiteDetailId=Brain_Frontal_Cortex_BA9&tissueSiteDetailId=Brain_Hippocampus&tissueSiteDetailId=Brain_Hypothalamus&tissueSiteDetailId=Brain_Nucleus_accumbens_basal_ganglia&tissueSiteDetailId=Brain_Putamen_basal_ganglia&tissueSiteDetailId=Brain_Spinal_cord_cervical_c-1&tissueSiteDetailId=Brain_Substantia_nigra&itemsPerPage=100000
```

Fields:

```text
data[*].snpId
data[*].variantId
data[*].pos
data[*].chromosome
data[*].gencodeId
data[*].geneSymbol
data[*].tissueSiteDetailId
data[*].nes
data[*].pValue
data[*].datasetId
```

Validated counts:

```text
all tissues: 748
brain tissues: 42
```

NOTE: This endpoint does not return ref/alt allele, MAF, or slope. It returns GTEx significant eQTL association fields.

---



## Section 9 - SNPs / ClinVar / OMIM / Open Targets



### 9.1 ClinVar ESearch

```text
GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=SREBF2%5Bgene%5D%20AND%20single_gene%5Bprop%5D&retmode=json&retmax=500&api_key={{ncbi_api_key}}
```

Fields:

```text
esearchresult.count
esearchresult.idlist
```



### 9.2 ClinVar ESummary

```text
GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={{clinvar_ids}}&retmode=json&api_key={{ncbi_api_key}}
```

Fields:

```text
result.{uid}.accession
result.{uid}.title
result.{uid}.obj_type
result.{uid}.variation_set[0].variation_name
result.{uid}.variation_set[0].measure_id
result.{uid}.variation_set[0].cdna_change
result.{uid}.variation_set[0].canonical_spdi
result.{uid}.variation_set[0].variation_loc[*].assembly_name
result.{uid}.variation_set[0].variation_loc[*].chr
result.{uid}.variation_set[0].variation_loc[*].start
result.{uid}.variation_set[0].variation_loc[*].stop
result.{uid}.genes[0].symbol
result.{uid}.genes[0].geneid
result.{uid}.germline_classification.description
result.{uid}.germline_classification.review_status
result.{uid}.germline_classification.last_evaluated
result.{uid}.trait_set[*].trait_name
result.{uid}.molecular_consequence_list[*]
result.{uid}.protein_change
```

NOTE: Use comma-separated IDs with no quotes/newlines. Long URLs should be chunked.

### 9.3 ClinVar EFetch Optional

```text
GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&id=4836289,4589479&rettype=vcv&retmode=xml&api_key={{ncbi_api_key}}
```

NOTE: Your failed EFetch had an empty/malformed ID list. ESummary is enough for the dossier table.

### 9.4 OMIM Search

```text
GET https://api.omim.org/api/entry/search?search=SREBF2&include=geneMap&format=json&apiKey={{omim_api_key}}
```



### 9.5 OMIM Entry

```text
GET https://api.omim.org/api/entry?mimNumber=600481&include=geneMap,clinicalSynopsis,text&format=json&apiKey={{omim_api_key}}
```

Fields:

```text
omim.entryList[*].entry.mimNumber
omim.entryList[*].entry.titles.preferredTitle
omim.entryList[*].entry.titles.alternativeTitles
omim.entryList[*].entry.geneMap.chromosome
omim.entryList[*].entry.geneMap.cytoLocation
omim.entryList[*].entry.geneMap.computedCytoLocation
omim.entryList[*].entry.geneMap.geneSymbols
omim.entryList[*].entry.geneMap.geneName
omim.entryList[*].entry.geneMap.geneIDs
omim.entryList[*].entry.geneMap.ensemblIDs
omim.entryList[*].entry.geneMap.phenotypeMapList
```

NOTE: OMIM confirmed SREBF2 MIM `600481`, but no strong OMIM disease/phenotype relationship was confirmed from the validated response.

### 9.6 Open Targets Disease Associations

```text
POST https://api.platform.opentargets.org/api/v4/graphql
```

Body:

```json
{
  "query": "query Srebf2DiseaseEvidence($ensemblId: String!) { target(ensemblId: $ensemblId) { id approvedSymbol approvedName associatedDiseases(page: { index: 0, size: 1000 }) { count rows { score disease { id name } datatypeScores { id score } datasourceScores { id score } } } } }",
  "variables": {
    "ensemblId": "ENSG00000198911"
  }
}
```

Fields:

```text
data.target.associatedDiseases.count
data.target.associatedDiseases.rows[*].disease.id
data.target.associatedDiseases.rows[*].disease.name
data.target.associatedDiseases.rows[*].score
data.target.associatedDiseases.rows[*].datatypeScores[*].id
data.target.associatedDiseases.rows[*].datatypeScores[*].score
data.target.associatedDiseases.rows[*].datasourceScores[*].id
data.target.associatedDiseases.rows[*].datasourceScores[*].score
```

Validated SREBF2 count:

```text
950 disease associations
```

NOTE: Current Open Targets counts may differ from older sample reports.

---



## Section 10 - Major Pathways



### 10.1 Reactome Pathways

```text
GET https://reactome.org/ContentService/data/mapping/UniProt/Q12772/pathways
```

Fields:

```text
[*].dbId
[*].stId
[*].stIdVersion
[*].displayName
[*].speciesName
[*].doi
[*].hasDiagram
[*].releaseDate
[*].lastUpdatedDate
[*].schemaClass
```

Build links:

```text
https://reactome.org/content/detail/{stId}
https://reactome.org/PathwayBrowser/#/{stId}&FLG=Q12772
```



### 10.2 Reactome Pathway Detail

```text
GET https://reactome.org/ContentService/data/query/R-HSA-1655829
```

Fields:

```text
stId
displayName
summation[*].text
literatureReference[*].pubMedIdentifier
```



### 10.3 WikiPathways Text Bulk File

```text
GET https://www.wikipathways.org/json/findPathwaysByText.json
```



### 10.4 WikiPathways Xref Bulk File

```text
GET https://www.wikipathways.org/json/findPathwaysByXref.json
```

Fields:

```text
pathwayInfo[*].id
pathwayInfo[*].url
pathwayInfo[*].name
pathwayInfo[*].species
pathwayInfo[*].revision
pathwayInfo[*].description
pathwayInfo[*].ncbigene
pathwayInfo[*].uniprot
pathwayInfo[*].ensembl
```

Matched identifiers:

```text
6721
Q12772
ENSG00000198911
SREBF2
```

NOTE: WikiPathways deprecated the old live query web services. Current workflow is bulk JSON download plus local filtering.

---



## Section 11 - Knockouts And Mouse Phenotypes



### 11.1 PubMed Literature Search

```text
GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=%28SREBF2%20OR%20SREBP-2%20OR%20Srebf2%29%20AND%20%28knockout%20OR%20deficient%20OR%20deletion%20OR%20conditional%20knockout%20OR%20hypomorphic%20OR%20transgenic%20OR%20RNA%20interference%29%20AND%20%28mouse%20OR%20mice%20OR%20Mus%20musculus%29&retmode=json&retmax=50&sort=relevance&api_key={{ncbi_api_key}}
```

Fields:

```text
esearchresult.count
esearchresult.idlist
```



### 11.2 PubMed Summaries

```text
GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={{pubmed_ids}}&retmode=json&api_key={{ncbi_api_key}}
```

Fields:

```text
result.{pmid}.uid
result.{pmid}.title
result.{pmid}.pubdate
result.{pmid}.source
result.{pmid}.authors[*].name
```



### 11.3 MouseMine Gene Lookup

```text
GET https://www.mousemine.org/mousemine/service/query/results?format=json&query=%3Cquery%20model%3D%22genomic%22%20view%3D%22Gene.primaryIdentifier%20Gene.symbol%20Gene.name%20Gene.organism.name%20Gene.ncbiGeneNumber%22%20sortOrder%3D%22Gene.primaryIdentifier%20asc%22%3E%3Cconstraint%20path%3D%22Gene.ncbiGeneNumber%22%20op%3D%22%3D%22%20value%3D%2220788%22%2F%3E%3C%2Fquery%3E
```

Validated:

```text
MGI:107585
Srebf2
20788
```



### 11.4 MouseMine Alleles

```text
GET https://www.mousemine.org/mousemine/service/query/results?format=json&query=%3Cquery%20model%3D%22genomic%22%20view%3D%22Allele.primaryIdentifier%20Allele.symbol%20Allele.name%20Allele.alleleType%20Allele.feature.primaryIdentifier%20Allele.feature.symbol%22%20sortOrder%3D%22Allele.symbol%20asc%22%3E%3Cconstraint%20path%3D%22Allele.feature.primaryIdentifier%22%20op%3D%22%3D%22%20value%3D%22MGI%3A107585%22%2F%3E%3C%2Fquery%3E
```



### 11.5 MouseMine Allele Phenotypes

```text
GET https://www.mousemine.org/mousemine/service/query/results?format=json&query=%3Cquery%20model%3D%22genomic%22%20view%3D%22Allele.primaryIdentifier%20Allele.symbol%20Allele.name%20Allele.alleleType%20Allele.ontologyAnnotations.ontologyTerm.identifier%20Allele.ontologyAnnotations.ontologyTerm.name%22%20sortOrder%3D%22Allele.symbol%20asc%22%3E%3Cconstraint%20path%3D%22Allele.feature.primaryIdentifier%22%20op%3D%22%3D%22%20value%3D%22MGI%3A107585%22%2F%3E%3C%2Fquery%3E
```

Fields:

```text
Allele.primaryIdentifier
Allele.symbol
Allele.name
Allele.alleleType
Allele.ontologyAnnotations.ontologyTerm.identifier
Allele.ontologyAnnotations.ontologyTerm.name
```

Validated result:

```text
15 phenotype annotations across 3 alleles
```



### 11.6 MouseMine Stocks / Carried By

```text
GET https://www.mousemine.org/mousemine/service/query/results?format=json&query=%3Cquery%20model%3D%22genomic%22%20view%3D%22Allele.primaryIdentifier%20Allele.symbol%20Allele.name%20Allele.alleleType%20Allele.carriedBy.primaryIdentifier%20Allele.carriedBy.symbol%20Allele.carriedBy.name%22%20sortOrder%3D%22Allele.symbol%20asc%22%3E%3Cconstraint%20path%3D%22Allele.feature.primaryIdentifier%22%20op%3D%22%3D%22%20value%3D%22MGI%3A107585%22%2F%3E%3C%2Fquery%3E
```

Fields:

```text
Allele.primaryIdentifier
Allele.symbol
Allele.name
Allele.alleleType
Allele.carriedBy.primaryIdentifier
Allele.carriedBy.name
```

NOTE: `Gene.alleles.phenotypeAnnotations` is not in the MouseMine model. Correct path is `ontologyAnnotations`.

---



## Section 12 - Major Labs Based On Publications



### 12.1 PubMed Major Author Search

```text
GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=%28SREBF2%20OR%20SREBP-2%20OR%20SREBP2%20OR%20%22sterol%20regulatory%20element%20binding%20factor%202%22%29%20AND%20%28cholesterol%20OR%20lipid%20OR%20sterol%20OR%20SREBP%29&retmode=json&retmax=200&sort=relevance&api_key={{ncbi_api_key}}
```



### 12.2 PubMed Summary Metadata

```text
GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={{pubmed_ids}}&retmode=json&api_key={{ncbi_api_key}}
```

Fields:

```text
result.{pmid}.uid
result.{pmid}.title
result.{pmid}.pubdate
result.{pmid}.source
result.{pmid}.fulljournalname
result.{pmid}.authors[*].name
result.{pmid}.lastauthor
```

Author aggregation:

```text
count all author appearances
count first-author appearances
count last-author appearances
```



### 12.3 PubMed Affiliation XML

```text
GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id={{selected_pubmed_ids}}&retmode=xml&rettype=abstract&api_key={{ncbi_api_key}}
```

Fields:

```text
PubmedArticle.MedlineCitation.PMID
Article.ArticleTitle
Article.AuthorList.Author.LastName
Article.AuthorList.Author.ForeName
Article.AuthorList.Author.AffiliationInfo.Affiliation
```



### 12.4 OpenAlex Works

```text
GET https://api.openalex.org/works?search=SREBF2%20SREBP2%20SREBP-2&filter=from_publication_date:1990-01-01&per-page=200&select=id,doi,title,publication_year,cited_by_count,authorships,primary_location
```

Fields:

```text
results[*].id
results[*].doi
results[*].title
results[*].publication_year
results[*].cited_by_count
results[*].authorships[*].author.display_name
results[*].authorships[*].author.id
results[*].authorships[*].institutions[*].display_name
results[*].authorships[*].institutions[*].country_code
```

Is it what we need? Yes, as an author/institution aggregation aid.

NOTE: Lab websites and emails are not reliably available from PubMed/OpenAlex. Add manually from institution pages.

---



## Section 13 - Commercial Antibodies

Purpose: discover commercial SREBF2/SREBP2 antibodies.

### 13.1 SerpAPI Google Search

```text
GET https://serpapi.com/search.json?engine=google&q=SREBF2%20OR%20SREBP2%20antibody%20catalog%20Abcam%20Novus%20R%26D%20Santa%20Cruz%20OriGene&api_key={{serpapi_api_key}}
```

Targeted searches:

```text
GET https://serpapi.com/search.json?engine=google&q=%22SREBP2%20Antibody%22%20%22catalog%22%20%22SREBF2%22&api_key={{serpapi_api_key}}
GET https://serpapi.com/search.json?engine=google&q=%22SREBP2%20Antibody%22%20%28Abcam%20OR%20Novus%20OR%20R%26D%20Systems%20OR%20LSBio%20OR%20Biorbyt%20OR%20Sino%20Biological%20OR%20Proteintech%29&api_key={{serpapi_api_key}}
```

Fields:

```text
organic_results[*].title
organic_results[*].source
organic_results[*].link
organic_results[*].snippet
```

Derived fields:

```text
antibody_product_name
vendor_name
catalog_number
antibody_description
product_url
```

Is it what we need? Yes for discovery.

NOTE: Vendor product pages are final provenance. Google/SerpAPI snippets are discovery metadata and can be stale.

NOTE: There is no single reliable official REST API for commercial antibodies.

---



## Section 14 - Patents



### 14.1 SerpAPI Google Patents

Broad search:

```text
GET https://serpapi.com/search.json?engine=google_patents&q=%28SREBF2%20OR%20SREBP2%20OR%20%22SREBP-2%22%29&api_key={{serpapi_api_key}}
```

Targeted searches:

```text
GET https://serpapi.com/search.json?engine=google_patents&q=%22SREBF2%22&api_key={{serpapi_api_key}}
GET https://serpapi.com/search.json?engine=google_patents&q=%22SREBP2%22&api_key={{serpapi_api_key}}
GET https://serpapi.com/search.json?engine=google_patents&q=%22SREBP-2%22&api_key={{serpapi_api_key}}
GET https://serpapi.com/search.json?engine=google_patents&q=%22sterol%20regulatory%20element-binding%20protein%202%22&api_key={{serpapi_api_key}}
```

Fields:

```text
organic_results[*].title
organic_results[*].publication_number
organic_results[*].patent_id
organic_results[*].assignee
organic_results[*].inventor
organic_results[*].priority_date
organic_results[*].filing_date
organic_results[*].publication_date
organic_results[*].grant_date
organic_results[*].snippet
organic_results[*].patent_link
organic_results[*].pdf
```

Build link:

```text
https://patents.google.com/{patent_id}
```

Is it what we need? Yes, for patent table discovery and summaries.

NOTE: Mark each patent with a relevance level: `direct`, `pathway`, `marker-list`, or `weak`.

NOTE: Some hits only mention SREBF2 in large gene lists. Do not overinterpret those as SREBF2-specific patents.

---



## Section 15 - NIH And ERC Grants



### 15.1 NIH RePORTER Exact Search

```text
POST https://api.reporter.nih.gov/v2/projects/search
```

Headers:

```text
Content-Type: application/json
Accept: application/json
```

Body:

```json
{
  "criteria": {
    "advanced_text_search": {
      "operator": "or",
      "search_field": "projecttitle,abstracttext,terms",
      "search_text": "SREBF2 SREBP2 SREBP-2 Srebf2 \"sterol regulatory element binding protein 2\" \"sterol regulatory element-binding protein 2\""
    },
    "include_active_projects": true,
    "exclude_subprojects": true
  },
  "include_fields": [
    "ProjectTitle",
    "ProjectNum",
    "CoreProjectNum",
    "OrganizationName",
    "BudgetStartDate",
    "BudgetEndDate",
    "AgencyICAdmin",
    "AgencyICFundings",
    "FiscalYear",
    "AwardAmount",
    "PrincipalInvestigators",
    "ProjectDetailUrl",
    "AbstractText",
    "Terms"
  ],
  "offset": 0,
  "limit": 100,
  "sort_field": "project_start_date",
  "sort_order": "desc"
}
```



### 15.2 NIH RePORTER Broader Pathway Search

```text
POST https://api.reporter.nih.gov/v2/projects/search
```

Body:

```json
{
  "criteria": {
    "advanced_text_search": {
      "operator": "or",
      "search_field": "projecttitle,abstracttext,terms",
      "search_text": "SREBP cholesterol mevalonate sterol lipid metabolism cholesterol biosynthesis"
    },
    "include_active_projects": true,
    "exclude_subprojects": true
  },
  "include_fields": [
    "ProjectTitle",
    "ProjectNum",
    "CoreProjectNum",
    "OrganizationName",
    "BudgetStartDate",
    "BudgetEndDate",
    "AgencyICAdmin",
    "AgencyICFundings",
    "FiscalYear",
    "AwardAmount",
    "PrincipalInvestigators",
    "ProjectDetailUrl",
    "AbstractText",
    "Terms"
  ],
  "offset": 0,
  "limit": 100,
  "sort_field": "project_start_date",
  "sort_order": "desc"
}
```

Fields:

```text
results[*].project_title
results[*].project_num
results[*].core_project_num
results[*].organization.org_name
results[*].budget_start
results[*].budget_end
results[*].agency_ic_admin.name
results[*].agency_ic_fundings[*].name
results[*].fiscal_year
results[*].award_amount
results[*].principal_investigators[*].full_name
results[*].project_detail_url
results[*].abstract_text
results[*].terms
```

Is it what we need? Yes, for NIH grant table.

NOTE: Use exact search first. Use broader pathway search only if the project title/abstract/terms clearly connect to SREBF2/SREBP biology.

### 15.3 ERC / European Data Portal

Attempted endpoint family:

```text
https://data.europa.eu/api/hub/search/
```

Search terms attempted:

```text
SREBF2
SREBP2
SREBP-2
"sterol regulatory element binding protein 2"
```

Status:

```text
No usable gene-level ERC grant records returned during manual validation.
```

Recommended report wording:

```text
Searches for SREBF2, SREBP2, SREBP-2, and "sterol regulatory element binding protein 2" through the European Data Portal route did not identify a validated SREBF2-related ERC grant.
```

NOTE: European Data Portal is a metadata catalog, not a direct gene-level ERC grant API. Do not treat generic catalog results as ERC grant hits.

---



## Recommended Normalized Output Tables



### General Gene Info


| Species | Entrez Gene ID | Gene Symbol | Gene Name | Ensembl ID | UniProt ID | Aliases |
| ------- | -------------- | ----------- | --------- | ---------- | ---------- | ------- |




### Pathways / PPI / Patents / Grants

Use one row per source record with:

```text
source
source_url
source_access_date
gene_symbol
gene_id
record_id
record_title_or_name
record_description
matched_field
matched_value
relevance_level
raw_response_pointer
```



### Provenance Rule

Every field in the dossier should trace back to:

```text
API endpoint
request parameters
response JSON/TSV/XML path
date accessed
```

NOTE: Store raw response files or response hashes later when moving to Python.