# Genomics Security

## Overall Design

A software tool that will:
- Take FastA input of a reference and subject sequence
- Perform local alignment of the subject to the reference
    - Take the first match found
    - Perform two analyses:
        - Track changes between purine/pyrimidine (NOT DONE!)
        - Track changes between Codons within a Protein

## Not Covered By This Design

- Bottom-up approach where entire subject gene sequence is utilized to identify overall structural and regional differences
- Instances in which gene duplicates found in alignment are analyzed independently (maybe we can offer just count?)

## File Size Estimation

Two factors would go into estimating the scale of data stored in the resultant format:
- Length of genome/gene/protein as it pertains to nucleotides
    - Two bits stored per nucleotide in order -- whether purine/pyrimidine and whether it has mutated from the reference
- Length of genome/gene/protein as it pertains to codons
    - Two bits for hyrophobicity, aromaticity, specialty -- true/false and whether mutated
    - Three bits for charge and whether mutated, as charge has effectively three states: positive, negative, none
    - ...so nine bits total for codon mutations
- Aside from these there is sequence metadata and potentially subject demographic metadata

Given the above, we can estimate file size roughly on magnitudes given size of sequence comparison saved (assuming ref/sample are similarly sized! this is optimism!):
- 10k chars in ref/sample = roughly 2.5kB for nucleotide mutations and roughly 3.7kB for codon mutations = roughly 6.2kB total (nucleotide count * 4 / 8; codon count / 3 * 9 / 8)
- 100k chars in ref/sample = roughly 25kB for nucleotides and roughly 37kB for codons
- 1billion chars in ref/sample = roughly 250mB for nucleotides and roughly 370mB for codons
- Human genome at 3.2billion = roughly 2gB of storage required per comparison against reference stored
