# Genomics Security

## Overall Design

A software tool that will:
- Take FastA input of a reference and subject sequence
- Perform local alignment of the subject to the reference
    - Take the first match found
    - Perform two analyses:
        - Track changes between purine/pyrimidine
        - Track changes between Codons within a Protein

## Not Covered By This Design

- Bottom-up approach where entire subject gene sequence is utilized to identify overall structural and regional differences
- Instances in which gene duplicates found in alignment are analyzed independently (maybe we can offer just count?)


## Need to:
- Find tata boxes
- Grab everything before the first by 500 characters and after the last by 500 characters

